#! /usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import random

try:
    import anvio.utils as u
except ModuleNotFoundError:
    print("Sorry :/ This script relies on some anvi'o modules, and therefore you need anvi'o to be installed on your system.")
    sys.exit(-1)
import anvio.fastalib as fastalib
import anvio.terminal as terminal

import reads_for_assembly.configs as configs

from reads_for_assembly.utils import simulate_errors

pp = terminal.pretty_print


class SingleReads(configs.SingleReadsConfiguration):
    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        configs.SingleReadsConfiguration.__init__(self, args)


    def generate(self):
        output = fastalib.FastaOutput(self.output_file)

        for index_fasta in range(0, len(self.fasta_files)):
            f = self.fasta_files_dict[self.fasta_files[index_fasta]]

            x = self.short_read_length
            c = f['coverage']

            self.progress.new('Working on file %d of %d (%s) with expected coverage of %d' % (index_fasta + 1, len(self.fasta_files), f['alias'], c))

            fasta = fastalib.SequenceSource(f['path'])
            total_num_errors = 0
            total_num_reads = 0
            while next(fasta):
                L = len(fasta.seq)

                av_num_short_reads_needed = int(L / x * c)
                total_num_reads += av_num_short_reads_needed

                for index_short_read in range(0, av_num_short_reads_needed):
                    if (index_short_read + 1) % 100 == 0:
                        self.progress.update('Entry %s :: %s nts :: reads %s of %s :: num errors: %s ...'\
                                                        % (pp(fasta.pos + 1), pp(len(fasta.seq)),
                                                           pp(index_short_read + 1), pp(av_num_short_reads_needed),
                                                           pp(total_num_errors)))

                    start_pos = random.randint(0, L - x)
                    short_read, num_errors = simulate_errors(self.error_rate, fasta.seq[start_pos:start_pos + x])
                    total_num_errors += num_errors

                    output.write_id('%s_%d|source:%s|start:%d|stop:%d' % (f['alias'], index_short_read, fasta.id, start_pos, start_pos + x))
                    output.write_seq(short_read)

            self.progress.end()
            self.run.info('%s w/ %d contigs' % (f['alias'], fasta.pos), '%s reads with %s errors (avg %.4f) for %sX avg cov.'\
                                        % (pp(total_num_reads),
                                           pp(total_num_errors),
                                           total_num_errors * 1.0 / (total_num_reads * x),
                                           pp(c),
                                           ))

        output.close()
        self.run.info('Fasta output', self.output_file)


class PairedEndReads(configs.PairedEndReadsConfiguration):
    def __init__(self, args, r=terminal.Run(), p=terminal.Progress()):
        self.args = args
        self.run = r
        self.progress = p

        configs.PairedEndReadsConfiguration.__init__(self, args)


    def generate(self):
        output_r1 = open(self.output_sample_name + '-R1.fastq', 'w')
        output_r2 = open(self.output_sample_name + '-R2.fastq', 'w')

        min_sequence_length = (self.short_read_length * 2) + (self.insert_size * 3)

        self.run.info('Read lenth', self.short_read_length)
        self.run.info('Insert size', self.insert_size)
        self.run.info('Insert size std', self.insert_size_std)
        self.run.info('Minimum sequence length REQUIRED', min_sequence_length)

        x = self.short_read_length
        self.Q_str = ''.join(['A'] * x)

        # FIRST WE SANITY CHECK FASTA FILES
        bad_fasta_files = []
        self.progress.new("Sanity checking FASTA files")
        for index_fasta in range(0, len(self.fasta_files)):
            self.progress.update('%d of %d' % (index_fasta + 1, len(self.fasta_files)))
            f = self.fasta_files_dict[self.fasta_files[index_fasta]]
            fasta = fastalib.SequenceSource(f['path'])

            bad_lengths = []
            while next(fasta):
                if len(fasta.seq) < min_sequence_length:
                    bad_lengths.append(len(fasta.seq))

            if len(bad_lengths):
                bad_fasta_files.append((f['alias'], len(bad_lengths)), )

        if len(bad_fasta_files):
            self.progress.reset()
            self.run.warning("Some of your FASTA files contain sequences that are shorter than the minimum seqeunce length "
                             "required given the short read length and insert size you have set in your config file. Briefly, each "
                             "sequence in your reference FASTA files must be longer than %d base pair. But the following "
                             "FASTA files have reads shorter than this value:" % (int(min_sequence_length)), header="PLEASE READ CAREFULLY")
            for alias, num_bad_reads in bad_fasta_files:
                self.run.info_single("%s has %d reads shorter than %d" % (alias, num_bad_reads, int(min_sequence_length)), mc="red")

            print()
            sys.exit(-1)
        self.progress.end()

        for index_fasta in range(0, len(self.fasta_files)):
            f = self.fasta_files_dict[self.fasta_files[index_fasta]]

            c = f['coverage']

            self.progress.new('Working on file %d of %d (%s) with expected coverage of %d' % (index_fasta + 1, len(self.fasta_files), f['alias'], c))

            fasta = fastalib.SequenceSource(f['path'])
            total_r1_errors = 0
            total_r2_errors = 0
            total_num_reads = 0
            while next(fasta):
                L = len(fasta.seq)

                av_num_short_reads_needed = int(L / x * c)
                total_num_reads += av_num_short_reads_needed

                av_num_pairs_needed = int(av_num_short_reads_needed / 2)

                for index_pair in range(0, av_num_pairs_needed):
                    if (index_pair + 1) % 100 == 0:
                        self.progress.update('Seq %s :: %s nts :: reads %s of %s :: num errors: %s ...'\
                                                        % (pp(fasta.pos + 1), pp(len(fasta.seq)),
                                                           pp(index_pair + 1), pp(av_num_pairs_needed),
                                                           pp(total_r1_errors + total_r2_errors)))

                    I = int(round(random.gauss(self.insert_size, self.insert_size_std)))

                    range_max = L - ((x * 2) + I)

                    start_pos = random.randint(0, range_max)

                    read_1_start = start_pos
                    read_1_stop = read_1_start + x

                    read_2_start = read_1_stop + I
                    read_2_stop = read_2_start + x

                    read_1, num_errors_r1 = simulate_errors(self.error_rate, fasta.seq[read_1_start:read_1_stop], bases=['A', 'T', 'C', 'G'])
                    read_2, num_errors_r2 = simulate_errors(self.error_rate, fasta.seq[read_2_start:read_2_stop], bases=['A', 'T', 'C', 'G'])

                    total_r1_errors += num_errors_r1
                    total_r2_errors += num_errors_r2

                    c1, c2 = random.randint(1, 10000), random.randint(1, 10000)
                    output_r1.write('@%s:23:B02CBACXX:8:2315:%d:%d 1:N:0:GATCAG\n' % (f['alias'], c1, c2))
                    output_r1.write(read_1 + '\n')
                    output_r1.write('+source:%s; start:%d; stop:%d; insert_size:%d; num_errors:%d\n' % (fasta.id, read_1_start, read_1_stop, I, num_errors_r1))
                    output_r1.write('%s\n' % self.Q_str)

                    output_r2.write('@%s:23:B02CBACXX:8:2315:%d:%d 2:N:0:GATCAG\n' % (f['alias'], c1, c2))
                    output_r2.write(u.rev_comp(read_2) + '\n')
                    output_r2.write('+source:%s; start:%d; stop:%d; insert_size:%d; num_errors:%d\n' % (fasta.id, read_2_start, read_2_stop, I, num_errors_r2))
                    output_r2.write('%s\n' % self.Q_str)

            self.progress.end()
            total_num_errors = total_r1_errors + total_r2_errors
            self.run.info('%s w/ %d contigs' % (f['alias'], fasta.pos),
                     '%s reads in %s pairs with %s errors (avg %.4f) for %sX avg cov.'\
                                        % (pp(total_num_reads),
                                           pp(total_num_reads / 2),
                                           pp(total_num_errors),
                                           total_num_errors * 1.0 / (total_num_reads * x) if total_num_errors else 0,
                                           pp(c),
                                           ))

        output_r1.close()
        output_r2.close()
        self.run.info('FASTQ R1',self.output_sample_name + '-R1.fastq')
        self.run.info('FASTQ R2',self.output_sample_name + '-R2.fastq')



