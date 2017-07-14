# Simulating Short Reads for Assembly, Mapping, and Binning experiments

Scripts to generate random short reads from a given set of contigs, so they can be assembled back (to test assemblers, or genome binning software, or more).

Using the config files and scripts,

* You can generate **single reads** in a FASTA-formatted output file.
* You can generate **paired-end reads** in FASTQ-formatted output files. You can define an **insert size with a standard deviation**.
* You can generate reads from a single FASTA file with a single contig, single FASTA file with multiple contigs, or you can **simulate a complex metagenomic mixture using multiple FASTA files each of with with multiple contigs**.
* You can specify a desired __average coverage__ for each entry, and define an expected __error rate__ for the sequencer.

The tool (both with respect to its functionality, and design) can be improved dramatically, but we only pushed it as much as what we needed from it. If you need something that is not here, please feel free to open an issue. 


# Generating single short reads in a FASTA file

_Note: I am running all these commands from within the source code directory_

---

For this, you will use the program `gen-single-reads` with a config file that looks like `gen-single-reads-example.ini`.

Say, you have multiple FASTA files with contigs from multiple genomes, and you
wish to create a single FASTA file that contains enough number of short reads
randomly created from these contigs to meet your expected coverage when they are
assembled.

For a little demonstration there is a `files` directory in the repository that includes 5 FASTA files:

    $ ls files/*fa
    files/fasta_01.fa files/fasta_02.fa files/fasta_03.fa files/fasta_04.fa files/fasta_05.fa

Each FASTA file is sampled from a different bacterial genome, and contains 5 contigs
that are about 7,500nts long.

There is also a sample configuration file that describes how to work on these FASTA files:

    $ cat gen-single-reads-example.ini
    [general]
    output_file = short_reads.fa
    short_read_length = 100
    error_rate = 0.05
    
    [files/fasta_01.fa]
    coverage = 100
    
    [files/fasta_02.fa]
    coverage = 200
    
    [files/fasta_03.fa]
    coverage = 50
    
    [files/fasta_04.fa]
    coverage = 150
    
    [files/fasta_05.fa]
    coverage = 50

The config file simply says;

> "Generage enough __100nt__ long short reads from all entries you find in
__sample/fasta_01.fa__, so when they are assembled and short reads mapped back
to resulting contigs, the average coverage would be about __100X__ for each contig.
While doing this, introduce some random base errors to meet __0.05__ error rate in
average. In fact do the same for every other entry in this config file with respect
to their required coverage values, and store all resulting reads in __short_reads.fa__".

Once you have your config file ready, this is how you run it:

    $ ./gen-single-reads gen-single-reads-example.ini
    fasta_01 .......: 38,000 reads w/ 189,975 errors (average rate of 0.0500) generated for 100X average coverage.
    fasta_02 .......: 76,000 reads w/ 379,624 errors (average rate of 0.0500) generated for 200X average coverage.
    fasta_03 .......: 19,000 reads w/ 94,498 errors (average rate of 0.0497) generated for 50X average coverage.
    fasta_04 .......: 57,000 reads w/ 284,756 errors (average rate of 0.0500) generated for 150X average coverage.
    fasta_05 .......: 19,000 reads w/ 95,257 errors (average rate of 0.0501) generated for 50X average coverage.
    Fasta output ...: short_reads.fa

There it is. So in theory, if you assemble short_reads.fa, you should recover
contigs that are about 50, 100, 150 and 200X coverage. Just for fun, I used
[velvet](https://www.ebi.ac.uk/~zerbino/velvet/) to assemble what is in the
short_reads.fa, then used [bowtie2](http://bowtie-
bio.sourceforge.net/bowtie2/index.shtml) to map reads back to resulting contigs,
and analyzed the mapping.

Here what I did step by step (excuse me for the crudeness of my demo and
directly copy-pasting from my command line):

    $ mkdir velvet_output/
    $ velveth velvet_output/ 17 short_reads.fa
    $ velvetg velvet_output/
    $ bowtie2-build velvet_output/contigs.fa contigs_ref
    $ bowtie2 -f short_reads.fa -x contigs_ref -S output.sam
    209000 reads; of these:
      209000 (100.00%) were unpaired; of these:
        4070 (1.95%) aligned 0 times
        204930 (98.05%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
    98.05% overall alignment rate
    $ samtools view -bS output.sam > output.bam

I further analyzed the BAM file in an _ad hoc_ manner with in-house scripts.
Here is the average coverage of each contig in `output.bam` file:

![average_coverage](https://raw.githubusercontent.com/meren/reads-for-assembly/master/files/average_coverage.png)

So randomly generated and assembled short reads did creaete expected coverage
profiles.

Here is the detailed coverage of one of those contigs that is at 50X coverage:

![average_coverage](https://raw.githubusercontent.com/meren/reads-for-assembly/master/files/50X.png)

And another example from a contig that is covered about 200X:

![average_coverage](https://raw.githubusercontent.com/meren/reads-for-assembly/master/files/200X.png)

# Generating paired-end reads in FASTQ R1/R2 files.

Everything is the same, except this time you will use the program `gen-paired-end-reads` with a config file that looks like `gen-paired-end-reads-example.ini`.

This is an example `gen-paired-end-reads-example.ini`:

``` ini
[general]
output_sample_name = test_sample
insert_size = 30
insert_size_std = 1
short_read_length = 100
error_rate = 0.05

[files/fasta_01.fa]
coverage = 10

[files/fasta_02.fa]
coverage = 20

[files/fasta_03.fa]
coverage = 5

[files/fasta_04.fa]
coverage = 15

[files/fasta_05.fa]
coverage = 5
```

`insert_size` is different than the real insert size, and describes the gap between the end of read 1 and the beginning of read 2.

Here is an example run:

``` bash
$ ./gen-paired-end-reads gen-paired-end-reads-example.ini
Read lenth ...............: 100
Insert size ..............: 30
Insert size std ..........: 1.0
fasta_01 .................: 380 paired-end reads w/ 18,845 errors (average rate of 0.2480) generated for 10X average coverage.
fasta_02 .................: 760 paired-end reads w/ 37,950 errors (average rate of 0.2497) generated for 20X average coverage.
fasta_03 .................: 190 paired-end reads w/ 9,635 errors (average rate of 0.2536) generated for 5X average coverage.
fasta_04 .................: 570 paired-end reads w/ 28,288 errors (average rate of 0.2481) generated for 15X average coverage.
fasta_05 .................: 190 paired-end reads w/ 9,665 errors (average rate of 0.2543) generated for 5X average coverage.
FASTQ R1 .................: test_sample-R1.fastq
FASTQ R2 .................: test_sample-R2.fastq
```

And this is how these files look like (those `A`'s are made up quality scores :/ it can easily be improved if you need more complex Q-score simulations):

``` bash
$ head -n 8 test_sample-R*
==> test_sample-R1.fastq <==
@fasta_01:23:B02CBACXX:8:2315:9436:8855 1:N:0:GATCAG
AAAGTGAGCGATACAAAGCTAAAGCCATCAACCAAAAAGCCCGCAACCCGCAAACCAGCTCCGCANACATCGGAGCTATCCAAGGACGATGACTNCCTGT
+source:640612206 slice:0-7600; start:5360; stop:5460; insert_size:32
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@fasta_01:23:B02CBACXX:8:2315:7851:9246 1:N:0:GATCAG
TGATGAGGCCTTCGTGGCCCAAGAGCGGACCAGCCGTAACTACGACCTGGTCGCCGTCTTTTAGGGCCTCGCTCATGGGCACCACGCGGTTTCCCCTATT
+source:640612206 slice:0-7600; start:2274; stop:2374; insert_size:31
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA

==> test_sample-R2.fastq <==
@fasta_01:23:B02CBACXX:8:2315:9436:8855 2:N:0:GATCAG
TAAATTCACCAATGACATGGACACTCTTTACGTCGGGCGCCCAAACGCAAATGCGCCAGCCGCGCGTGCGCTGCTGCGTGTCGGGATTCGCGCCGATCTT
+source:640612206 slice:0-7600; start:5492; stop:5592; insert_size:32
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
@fasta_01:23:B02CBACXX:8:2315:7851:9246 2:N:0:GATCAG
GCGATACGGCAGATCCTCGGACCGTGCAACAGTACCTGCTGCGCATGGATGACTTTGCCCGGGTGCTTTCACAGGACGGTCAGTTTGTGCCGTTGGCAAA
+source:640612206 slice:0-7600; start:2405; stop:2505; insert_size:31
AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
```


# Need more from this tool?

It is indeed going to be useful for my purposes, but I would like it to be
useful to you as well. So if this is what you need to do benchmarks or to test
your stuff, and if its missing something for your purposes, please don't
hesitate to get in touch with [me](http://meren.org).
