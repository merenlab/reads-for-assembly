# Short Reads for Assembly

This is a script to generate randomly-spliced short reads from a given set of
contigs, so they can be assembled back (to test assemblers, or genome binning
software, or more).

While doing this you can specify a desired __average coverage__ for each entry,
and define an expected __error rate__.

# An example

_Note: I am running all these commands from within the source code directory_

---

Say, you have multiple FASTA files with contigs from -draft- genomes, and you
wish to create a single FASTA file that contains enough number of short reads
randomly created from these contigs to meet your expected coverage when they are
assembled.

For a little demonstration I put a `files` directory with 5 FASTAS files:


    ls files/*fa
    files/fasta_01.fa files/fasta_02.fa files/fasta_03.fa files/fasta_04.fa files/fasta_05.fa

Each FASTA file is sampled from a different bacterial genome, contains 5 contigs
that are about 7,500nts long.

To generate short reads from this I have this config file:


    $ cat sample-config.ini
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
__sample/fasta_01.fa__, so when I assemble them back, the average coverage would
be about __100X__ throughout the contig. Oh, while doing this, introduce some random
base errors to meet __0.05__ error rate in average. In fact do the same for every
other entry in this config file with respect to their coverage values, and store
all these results in __short_reads.fa__".

Once you have your config file ready, this is how you run it:

    $ ./gen-short-reads sample-config.ini
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

# Need more from this tool?

It is indeed going to be useful for my purposes, but I would like it to be
useful to you as well. So if this is what you need to do benchmarks or to test
your stuff, and if its missing something for your purposes, please don't
hesitate to get in touch with [me](http://meren.org).
