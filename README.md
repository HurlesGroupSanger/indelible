# Indelible
#### Genomic Structural Variant Caller by Adaptive Training

## Authors

Alejandro Sifrim (Creator, Developer)

Diana Rajan (Experimental validation)

Elena Prigmore (Experimental validation)

Matthew Hurles (Group Leader)

We are affiliated with the [Wellcome Trust Sanger Institute](http://www.sanger.ac.uk/science/groups/hurles-group), Cambridge, United Kingdom

## Abstract

TBD

## Installation

## Configuration

## Usage

The main help page of the program can be accessed by executing the indelible script with the `-h` flag as follows:

```
./indelible.py -h
usage: indelible [-h] <command> ...

positional arguments:
  <command>   One of the following commands:
    fetch     fetch reads from BAM file
    aggregate
              aggregate information per position
    score     score positions using Random Forest model
    blast     blast clipped sequences
    annotate  annotate positions with additional information
    denovo    searches for de novo events
    train     trains the Random Forest model on a bunch of examples
    complete  Performs the complete Indelible analysis

optional arguments:
  -h, --help  show this help message and exit
```

The Indelible variant calling process follows several steps:

1. Soft-clipped reads are extracted from the BAM files
2. Information is aggregated across reads to find positions where multiple reads are clipped.
3. Positions are scored using a Random Forest model taking into account the number/quality of clipped reads and the sequence context
4. Longer clipped segments are blasted against the human genome and repeat databases
4. Putative SVs are annotated with additional information (e.g. gene annotations) and the blast results from the previous step
5. If trio information is available, de novo events are called and inheritance information is appended.

These different steps can be performed by the different sub-commands:

The **fetch** command extracts the reads from the BAM file, it takes 2 arguments:

* `--i` : path to the input BAM file.
* `--o` : path to output the clipped reads to.

```
./indelible.py fetch --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.sc_reads
```

The **aggregate** merges information across reads towards a position-level view of the data:

* `--i` : path to the input file (the output of the *fetch* command from previous step).
* `--b` : path to the BAM file used to generate the input file.
* `--o` : the path to the output file.
* `--r` : path to reference genome 

```
./indelible.py aggregate --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.sc_reads --b test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts --r /lustre/scratch113/projects/ddd/resources/v1.2/hs37d5.fasta
```

The **score** command scores positions based on the read information and sequence context:

* `--i` : path to the input file (the output of the *aggregate* command from previous step).
* `--o` : the path to the output file.

```
./indelible.py score --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored
```

The **blast** command blasts longer clipped segments (20+ bp) to find matches elsewhere in the human genome and/or repeat database:
* `--i` :  path to the input file (the output of the *score* command from previous step).

This will automatically generate 3 files: a fasta file with the sequences to be blasted, 2 results files (one for repeat sequences, one for non-repeat sequences). These files are used in the next step.

```
./indelible.py blast --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored
```

The **annotate** command enriches the result with gene/exon annotations and merges the blast results with the position file:

* `--i` : path to the input file (output of score command after running the blast command).
* `--o` : path to output the annotated file.

```
./indelible.py annotate --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored.annotated
```

All the previous commands can be performed in succession automatically with the **complete** command:

* `--i` : path to the input BAM file.
* `--o` : path to directory where output files will be generated.
* `--r` : path to reference genome
* `--keeptmp` : If this flag is given, intermediate files are kept. Otherwise these files will be removed once the analysis is finished.

```
./indelible.py complete --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o test_data/ --r /lustre/scratch113/projects/ddd/resources/v1.2/hs37d5.fasta --keeptmp
```

One can also look for *de novo* mutation events using the **denovo** command:

***Note***: If maternal and/or paternal bam files are not supplied, *denovo* filtering will not be performed. This behaviour is intended to format non-trio data identically to trio data. If one of maternal **or** paternal bam is provided, Indelible will count coverage within that bam.

* `--c` : path to scored/annotated calls in the proband.
* `--m` : path to maternal BAM file.
* `--p` : path to paternal BAM file.
* `--o` : path to output file.

```
./indelible.py denovo --c test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.indelible.tsv --m maternal.bam --p paternal.bam --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.indelible.denovo.tsv
```

