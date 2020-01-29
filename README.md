# InDelible
#### Genomic Structural Variant Caller by Adaptive Training

## Authors

Alejandro Sifrim (Creator, Developer)

Eugene Gardner (Developer)

Diana Rajan (Experimental validation)

Elena Prigmore (Experimental validation)

Sarah Lindsay (Experimental validation)

Matthew Hurles (Group Leader)

We are affiliated with the [Wellcome Sanger Institute](http://www.sanger.ac.uk/science/groups/hurles-group), Cambridge, United Kingdom

## About InDelible

### Abstract

TBD

### What is InDelible for?

### What is InDelible not for?

## Installation

### Required Software Dependencies

InDelible is written for Python2.7.* or Python3.7.*

InDelible requires the following software to be installed and in `$PATH`:

* [bedtools](https://bedtools.readthedocs.io/en/latest/)
  * **Note**: If using CRAM formated files with InDelible, bedtools v2.28 or later is required.
* [tabix](http://www.htslib.org/download/) 
* [bgzip](http://www.htslib.org/download/)

The python package Biopython requires a local install of [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) in `$PATH` in order to function. This needs to be installed prior to [installing InDelible](#installing-InDelible).

### Installing InDelible

To install InDelible:

1. Clone the git repo:

```
git clone https://github.com/eugenegardner/Indelible.git
cd Indelible/
```

2. Create a virtual environment and activate it:

* for Python2:
```
virtualenv venv/
source venv/bin/activate
```

* for Python3:
```
python3 -m venv venv/
source venv/bin/activate
```

3. Install cython and other required packages:

```
pip install cython
pip install -r requirements.txt
```

**Note**: If you get error(s) about pysam not being able to load specific libraries (like openssl, libbz2, etc.) that is a pysam problem related to htslib. Please see the pysam website to get help.

Indelible was tested with the following version of the packages in requirements.txt:

- cython v0.29.13
- numpy v1.17.2
- pandas v0.25.1
- pybedtools v0.8.0
- pyfaidx v0.5.5.2
- pysam v0.15.3
- scipy v1.3.1
- scikit-learn v0.21.3
  - **Note**: The random forest model provided in `data.zip` will likely only work with v0.21.3 of scikit-learn (see [this](https://scikit-learn.org/stable/modules/model_persistence.html#security-maintainability-limitations) link for an explanation why). Thus, if using a different version of scikit-learn, it is necessary to re-train the model with the provided test set included with this repository (`data/observation_data.DDD.17IX2019.txt`). Please see documentation [below](#_train) on how to train the random forest used by InDelible.
- PyYAML v5.1.2
- Biopython v1.74
- intervaltree v3.0.2

4. Unzip required data files:

```
cd data/
unzip data.zip
```


5. Download required blast resources:

```
## Download windowmasker:
wget ftp://ftp.ncbi.nlm.nih.gov/blast/windowmasker_files/9606/wmasker.obinary

## Download the GRCh37 human reference and create the blast db:
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai
makeblastdb -in hs37d5.fa -dbtype nucl
makeblastdb -in repeats.fasta -dbtype nucl
```

6. Edit the config file to point to required data files and edit any default parameters:

An example config.yml file (`example_config.yml`) is included in the top level directory. Edit this and save as `config.yml`.

```
cd ..

## Provide the path to each required file at the top of the config.yml
## Edit parameters at bottom of file.
vim example_config.yml

## cp/mv to config.yml
cp example_config.yml config.yml
```

Parameters in `config.yml` are as follows:

```
HC_THRESHOLD: Quality threshold at which to hard clip from the ends
MINIMUM_LENGTH_SPLIT_READ: Minimum number of bases that need to be soft clipped for a read to be included
MINIMUM_MAPQ: Minimum mapping quality for read to be included
MININUM_AVERAGE_BASE_QUALITY_SR: Minimum average base quality for the soft-clipped segment for read to be included
SHORT_SR_CUTOFF: The cutoff at which a soft-clipped segment is considered "short".
MINIMUM_SR_COVERAGE: Minimum number of reads with soft-clipped segments at a position for that position to be outputted.
SCORE_THRESHOLD: Minimum InDelible score to be considered for denovo calling
SR_THRESHOLD: Maximum number of clipped reads in parental samples to be considered inherited
COV_THRESHOLD: Minimum parental coverage to be able to call event as denovo
WINDOW_SIZE: window around position to look for indels/clipped reads (window_size/2 to the left and to the right)
```

**Note**: The above data resources will only work if you run InDelible using the the human GRCh37 reference. 

### Using vr-runner

InDelible comes packages with a script for use with the vr-runner packaged which automates the analysis of multiple 
samples. For information on how to install and use vr-runner, please see [this link](https://github.com/VertebrateResequencing/vr-runner).

The InDelible-specific vr-runner script is located at `./vr_runner/run-indelible`. This script requires a list of tab-delimited bam paths, where columns are:

1. child bam
2. mum bam
3. dad bam 

For example:

```
/path/to/child.bam  /path/to/mum.bam    /path/to/dad.bam
```

If mum or dad bams are not available, substitute a "-" like:

```
/path/to/child.bam  -   -
/path/to.child.bam  /path/to/mum.bam    -
/path/to.child.bam  - /path/to/dad.bam
```

A simple command for running InDelible with vr-runner is as follows:

```
./Indelible/vr_runner/run-indelible -o ./output/ -b crams.txt +maxjobs 1000 +loop 100 +retries -2
```

Resulting output will be in the `./output/` folder.

Additional commands for run-indelible:

```
   -b, --bams-list <file>      File with bam files
   -o, --outdir <dir>          Output directory
   -d, --database              Path to a precomputed MAF database [null]
```

***Note:*** If a filepath to a prebuilt MAF database (e.g. `./Indelible/data/`) is not provided for `-d/--database`, InDelible will rebuild the MAF database based only on samples included in `-b`!

## Usage

The main help page of the program can be accessed by executing the InDelible script with the `-h` flag as follows:

```
./indelible.py -h
usage: indelible [-h] <command> ...

positional arguments:
  <command>   One of the following commands:
    fetch     fetch reads from BAM file
    aggregate aggregate information per position
    score     score positions using Random Forest model
    blast     blast clipped sequences
    annotate  annotate positions with additional information
    denovo    searches for de novo events
    complete  Performs the complete InDelible analysis
    database  build SR allele frequency database
    train     trains the Random Forest model on a bunch of examples

optional arguments:
  -h, --help  show this help message and exit
```

### Primary SV Calling Pipeline

The InDelible variant calling process follows several steps:

1. [Fetch](#1-fetch) – Soft-clipped reads are extracted from the BAM files 
2. [Aggregate](#2-aggregate) - Information is aggregated across reads to find positions where multiple reads are clipped.
3. [Score](#3-score) - Positions are scored using a Random Forest model taking into account the number/quality of clipped reads and the sequence context
4. [Blast](#4-blast) - Longer clipped segments are blasted against the human genome and repeat databases
5. [Annotate](#5-annotate) - Putative SVs are annotated with additional information (e.g. gene annotations) and the blast results from the previous step
6. [_denovo_](#6-denovo) - _de novo_ events are called and inheritance information is appended.

These different steps can be performed by the different sub-commands:

#### 1. Fetch

The **fetch** command extracts the reads from the BAM file, it takes 2 arguments:

* `--i` : path to the input CRAM/BAM file.
* `--o` : path to output the clipped reads to.

```
./indelible.py fetch --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.sc_reads
```

#### 2. Aggregate

The **aggregate** merges information across reads towards a position-level view of the data:

* `--i` : path to the input file (the output of the *fetch* command from previous step).
* `--b` : path to the CRAM/BAM file used to generate the input file.
* `--o` : the path to the output file.
* `--r` : path to reference genome 

```
./indelible.py aggregate --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.sc_reads --b test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts --r hs37d5.fasta
```

**Note**: It is recommended to [retrain](#train) the RandomForest following this step with data that you have manually inspected. 

#### 3. Score

The **score** command scores positions based on the read information and sequence context:

* `--i` : path to the input file (the output of the *aggregate* command from previous step).
* `--o` : the path to the output file.

```
./indelible.py score --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored
```

**Note**: It is highly recommended if analysing a large amount of data at once to rebuild the InDelible frequency database. Please see the [Database](#database) command below for instructions.

#### 4. Blast

The **blast** command blasts longer clipped segments (20+ bp) to find matches elsewhere in the human genome and/or repeat database:
* `--i` :  path to the input file (the output of the *score* command from previous step).

This will automatically generate 3 files: a fasta file with the sequences to be blasted, 2 results files (one for repeat sequences, one for non-repeat sequences). These files are used in the next step.

```
./indelible.py blast --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored
```

#### 5. Annotate

The **annotate** command enriches the result with gene/exon annotations and merges the blast results with the position file:

* `--i` : path to the input file (output of score command after running the blast command).
* `--o` : path to output the annotated file.
* `--d` : path to the InDelible frequency database. A default frequency database from the Deciphering Developmental Disorders WES data is availible at `./Indelible/data/Indelible_db_10k.bed`. 

```
./indelible.py annotate --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored.annotated --d data/InDelible_db_10k.bed
```

#### 6. Denovo

One can then look for *de novo* mutation events using the **denovo** command:

***Note***: If maternal and/or paternal bam files are not supplied, *denovo* filtering will not be performed. This behaviour is intended to format non-trio data identically to trio data. If one of maternal **or** paternal bam is provided, InDelible will count coverage within that sample.

* `--c` : path to scored/annotated calls in the proband.
* `--m` : path to maternal CRAM/BAM file. [optional]
* `--p` : path to paternal CRAM/BAM file. [optional]
* `--o` : path to output file.

```
./indelible.py denovo --c test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.InDelible.tsv --m maternal.bam --p paternal.bam --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.InDelible.denovo.tsv
```

### Additional Commands

InDelible also includes several helper commands:

#### InDelible Complete

All steps in the InDelible calling pipeline can be performed in succession automatically with the **complete** command:

* `--i` : path to the input CRAM/BAM file.
* `--o` : path to directory where output files will be generated.
* `--r` : path to reference genome
* `--d` : path to the InDelible frequency [database](#database).
* `--m` : path to maternal CRAM/BAM file. [optional]
* `--p` : path to paternal CRAM/BAM file. [optional]
* `--keeptmp` : If this flag is given, intermediate files are kept. Otherwise these files will be removed once the analysis is finished.

```
./indelible.py complete --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o test_data/ --r hs37d5.fasta --keeptmp --m maternal.bam --p paternal.bam
```

For InDelible _de novo_ discovery behaviour when not providing maternal or paternal bams, please see [de novo](#6-denovo).

#### Database

The **database** command will generate the database required for the step [Annotate](#5-annotate). If analysing a large amount of data, it is highly recommended to run this command following running the [Score](#3-score) command.

* `--f` : file of files to merge to generate split read "allele frequencies"
* `--o` : output file to generate

```
ls InDelible_files/*.scored > fofn.txt
./indelible.py database --f fofn.txt --o InDelible_db.tsv
```

#### Train

The **train** command will run the active learning RandomForest on training data with annotations provided by manual inspection of data.

* `--i input training data`
* `--o output random forest`
* `--k number of samples to use for initial training and subsequent active selection [50].`
* `--s convergence parameter to stop learning at when accuracy does not improve by X% [0.01]`

The input data has identical columns to the output generated by [Aggregate](#2-aggregate), except with an additional column appended to the end of the file named "annotation". **Column names need to be identical!** See an example training file which can be used to regenerate the RandomForest used on the DDD data in `./data/`.

```
./indelible.py train --i input_training_data.txt --o output.pkl
```

**Note**: Remember to change the path to the RandomForest from `--o` in the ./config.yml file!