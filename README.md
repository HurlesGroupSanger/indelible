# InDelible: Genomic Structural Variant Caller by Adaptive Training

## Table of Contents

1. [Development Team](#development-team)
2. [About InDelible](#about-indelible)
3. [Installation](#installation)

    a. [Required Software Dependencies](#required-software-dependencies)
    
    b. [Local Installation](#installing-indelible-on-a-local-machine)
    
    c. [Using vr-runner](#using-vr-runner) 
    
    d. [Docker & Singularity](#using-singularity-or-docker)
    
4. [Usage](#usage)
5. [Output](#output)
    
    a. [Primary TSV File](#primary-tsv-file)
    
    b. [Recommended Filtering](#recommended-filtering)

## Development Team

Alejandro Sifrim (Creator, Developer)

Eugene Gardner (Developer)

Diana Rajan (Experimental validation)

Elena Prigmore (Experimental validation)

Sarah Lindsay (Experimental validation)

Matthew Hurles (Group Leader)

We are affiliated with the [Wellcome Sanger Institute](http://www.sanger.ac.uk/science/groups/hurles-group), Cambridge, United Kingdom

## About InDelible

### Abstract

Structural Variations (SVs) are genetic differences greater than 50bps in size and have been identified as causative of 
diseases such as rare developmental disorders (DD). Patients presenting with DD are typically referred for chromosomal 
microarray (CMA) to identify large copy-number variants or for single gene or panel tests based on presenting symptoms.
Increasingly, patients for which a diagnosis is not forthcoming are additionally referred for whole exome sequencing (WES)
which is used to identify single nucleotide variants or small insertions/deletions (InDels). This leaves a class of
intermediate size deletions that are technically difficult to identify, hence patients with SVs undetectable by
conventional CMA or WES analysis often remain undiagnosed. To this end, we have developed a novel SV discovery approach,
‘InDelible’, and applied it to 13,438 probands with severe DD recruited as part of the Deciphering Developmental Disorders
(DDD) study. InDelible queries WES data to identify split read clusters within a gene set of interest and performs variant
quality-control utilizing an active learning methodology. Using InDelible we were able to find 59 previously undetected
variants among DDD probands, of which 89.8% (53) were in genes previously associated with DD, had phenotypes which
putatively match the conditions of the patient in which they were found, and were thus reported to the referring 
clinician. InDelible was particularly effective at ascertaining variants between 20-500 bps in size, and increased
the total number of causal variants identified by DDD in this size range by 46.4% (n = 26 variants). Of particular
interest were seven confirmed de novo SVs in the gene MECP2; these variants represent 35.0% of all de novo PTVs
in MECP2 among DDD patients and represent an enrichment of large, causal variants compared to other DD-associated
genes. InDelible provides a rapid framework for the discovery of likely pathogenic SVs and has the potential to
improve the diagnostic yield of WES.

### What is InDelible for?

InDelible was originally designed for the ascertainment of large InDels (>20bp) and Structural Variants from Whole Exome 
Sequencing (WES) data for which other ascertainment approaches have proven refractory. To reduce the search space, InDelible
also makes use of a "target gene list" containing genes that the end user is interested in.

### Potential Limitations of InDelible

InDelible is likely to be adaptable to a wide range of sequencing technologies, genetic architectures, and disease/conditions. 
However, we have not specifically tested InDelible with:

1. _Whole Genome Sequencing (WGS) Data_
    - While InDelible should in theory be able to identify variants from WGS, the number of 
reads that InDelible has to process would likely result in significantly increased run-times. 

2. _Other Diseases_
    - InDelible was designed as part of the Deciphering Developmental Disorders (DDD) study and, as such, was targeted to
    genes known to contribute to dominant developmental disorders. We have provided for the possibility that end-users may
    want to identify variants within other gene sets, but have not specifically performed variant discovery among a patient
    cohort to test this functionality.
    
3. _Other Genetic Architectures_
    - While causal _de novo_ variants play an important role in the genetic architecture of severe DD, recessive causes of DD are likewise 
    a major contributor. While we have focused our primary analysis on _de novo_ variation to try to maximise our discovery potential,
    we do not preclude the possibility that InDelible could also be used to identify homozygous or compound heterozygous 
    variants which could be plausibly linked to a patient's symptoms.
    
4. _Higher Allele Frequency Variants_
    - Likewise, since InDelible was developed as part of DDD, where the largest contributor to patient symptoms is high-penetrance
    _de novo_ variation, we focused our filtering to such variants. However, InDelible does report all high-confidence variants
    identified for each patient and could, in theory, be used to identify population level variation.    

### How to Cite InDelible

Eugene J. Gardner, Alejandro Sifrim, Sarah J. Lindsay, Elena Prigmore, Diana Rajan, Petr Danecek, Giuseppe Gallone, Ruth Y. Eberhardt, Hilary C. Martin, Caroline F. Wright, David R. FitzPatrick, Helen V. Firth, Matthew E. Hurles.
**InDelible: Detection and Evaluation of Clinically-relevant Structural Variation from Whole Exome Sequencing.** medRxiv (2020).

## Installation

### Required Software Dependencies

InDelible is written for Python2.7.* or Python3.7.*

InDelible requires the following software to be installed and in `$PATH`:

* [bedtools](https://bedtools.readthedocs.io/en/latest/)
  * **Note**: If using CRAM formated files with InDelible, bedtools v2.28 or later is required.
* [tabix](http://www.htslib.org/download/) 
* [bgzip](http://www.htslib.org/download/)

The python package Biopython requires a local install of [blast](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) in `$PATH` in order to function. This needs to be installed prior to [installing InDelible](#installing-InDelible).

### Installing InDelible on a Local Machine

To install InDelible:

1. Clone the git repo:

```
git clone https://github.com/eugenegardner/indelible.git
cd indelible/
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
## Download windowmasker (*hg19 ONLY*):
wget ftp://ftp.ncbi.nlm.nih.gov/blast/windowmasker_files/9606/wmasker.obinary

## Download the GRCh37 human reference and create the blast db:
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz.fai
makeblastdb -in hs37d5.fa -dbtype nucl
makeblastdb -in repeats.fasta -dbtype nucl
```

You can also make your own windowmasker resources. windowmasker binaries are provided as part of the blast+ distrobution. To make
the `*.obinary` format file one can use the commands: 

```
## Generate counts file:
windowmasker -in ref.fa -infmt blastdb -mk_counts -parse_seqids -out ref.counts

## Convert to .obinary:
windowmasker -convert -in ref.asnb -out ref.obinary -sformat obinary
```

See the [windowmasker documentation](https://www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/source/src/app/winmasker/) for more information.

6. Edit the config file to point to required data files and edit any default parameters:

An example config.yml file (hg19: `example_config.hg19.yml`, hg38: `example_config.hg38.yml`) is included in the top level directory. Edit this and save as `config.<VERSION>.yml`.

```
cd ..

## Provide the path to each required file at the top of the config.<VERSION>.yml
## Edit parameters at bottom of file.
vim example_config.yml

## cp/mv to config.<VERSION>.yml
cp example_config.<VERSION>.yml config.<VERSION>.yml
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

**Note**: The above data resources will only work if you run InDelible using the the human GRCh37/38 reference.

A brief note on Hg38 resources: 

* The InDelible MAF database provided with this distribution is a [liftOver](https://genome.ucsc.edu/cgi-bin/hgLiftOver)
of project resources generated using version hg19 of the human genome. As such, the position accuracy of variants contained there-in
should not be considered 100% accurate.
* The training data packaged with InDelible is technically agnostic to genome build. However, differences between alignment methods and/or genome builds
may result in slight differences in covariate importance. As such, it may be a good idea to train the random forest with data from study being assessed, if sufficiently large enough.
* As gnomAD has not generated pLI scores specifically using Hg38-aligned genomes, pLI values are from gnomADv2.1.1.

**Note**: All commands take the command-line option `--config`. The user **must** provide the path to a valid config yml to 
InDelible at runtime.

### Using vr-runner

InDelible comes packages with a script for use with the vr-runner packaged which automates the analysis of multiple 
samples. For information on how to install and use vr-runner, please see [this link](https://github.com/VertebrateResequencing/vr-runner).

The InDelible-specific vr-runner script is located at `./vr_runner_scripts/run-indelible`. This script requires a list of tab-delimited bam paths, where columns are:

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
/path/to/child.bam  /path/to/mum.bam    -
/path/to/child.bam  - /path/to/dad.bam
```

First setup the vr-runner config file:

```
./indelible/vr_runner_scripts/run-indelible +sampleconf > my.conf
## Change ALL parameters to point to the correct files:
## bams = list of bams made above
## ref = reference genome.fa
## config = config.yml
## indelible = path to the main indelible.py script
## database = path to the indelible db
```

A simple command for running InDelible with vr-runner is as follows:

```
./indelible/vr_runner_scripts/run-indelible -o ./output/ +config my.conf +maxjobs 1000 +loop 100 +retries -2
```

Resulting output will be in the `./output/` folder.

Additional commands for run-indelible:

```
   -b, --bams-list <file>      File with bam files
   -o, --outdir <dir>          Output directory
   -d, --database              Path to a precomputed MAF database [null]
```

***Note:*** If a filepath to a prebuilt MAF database (e.g. `./Indelible/data/`) is not provided for `-d/--database`, InDelible will rebuild the MAF database based only on samples included in `-b`!

### Using Singularity or Docker

We have developed both a [Docker](https://www.docker.com/) and derived [Singularity](https://sylabs.io/docs/) VM image to enable quick deployment of InDelible to both local and cloud compute platforms.
The InDelible docker image is available through [Docker Hub](https://hub.docker.com/) here: 

https://hub.docker.com/repository/docker/mercury/indelible
 
Please see Docker's documentation for information on running InDelible via the Docker framework.

We have also developed a separate github repository which hosts the associated InDelible Dockerfile as well as instructions for converting this Dockerfile into
a Singularity image:

https://github.com/wtsi-hgi/indelible-docker/tree/master

This Docker/Singularity image contains both an install of InDelible running on Python3.7 as well as all of the required 
library files as described in the [local install](#installing-indelible-on-a-local-machine) instructions. We describe here 
a basic command-line for running InDelible via Singularity. Please see the [Singularity documentation](https://sylabs.io/docs/) for more detail instructions
for running Singularity itself.

```
/software/singularity-v3.5.3/bin/singularity exec \
    -c \
    --cleanenv \
    --bind /path/to/local/storage/device/ \
    --pwd /usr/src/app/Indelible/ \
    /path/to/indelible_singularity_container.sif \
    indelible.py complete \
    --i /path/to/local/storage/device/proband.cram \
    --o /path/to/local/storage/device/ \
    --r data/hs37d5.fa \
    --d data/Indelible_db_10k.bed \
    --m /path/to/local/storage/device/mum.cram \
    --p /path/to/local/storage/device/dad.cram
```

The above command will run the [complete](#indelible-complete) InDelible SV discovery pipeline on the file `proband.cram`. 
Additional notes on the above command:

1. `-c` and `--cleanenv` clean your enviornment variables and home directory prior to running InDelible. These may not be strictly necessary, but are recommended to prevent any conflicts with the environment internal to the Singularity image.
2. `--bind` is meant to be used to mount your local storage folder to the Singularity image and should be pointed to the location on your local machine/cluster where your sequence files are stored.
3. `--pwd` is the location of the InDelible directory **within** the Singularity image. This part of the command line **must not** be changed.
4.  `--r` and `--d` point to reference files within the Singularity image. The other references files are also located within the Singularity instance at `4.  `--r` and `--d` point to reference files within the Singularity image. Prebuilt references files are also located within the Singularity instance at /usr/src/app/Indelible/data`.
  
**Big Note**: If using InDelible for anything beyond the built-in data files (i.e. on files aligned to other human references), 
you will need to generate a local version of the `config.hg19.yml` file, edit it to point to your own resource files, and then point 
indelible.py to it with `--config /path/to/config.new.yml`. The paths in `config.new.yml` can reflect a mix of both local paths and 
paths already within the Singularity instance.

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

**Note**: All commands also take the command-line option `--config` which overrides the default config.yml path. The user 
can either change the default config.yml, or provide a path to a different file with this option. See the section above for configuring the `config.yml` file.

These different steps can be performed by the different sub-commands:

#### 1. Fetch

The **fetch** command extracts the reads from the BAM file, it takes 2 arguments:

* `--i` : path to the input CRAM/BAM file.
* `--o` : path to output the clipped reads to.
* `--config` : path to the config.yml file.

```
./indelible.py fetch --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.sc_reads
```

#### 2. Aggregate

The **aggregate** merges information across reads towards a position-level view of the data:

* `--i` : path to the input file (the output of the *fetch* command from previous step).
* `--b` : path to the CRAM/BAM file used to generate the input file.
* `--o` : the path to the output file.
* `--r` : path to reference genome .
* `--config` : path to the config.yml file.

```
./indelible.py aggregate --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.sc_reads --b test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts --r hs37d5.fasta
```

**Note**: It is recommended to [retrain](#train) the RandomForest following this step with data that you have manually inspected. 

#### 3. Score

The **score** command scores positions based on the read information and sequence context:

* `--i` : path to the input file (the output of the *aggregate* command from previous step).
* `--o` : the path to the output file.
* `--config` : path to the config.yml file.

```
./indelible.py score --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored
```

**Note**: It is highly recommended if analysing a large amount of data at once to rebuild the InDelible frequency database. Please see the [Database](#database) command below for instructions.

#### 4. Blast

The **blast** command blasts longer clipped segments (20+ bp) to find matches elsewhere in the human genome and/or repeat database:
* `--i` :  path to the input file (the output of the *score* command from previous step).
* `--config` : path to the config.yml file.

This will automatically generate 3 files: a fasta file with the sequences to be blasted, 2 results files (one for repeat sequences, one for non-repeat sequences). These files are used in the next step.

```
./indelible.py blast --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored
```

#### 5. Annotate

The **annotate** command enriches the result with gene/exon annotations and merges the blast results with the position file:

* `--i` : path to the input file (output of score command after running the blast command).
* `--o` : path to output the annotated file.
* `--d` : path to the InDelible frequency database. A default frequency database from the Deciphering Developmental Disorders WES data is availible at `./Indelible/data/Indelible_db_10k.bed`. 
* `--config` : path to the config.yml file.

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
* `--config` : path to the config.yml file.

```
./indelible.py denovo --c test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.counts.scored.annotated --m maternal.bam --p paternal.bam --o test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam.indelible.denovo.tsv
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
* `--config` : path to the config.yml file.

```
./indelible.py complete --i test_data/DDD_MAIN5194229_Xchrom_subset_sorted.bam --o test_data/ --r hs37d5.fasta --keeptmp --m maternal.bam --p paternal.bam
```

For InDelible _de novo_ discovery behaviour when not providing maternal or paternal bams, please see [de novo](#6-denovo).

#### Database

The **database** command will generate the database required for the step [Annotate](#5-annotate). If analysing a large amount of data, it is highly recommended to run this command following running the [Score](#3-score) command.

* `--f` : file of files to merge to generate split read "allele frequencies"
* `--o` : output file to generate
* `--config` : path to the config.yml file.

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
* `--config` : path to the config.yml file.

The input data has identical columns to the output generated by [Aggregate](#2-aggregate), except with an additional column appended to the end of the file named "annotation". **Column names need to be identical!** See an example training file which can be used to regenerate the RandomForest used on the DDD data in `./data/`.

```
./indelible.py train --i input_training_data.txt --o output.pkl
```

**Note**: Remember to change the path to the RandomForest from `--o` in the ./config.yml file!

## Output

### Primary TSV File

The primary output from InDelible is the output file from the [_denovo_](#6-denovo) command (e.g. `<BAM_NAME>.bam.indelible.denovo.tsv`). This file has 41 columns, of which most are only relevant to how InDelible handles internal filtering with the adaptive learning model. What each column contains is listed in the table below:

| Column Name | Column # | Description |
|-------------|----------|-------------|
|chrom| 1 | Chromosome of breakpoint |
|position| 2 | Position of breakpoint |
|coverage| 3 | total number of reads covering breakpoint | 
|insertion_context| 4 | total number of insertions (cigar "I") in reads overlapping this breakpoint | 
|deletion_context| 5 | total number of deletions (cigar "D") in reads overlapping this breakpoint |
|sr_total| 6 | total number of split reads (cigar "S") in reads overlapping this breakpoint |
|sr_total_long| 7 | Number of reads with SR length ≥ MINIMUM_LENGTH_SPLIT_READ |
|sr_total_short| 8 | Number of reads with SR length < MINIMUM_LENGTH_SPLIT_READ |
|sr_long_5| 9 | sr_total_long for 5' end of reads |
|sr_short_5| 10 | sr_total_short for 5' end of reads |
|sr_long_3| 11 | sr_total_long for 3' end of reads |
|sr_short_3| 12 | sr_total_short for 3' end of reads |
|sr_entropy| 13 | Sequence entropy of the longest SR sequence given by the formula from Schmitt and Herzel (1997) |
|context_entropy| 14 | Sequence entropy of the ±20bp from the breakpoint position |
|entropy_upstream| 15 | Sequence entropy of the +20bp from the breakpoint position |
|entropy_downstream| 16 | Sequence entropy of the -20bp from the breakpoint position |
|sr_sw_similarity| 17 | Smith-Waterman based similarity of split reads from the breakpoint |
|avg_avg_sr_qual| 18 | Average sequence quality of split bases |
|avg_mapq| 19 | Average mapping quality of reads supporting the breakpoint |
|seq_longest| 20 | longest split sequence |
|pct_double_split| 21 | Number of reads with both 5' and 3' split reads |
|prob_N| 22 | Probability of the breakpoint being a false positive based on the adaptive learning model (1 - prob_Y) | 
|prob_Y| 23 | Probability of the breakpoint being a true positive based on the adaptive learning model |
|predicted| 24 | Is prob_Y > prob_N? | 
|ddg2p| 25 | Does this breakpoint intersect any genes given by `ddg2p_bed` file in config.yml |
|hgnc| 26 | Does this breakpoint intersect any genes given by `hgnc_file` in config.yml |
|hgnc_constrained| 27 | Does this breakpoint intersect any genes given by `hgnc_constrained` in config.yml | 
|exonic| 28 | Does this breakpoint intersect any exons given by `ensembl_exons` in config.yml |
|transcripts| 29 | What transcripts does this breakpoint intersect? If > 10 transcripts, will return 'multiple_transcripts' |
|maf| 30 | "Allele Frequency" based on the InDelible database provided with `--d`
|blast_hit| 31 | The coordinate given by BLAST for the longest SR of this breakpoint. If multiple BLAST hits, will be 'multi_hit', if overlaps a region from `repeatdb` in config.yml, will be repeats hit. |
|blast_strand| 32 | Strand of blast_hit |
|blast_identity| 33 | Percent identity of blast_hit |
|blast_dist| 34 | Distance to blast_hit from 'position' |
|blast_hgnc| 35 | Gene overlap of blast_hit |
|otherside| 36 | Coordinate of likely 5'/3' breakpoint if present |
|sv_type| 37 | If otherside found or blast_hit = "repeats_hit" potential SV type. Possible values are DUP (duplication), DEL (deletion), INS_<CLASS> (mobile element insertion), or SEGDUP_TRANS (segmental duplication or translocation). For INS, this will list the likely type of templated insertion from the *.fasta.hits_repeats file. SEGDUP_TRANS represents either a segmenetal duplication or translocation. As we cannot discern with short read data between a SEGDUP or translocation, we list both here. |
|mum_sr| 38 | Number of SRs in the bam/cram provided to --m with the same 'position' |
|dad_sr| 39 | Number of SRs in the bam/cram provided to --d with the same 'position' |
|mum_indel_context| 40 | Number of reads in the bam/cram provided to --m with cigar 'I/D' values |
|dad_indel_context| 41 | Number of reads in the bam/cram provided to --d with cigar 'I/D' values | 
|mum_cov| 42 | Coverage in in the bam/cram provided to --m |
|dad_cov| 43 | Coverage in in the bam/cram provided to --d |

### Recommended Filtering

As we have described in the [InDelible manuscript](#how-to-cite-indelible), we follow a strict filtering regimen to further refine our TSV output. 
Filters are as follows, and are meant to "funnel" variants down (i.e. each successive filter is applied following the 
previous filter). Filter names are as above to avoid confusion. Parentheses are applied as pseudocode to demonstrate order 
of operations.

1. maf ≤ 0.0004 **AND** avg_mapq ≥ 20
    * *Note*: The maf filter is dependent on using the maf data provided in one of the *data.zip files
2. (pct_double_split > 0.1 **AND** blast_hit != "no_hit") **OR** pct_double_split ≤ 0.1 
3. ddg2p != "NA" **AND** sr_total ≥ 5 **AND** exonic = "True"
4. mom_sr < 2 **AND** dad_sr < 2



