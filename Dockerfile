FROM ubuntu:20.04

WORKDIR /usr/src/app

# install apt dependencies
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && \
    apt-get -y install \
    ncbi-blast+ \
    python3-dev \
    gcc make libcurl4-gnutls-dev zlib1g-dev libncurses5-dev pkg-config \
    libncursesw5-dev liblzma-dev libz-dev g++ unzip gzip bwa libssl-dev \
    libbz2-dev liblzma-dev build-essential samtools bedtools tabix git python3-pip \
    libblas-dev libmkl-dev liblapack-dev

RUN apt-get -y install gfortran

# install repo and pip requirements
RUN git clone https://github.com/HurlesGroupSanger/indelible.git
WORKDIR /usr/src/app/indelible
RUN pip install 'numpy==1.17.2' 'cython==0.29.13' \
    && pip install -r requirements.txt

# install other dependencies: bedtools tabix htslib
# RUN conda install -c bioconda htslib bedtools cython tabix -y
#RUN conda install -c htslib bedtools tabix -y

# Unzip required data files
WORKDIR /usr/src/app/indelible/data/
RUN unzip data_hg19.zip \
    && unzip -n data_hg38.zip

# indelible.py looks for python ("#!/usr/bin/env python"), but it requires python3
RUN ln -s /usr/bin/python3 /usr/bin/python && chmod a+rx /usr/bin/python

## Download the GRCh37 human reference and create blastdb.
ADD https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz hs37d5.fa.gz
# Have to do exit 0 because of odd exit behaviour from gunzip
RUN gunzip -q hs37d5.fa.gz; exit 0
RUN samtools faidx hs37d5.fa
RUN bwa index hs37d5.fa

# Build repeats DB
RUN makeblastdb -in repeats.fasta -dbtype nucl

## Now do the same for GRCh38
ADD https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz GRCh38_full_analysis_set.fa.gz
RUN gunzip GRCh38_full_analysis_set.fa.gz \
    && samtools faidx GRCh38_full_analysis_set.fa
RUN bwa index GRCh38_full_analysis_set.fa

# Check everything is in the right place
RUN ls -ltra

# add config and set final WORKDIR
WORKDIR /usr/src/app/indelible
RUN mv example_config.hg19.yml config.hg19.yml \
    && mv example_config.hg38.yml config.hg38.yml

ENV PATH="/usr/src/app/indelible:${PATH}"
CMD [ "python3", "./indelible.py" ]
