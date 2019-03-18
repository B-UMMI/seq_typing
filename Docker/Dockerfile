FROM ubuntu:16.04
MAINTAINER Miguel Machado <mpmachado@medicina.ulisboa.pt>
LABEL version="2.3-dev"

WORKDIR /NGStools/

RUN apt-get update

# -- General Dependencies ---
RUN apt-get install -y wget git python python3 python3-pip && \
    pip3 install setuptools

# -- seq_typing Dependencies ---

# --- Blast+ ---
# TODO: maybe return to v2.8.1 after fixing the headers name
#RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.8.1/ncbi-blast-2.8.1+-x64-linux.tar.gz && \
#    tar -xf ncbi-blast-2.8.1+-x64-linux.tar.gz && \
#    rm ncbi-blast-2.8.1+-x64-linux.tar.gz
#ENV PATH="/NGStools/ncbi-blast-2.8.1+/bin:${PATH}"
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    tar -xf ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    rm ncbi-blast-2.7.1+-x64-linux.tar.gz
ENV PATH="/NGStools/ncbi-blast-2.7.1+/bin:${PATH}"

# --- ReMatCh ---
RUN git clone https://github.com/B-UMMI/ReMatCh.git && \
    cd ReMatCh && \
    python3 setup.py install && \
    cd /NGStools
ENV PATH="/NGStools/ReMatCh/ReMatCh/src/samtools-1.3.1/bin:/NGStools/ReMatCh/ReMatCh/src/bcftools-1.3.1/bin:/NGStools/ReMatCh/ReMatCh/src/bowtie2-2.2.9:${PATH}"

# -- seq_typing --
# TODO: remove checkout
RUN git clone https://github.com/B-UMMI/seq_typing.git && \
    cd seq_typing && \
    git checkout gbs_types && \
    python3 setup.py install

# -- Clean --
RUN pip3 uninstall setuptools && \
    apt-get remove -y python3-pip && \
    apt-get autoclean -y

WORKDIR /data/
