FROM ubuntu:16.04
MAINTAINER Miguel Machado <mpmachado@medicina.ulisboa.pt>
LABEL version="2.2-02"

WORKDIR /NGStools/

RUN apt-get update

# -- General Dependencies ---
RUN apt-get install -y wget git python python3 python3-pip && \
    pip3 install setuptools

# -- seq_typing Dependencies ---
# --- Blast+ ---
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    tar -xf ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    rm ncbi-blast-2.7.1+-x64-linux.tar.gz
ENV PATH="/NGStools/ncbi-blast-2.7.1+/bin:${PATH}"
# --- ReMatCh ---
RUN git clone https://github.com/B-UMMI/ReMatCh.git && \
    cd ReMatCh && \
    # git checkout bowtie_algo && \
    python3 setup.py install && \
    cd /NGStools

# -- seq_typing --
RUN git clone https://github.com/B-UMMI/seq_typing.git && \
    cd seq_typing && \
    # git checkout bowtie_index && \
    python3 setup.py install

# -- Clean --
RUN pip3 uninstall setuptools && \
    apt-get remove -y python3-pip && \
    apt-get autoclean -y

WORKDIR /data/
