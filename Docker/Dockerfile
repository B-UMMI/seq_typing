FROM continuumio/miniconda3:4.8.2-alpine
LABEL version="2.3-dev-2021"

# -- General Dependencies ---

RUN /opt/conda/bin/conda install -y git && \
    /opt/conda/bin/pip install setuptools && \
    /opt/conda/bin/conda clean -y --all

ENV PATH="/opt/conda/bin:${PATH}"

# -- seq_typing Dependencies ---

USER anaconda

WORKDIR /NGStools/

RUN ls -l /
RUN ls -l /NGStools/
RUN whoami

# --- Blast+ ---
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    tar -xf ncbi-blast-2.7.1+-x64-linux.tar.gz && \
    rm ncbi-blast-2.7.1+-x64-linux.tar.gz

# --- ReMatCh ---
# TODO: remove checkout
RUN git clone https://github.com/B-UMMI/ReMatCh.git && \
    cd ReMatCh && \
    git checkout dev && \
    python3 setup.py install && \
    cd /NGStools

# -- seq_typing --

# TODO: remove checkout
RUN git clone https://github.com/B-UMMI/seq_typing.git && \
    cd seq_typing && \
    git checkout gbs_types && \
    python3 setup.py install

ENV PATH="/NGStools/ReMatCh/ReMatCh/src/samtools-1.3.1/bin:/NGStools/ReMatCh/ReMatCh/src/bcftools-1.3.1/bin:/NGStools/ReMatCh/ReMatCh/src/bowtie2-2.2.9:/NGStools/ncbi-blast-2.7.1+/bin:${PATH}"

WORKDIR /data/
