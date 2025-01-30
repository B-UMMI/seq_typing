FROM continuumio/miniconda3:4.8.2-alpine

# -- Dependencies ---

RUN /opt/conda/bin/conda install -y -c bioconda python=3.9 blast=2.10.1 biopython && \
    /opt/conda/bin/pip install setuptools && \
    /opt/conda/bin/conda clean -y --all

ENV PATH="/opt/conda/bin:${PATH}"

# -- seq_typing Dependencies ---

WORKDIR /NGStools/

USER root
RUN apk add bash git libstdc++

# --- ReMatCh ---
RUN git clone https://github.com/B-UMMI/ReMatCh.git && \
    cd ReMatCh && \
    git checkout feature/save_ref_length && \
    python3 setup.py install && \
    cd /NGStools

# -- seq_typing --
RUN git clone https://github.com/B-UMMI/seq_typing.git && \
    cd seq_typing && \
    python3 setup.py install

ENV PATH=/NGStools/ReMatCh/ReMatCh/src/samtools-1.3.1/bin:/NGStools/ReMatCh/ReMatCh/src/bcftools-1.3.1/bin:/NGStools/ReMatCh/ReMatCh/src/bowtie2-2.2.9:${PATH}

WORKDIR /data/

CMD ["/bin/bash"]
