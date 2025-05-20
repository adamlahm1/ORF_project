FROM python:3.12-slim

RUN pip install biopython matplotlib

RUN apt-get update && \
    apt-get install -y wget libgomp1 

WORKDIR /app

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    tar -xzvf ncbi-blast-2.16.0+-x64-linux.tar.gz && \
    mv ncbi-blast-2.16.0+ /usr/local/ && \
    rm ncbi-blast-2.16.0+-x64-linux.tar.gz

ENV PATH="/usr/local/ncbi-blast-2.16.0+/bin:$PATH"

RUN wget https://ftp.ncbi.nlm.nih.gov/blast/db/v5/swissprot.tar.gz && \
    mkdir -p /db/swissprot && mv swissprot* /db/swissprot && \
    tar -xzvf /db/swissprot/swissprot.tar.gz -C /db/swissprot

COPY . /app
