FROM ubuntu:20.04
FROM python:3.9
MAINTAINER Melissa Meredith, mmmeredi@ucsc.edu

RUN mkdir -p /home/apps


RUN cd /home/apps && \
    apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y vim git wget make build-essential cmake \
    tabix python3.8 python3.8-dev python3-pip \
    protobuf-compiler pkg-config libprotobuf-dev libjansson-dev libhts-dev libncurses-dev \
    libbz2-dev liblzma-dev zlib1g-dev autoconf libcurl4-openssl-dev curl libomp-dev libssl-dev python3-tk && \
    pip3 install modbamtools

WORKDIR /home/apps
RUN wget https://github.com/samtools/samtools/releases/download/1.15.1/samtools-1.15.1.tar.bz2 && \
    tar xvf samtools-1.15.1.tar.bz2 && \
    rm samtools-1.15.1.tar.bz2 && \
    cd samtools-1.15.1 && \
    ./configure && \
    make 

WORKDIR /home/apps
RUN wget https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24.tar.bz2 && \
    tar xvf minimap2-2.24.tar.bz2 && \
    rm minimap2-2.24.tar.bz2 && \
    cd minimap2-2.24 && \
    make


WORKDIR /home/apps
RUN wget https://github.com/attractivechaos/k8/releases/download/0.2.5/k8-0.2.5.tar.bz2 && \
    tar xvf k8-0.2.5.tar.bz2 && \
    rm k8-0.2.5.tar.bz2 && \
    cp k8-0.2.5/k8-`uname -s` k8-0.2.5/k8


WORKDIR /home/apps
RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.2.0/modkit_v0.2.0_centos7_x86_64.tar.gz && \
    tar -xvzf modkit_v0.2.0_centos7_x86_64.tar.gz

WORKDIR /home/apps
RUN git clone --recursive https://github.com/meredith705/methyl_read_grouping.git

ENV PATH="/home/apps/samtools-1.15.1:${PATH}"
ENV PATH="/home/apps/minimap2-2.24:/opt/minimap2/minimap2-2.24/misc:$PATH"
ENV PATH="/home/apps/k8-0.2.5:$PATH"
ENV PATH="/home/apps/dist:${PATH}"
ENV PATH="/home/apps/methyl_read_grouping/methlotype:$PATH"


RUN mkdir -p /data
WORKDIR /data

