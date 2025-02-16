# Use the official Ubuntu 20.04 image
FROM ubuntu:20.04

# Set maintainer label (use LABEL instead of MAINTAINER)
LABEL maintainer="Shloka Negi, shnegi@ucsc.edu"

# Set environment variables to suppress timezone prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV DEBCONF_NONINTERACTIVE_SEEN=true
ENV TZ=America

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    wget \
    git \
    python \
    gcc \
    g++ \
    make \
    bzip2 \
    tabix \
    libxml2-dev \
    libssl-dev \
    libmariadb-dev \
    libcurl4-openssl-dev \
    apt-transport-https \
    software-properties-common \
    dirmngr \
    gpg-agent \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    && rm -rf /var/lib/apt/lists/*

# Set the timezone explicitly (adjust as needed)
RUN ln -fs /usr/share/zoneinfo/$TZ /etc/localtime && \
    dpkg-reconfigure -f noninteractive tzdata

# Install BedTools
WORKDIR /build/apps
ADD . /build/apps/
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz && \
    tar -zxvf bedtools-2.29.1.tar.gz && \
    cd bedtools2 && \
    make
# Add BedTools to the PATH
ENV PATH="${PATH}:/build/apps/bedtools2/bin"

# Install R and required R packages
WORKDIR /
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    r-base \
    && rm -rf /var/lib/apt/lists/*

RUN R -e "install.packages(c('dplyr', 'optparse', 'data.table'))"

# Create scripts directory and add the R script
ADD getMethylStats.R /

WORKDIR /