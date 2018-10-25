# Dockerfile for RNA-seq pipeline
FROM ubuntu:16.04
MAINTAINER Farrah Roy

RUN apt-get update && apt-get install -y software-properties-common && add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-7-jdk \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        unzip \
        tar \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*


#-----------------------------
# Pipeline components
#-----------------------------

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.5/htslib-1.5.tar.bz2 && \
    tar -xf htslib-1.5.tar.bz2 && rm htslib-1.5.tar.bz2 && cd htslib-1.5 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 && \
    tar -xf samtools-1.5.tar.bz2 && rm samtools-1.5.tar.bz2 && cd samtools-1.5 && \
    ./configure --with-htslib=/opt/htslib-1.5 && make && make install && make clean

# bamtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/pezmaster31/bamtools/archive/v2.4.1.tar.gz && \
    tar -xf v2.4.1.tar.gz && rm v2.4.1.tar.gz && cd bamtools-2.4.1 && mkdir build && cd build && cmake .. && make && make install && make clean
ENV LD_LIBRARY_PATH /usr/local/lib/bamtools:$LD_LIBRARY_PATH

#bedtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz && \
    tar -xf bedtools-2.26.0.tar.gz && \
    cd bedtools2 && \
    make

# HISAT2
RUN wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-Linux_x86_64.zip && \
    unzip hisat2-2.1.0-Linux_x86_64.zip
RUN cp -p hisat2-2.1.0/hisat2 hisat2-2.1.0/hisat2-* /usr/bin
ENV PATH /usr/bin/hisat2-2.1.0/hisat2 hisat2-2.1.0/hisat2-*:$PATH

# python modules
RUN pip3 install --upgrade pip
RUN pip install HTSeq

#Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-Src-0.36.zip && \
 unzip Trimmomatic-Src-0.36.zip && \
 cp -r trimmomatic-0.36 /usr/bin 

#install FastQC 0.11.5
RUN wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && unzip fastqc_v0.11.5.zip && chmod +x FastQC/fastqc 

# Install OpenJDK-8
RUN apt-get update && \
    apt-get install -y openjdk-8-jdk && \
    apt-get install -y ant && \
    apt-get clean;

# Fix certificate issues for OpenJDK
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f;

# Setup JAVA_HOME -- useful for docker commandline
ENV JAVA_HOME /usr/bin/java-8-openjdk-amd64/:$PATH
RUN export JAVA_HOME

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
    
ADD Trimming_protocol.sh /root/Trimming_protocol.sh
RUN chmod +x Trimming_protocol.sh

ENTRYPOINT ["/root/Trimming_protocol.sh"]