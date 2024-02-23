ARG SAMTOOLS_VERSION="1.17"
ARG BCFTOOLS_VERSION="1.17"
ARG HTSLIB_VERSION="1.17"
ARG BEDTOOLS_VERSION="2.31.0"
ARG IVAR_VERSION="1.4.2"
ARG VARSCAN_VERSION="2.4.6"
ARG FREEBAYESS_VERSION="1.3.6"
ARG SNPEFF_VERSION="4_3t"
ARG PANGOLIN_VERSION="4.3.1"
ARG NEXTCLADE_VERSION="2.14.0"
ARG FASTQC_VERSION="0.12.1"
ARG TRIMMOMATIC_VERSION="0.39"
ARG BWA_VERSION="0.7.17"
ARG PICARD_VERSION="2.27.5"
ARG MINIMAP2_VERSION="2.26"
ARG GOFASTA_VERSION="1.2.1"
ARG USHER_VERSION="0.6.2"
ARG ISA_VERSION="2.30.0"
ARG TBB_VERSION="2019_U9"
ARG TBB_VERSION2="v2021.10.0"

# This is multistage, parallel build. We first build a base image.
# Then, in parallel, we build all the tools we need.
# And then in the end we copy only necessary files to production image.

# builder-base
FROM python:3.11.8-bookworm AS builder-base

ENV LANG en_US.UTF-8
ENV VIRTUAL_ENV=/opt/venv
ENV PATH=/$VIRTUAL_ENV/bin:$PATH:usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

# We install some basic tools, build tools and libraries
RUN apt update && \
    apt install -y --no-install-recommends vim \
                                           unzip \
                                           bc \
                                           less \
                                           tree \
                                           ncdu \
                                           parallel  \
                                           curl && \
    apt install -y --no-install-recommends git \
                                           build-essential \
                                           cmake \
                                           autoconf \
                                           meson \
                                           ninja-build \
                                           autotools-dev && \
    apt install -y --no-install-recommends libncurses5-dev \
                                           libbz2-dev \
                                           liblzma-dev \
                                           libcurl4-gnutls-dev \
                                           libgsl-dev \
                                           libvcflib-tools \
                                           libvcflib-dev \
                                           libseqlib-dev \
                                           zlib1g-dev

# Here we create a new python virtual environment. We can later just copy it to production container
# It is "activated" by adding it to PATH
RUN python3 -m venv $VIRTUAL_ENV

# builder-samtools
FROM builder-base AS builder-samtools
ARG SAMTOOLS_VERSION
WORKDIR /downloads
RUN curl -fsSL "https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2" \
         -o "samtools-${SAMTOOLS_VERSION}.tar.bz2" && \
    tar -xf samtools-${SAMTOOLS_VERSION}.tar.bz2 && \
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2
WORKDIR /downloads/samtools-${SAMTOOLS_VERSION}
RUN ./configure --prefix=/opt/samtools && \
    make -j `nproc` && \
    make -j `nproc` install
ENV PATH=/opt/samtools/bin:${PATH}

# builder-bcftools
FROM builder-base as builder-bcftools
ARG BCFTOOLS_VERSION
WORKDIR /downloads
RUN curl -fsSL "https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2" -o "bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    tar -xf "bcftools-${BCFTOOLS_VERSION}.tar.bz2" && \
    rm "bcftools-${BCFTOOLS_VERSION}.tar.bz2"
WORKDIR /downloads/bcftools-${BCFTOOLS_VERSION}
RUN ./configure --prefix=/opt/bcftools && \
    make -j `nproc` && \
    make -j `nproc` install
ENV PATH=/opt/bcftools/bin:${PATH}


# builder-htslib
FROM builder-base as builder-htslib
ARG HTSLIB_VERSION
WORKDIR /downloads
RUN curl -fsSL "https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2" -o "htslib-${HTSLIB_VERSION}.tar.bz2" && \
    tar -xf htslib-${HTSLIB_VERSION}.tar.bz2 && \
    rm htslib-${HTSLIB_VERSION}.tar.bz2
WORKDIR /downloads/htslib-${HTSLIB_VERSION}
RUN ./configure --prefix="/opt/htslib"
RUN make -j `nproc` && \
    make -j `nproc` install
ENV PATH=/opt/htslib/bin:${PATH}

# builder-bedtools
FROM builder-base as builder-bedtools
ARG BEDTOOLS_VERSION
WORKDIR /downloads
RUN curl -fsSL "https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz" -o "bedtools-${BEDTOOLS_VERSION}.tar.gz" && \
    tar -zxf bedtools-${BEDTOOLS_VERSION}.tar.gz && \
    rm bedtools-${BEDTOOLS_VERSION}.tar.gz
WORKDIR /downloads/bedtools2
RUN make -j `nproc` && \
    mkdir -p /opt/bedtools2/ && \
    cp -r /downloads/bedtools2/bin /opt/bedtools2/bin
ENV PATH=/opt/bedtools2/bin:${PATH}

# builder-ivar
FROM builder-htslib as builder-ivar
ARG IVAR_VERSION
ENV HTSLIB_PATH="/opt/htslib"
WORKDIR /downloads
RUN curl -fsSL "https://github.com/andersen-lab/ivar/archive/refs/tags/v${IVAR_VERSION}.tar.gz" -o "v${IVAR_VERSION}.tar.gz" && \
    tar -zxf v${IVAR_VERSION}.tar.gz && \
    rm v${IVAR_VERSION}.tar.gz
WORKDIR /downloads/ivar-${IVAR_VERSION}
RUN ./autogen.sh && \
    CPPFLAGS=-I/opt/htslib/include ./configure --with-hts=${HTSLIB_PATH} --prefix=/opt/ivar && \
    make -j `nproc` && \
    make -j `nproc` install
ENV LD_LIBRARY_PATH=/opt/htslib/lib:${LD_LIBRARY_PATH}
ENV PATH=/opt/ivar/bin:${PATH}

# builder-varscan
FROM builder-base as builder-varscan
ARG VARSCAN_VERSION
WORKDIR /opt/varscan
RUN curl -fsSL "https://github.com/dkoboldt/varscan/releases/download/v${VARSCAN_VERSION}/VarScan.v${VARSCAN_VERSION}.jar" -o "VarScan.v${VARSCAN_VERSION}.jar" && \
    ln -s VarScan.v${VARSCAN_VERSION}.jar varscan.jar


# builder-lofreq
FROM builder-base as builder-lofreq
# TODO: git is always downloading latest version. We should use specific version, tag or commit.
ARG HTSLIB_PATH
WORKDIR /downloads/lofreq
RUN apt install -y --no-install-recommends libgsl27 && \
    git clone https://github.com/CSB5/lofreq.git /downloads/lofreq && \
    ./bootstrap && \
    ./configure --with-htslib=${HTSLIB_PATH} --prefix=/opt/lofreq && \
    make -j `nproc` && \
    make -j `nproc` install
ENV PATH="/opt/lofreq/bin:$PATH"

# builder-vcftools
FROM builder-base as builder-vcftools
# TODO: git is always downloading latest version. We should use specific version, tag or commit.
WORKDIR /downloads/vcftools
RUN git clone https://github.com/vcftools/vcftools.git /downloads/vcftools && \
    ./autogen.sh && \
    ./configure --prefix=/opt/vcftools && \
    make -j `nproc` && \
    make -j `nproc` install
ENV PATH="/opt/vcftools/bin:$PATH"

# builder-freebayes
FROM builder-base as builder-freebayes
ARG FREEBAYESS_VERSION
WORKDIR /opt/freebayes/bin
RUN curl -fsSL "https://github.com/freebayes/freebayes/releases/download/v${FREEBAYESS_VERSION}/freebayes-${FREEBAYESS_VERSION}-linux-amd64-static.gz" -o "freebayes-${FREEBAYESS_VERSION}-linux-amd64-static.gz" && \
    gunzip freebayes-${FREEBAYESS_VERSION}-linux-amd64-static.gz && \
    chmod 755 freebayes-${FREEBAYESS_VERSION}-linux-amd64-static && \
    ln -s freebayes-${FREEBAYESS_VERSION}-linux-amd64-static freebayes
ENV PATH="/opt/freebayes/bin:$PATH"

# builder-snpeff
FROM builder-base as builder-snpeff
WORKDIR /opt
ADD third-party/snpEff.tar.gz /opt
COPY data/snpEff /opt/snpEff

# builder-python-stage1-pangolin
FROM builder-base as builder-python-stage1-pangolin
ARG PANGOLIN_VERSION
COPY requirements.txt /
RUN pip install -r /requirements.txt
WORKDIR /downloads
RUN wget https://github.com/cov-lineages/pangolin/archive/refs/tags/v${PANGOLIN_VERSION}.tar.gz  && \
    tar -zxvf v${PANGOLIN_VERSION}.tar.gz && \
    rm v${PANGOLIN_VERSION}.tar.gz
WORKDIR /downloads/pangolin-${PANGOLIN_VERSION}
RUN python setup.py install && \
    pangolin --update

# builder-nextclade
FROM builder-base as builder-nextclade
ARG NEXTCLADE_VERSION
WORKDIR /opt/nextclade/bin
RUN curl -fsSL "https://github.com/nextstrain/nextclade/releases/download/${NEXTCLADE_VERSION}/nextclade-x86_64-unknown-linux-gnu" -o "nextclade" && \
    chmod +x nextclade
RUN mkdir -p /SARS-CoV2/
# Switch comments to download nextclade dataset insted of copying it from local
ADD data/nextclade/SARS-CoV2.tar.gz /home
RUN nextclade dataset get --name='sars-cov-2' --output-dir=/SARS-CoV2/nextclade/
ENV PATH="/opt/nextclade/bin:$PATH"

# builder-fastqc
FROM builder-base as builder-fastqc
ARG FASTQC_VERSION
WORKDIR /opt
RUN curl -fsSL "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v${FASTQC_VERSION}.zip" -o "fastqc_v${FASTQC_VERSION}.zip" && \
    unzip "fastqc_v${FASTQC_VERSION}.zip" && \
    rm "fastqc_v${FASTQC_VERSION}.zip"
ENV PATH="/opt/FastQC:$PATH"

# builder-trimmomatic
FROM builder-base as builder-trimmomatic
ARG TRIMMOMATIC_VERSION
WORKDIR /opt
RUN curl -fsSL "https://github.com/usadellab/Trimmomatic/files/5854859/Trimmomatic-${TRIMMOMATIC_VERSION}.zip" -o "Trimmomatic-${TRIMMOMATIC_VERSION}.zip" && \
    unzip "Trimmomatic-${TRIMMOMATIC_VERSION}.zip" && \
    rm "Trimmomatic-${TRIMMOMATIC_VERSION}.zip" && \
    mv Trimmomatic-${TRIMMOMATIC_VERSION} trimmomatic && \
    ln -s /opt/trimmomatic/trimmomatic-${TRIMMOMATIC_VERSION}.jar /opt/trimmomatic/trimmomatic.jar

# builder-bwa
FROM builder-base as builder-bwa
#TODO Add commit as argument
WORKDIR /download
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /download/bwa
RUN git checkout 139f68f && \
    make -j `nproc` && \
    mkdir -p /opt/bwa/bin && \
    cp bwa /opt/bwa/bin
ENV PATH="/opt/bwa/bin:$PATH"

# builder-picard
FROM builder-base as builder-picard
ARG PICARD_VERSION
WORKDIR /opt/picard
RUN curl -fsSL "https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar" -o "picard.jar"

# builder-minimap2
FROM builder-base as builder-minimap2
ARG MINIMAP2_VERSION
WORKDIR /download
RUN curl -fsSL "https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}.tar.bz2" -o "minimap2-${MINIMAP2_VERSION}.tar.bz2" && \
    tar -jxvf "minimap2-${MINIMAP2_VERSION}.tar.bz2" && \
    rm "minimap2-${MINIMAP2_VERSION}.tar.bz2"
WORKDIR /download/minimap2-${MINIMAP2_VERSION}
RUN make -j `nproc` && \
    mkdir -p /opt/minimap2/bin && \
    cp minimap2 /opt/minimap2/bin
ENV PATH="/opt/minimap2/bin:$PATH"

# builder-gofasta
FROM builder-base as builder-gofasta
ARG GOFASTA_VERSION
WORKDIR /opt/gofasta/bin
RUN curl -fsSL "https://github.com/virus-evolution/gofasta/releases/download/v${GOFASTA_VERSION}/gofasta-linux-amd64" -o "gofasta" && \
    chmod +x gofasta
ENV PATH="/opt/gofasta/bin:$PATH"

# builder-faToVcf
# TODO: We consider this as unstable resource, so we should store it locally in third-party directory.
FROM builder-base as builder-faToVcf
WORKDIR /opt/blat/bin
RUN curl -fsSL "https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faToVcf" -o "faToVcf" && \
    chmod +x faToVcf
ENV PATH="/opt/blat/bin:$PATH"

# builder-isa-l
# This lib is a dependency of one of Usher binaries, we don't use it but it is still needed for building usher.
# We do not copy that unnecesary binary to production image, and also we don't copy this i-sal lib, either.
FROM builder-base as builder-isa-l
ARG ISA_VERSION
WORKDIR /downloads
RUN apt install -y --no-install-recommends nasm && \
    curl -fsSL "https://github.com/intel/isa-l/archive/refs/tags/v${ISA_VERSION}.tar.gz" -o "v${ISA_VERSION}.tar.gz" && \
    tar -zxf v${ISA_VERSION}.tar.gz && \
    rm v${ISA_VERSION}.tar.gz
WORKDIR /downloads/isa-l-${ISA_VERSION}
RUN ./autogen.sh && \
    ./configure --prefix "/opt/isa-l" && \
    make -j $(nproc) && \
    make -j $(nproc) install
WORKDIR /home
RUN rm -rf /downloads

# builder-tbb
FROM builder-base as builder-tbb
ARG TBB_VERSION2
WORKDIR /downloads
RUN git clone --depth 1 --branch ${TBB_VERSION2} https://github.com/oneapi-src/oneTBB.git
WORKDIR /downloads/oneTBB/build
RUN cmake -DCMAKE_INSTALL_PREFIX=/opt/oneTBB -DTBB_TEST=OFF .. && \
    make -j $(nproc) && \
    make -j $(nproc) install
CMD bash

# builder-usher
FROM builder-base as builder-usher
ARG USHER_VERSION
ARG TBB_VERSION
COPY --from=builder-isa-l /opt/isa-l /opt/isa-l
WORKDIR /downloads
#RUN apt -y --no-install-recommends install libprotoc-dev \
#                                           libopenmpi-dev \
#                                           libboost-dev \
#                                           libboost-program-options-dev \
#                                           libboost-iostreams-dev \
#                                           libboost-filesystem-dev \
#                                           libboost-date-time-dev
RUN git clone --depth 1 --branch v${USHER_VERSION} https://github.com/yatisht/usher.git && \
    cd usher && \
    mkdir build
WORKDIR /downloads/usher/build
RUN curl -fsSL "https://github.com/oneapi-src/oneTBB/archive/${TBB_VERSION}.tar.gz" -o "${TBB_VERSION}.tar.gz" && \
    tar -zxf ${TBB_VERSION}.tar.gz && \
    rm ${TBB_VERSION}.tar.gz

RUN apt -y --no-install-recommends install gfortran-12 \
                                          ibverbs-providers \
                                          libboost1.74-dev \
                                          libboost-date-time1.74-dev \
                                          libboost-date-time1.74.0 \
                                          libboost-date-time-dev \
                                          libboost-dev \
                                          libboost-filesystem1.74-dev \
                                          libboost-filesystem1.74.0 \
                                          libboost-filesystem-dev \
                                          libboost-iostreams1.74-dev \
                                          libboost-iostreams1.74.0 \
                                          libboost-iostreams-dev \
                                          libboost-program-options1.74-dev \
                                          libboost-program-options1.74.0 \
                                          libboost-program-options-dev \
                                          libboost-regex1.74-dev \
                                          libboost-regex1.74.0 \
                                          libboost-serialization1.74-dev \
                                          libboost-serialization1.74.0 \
                                          libboost-system1.74-dev \
                                          libboost-system1.74.0 \
                                          libfabric1 \
                                          libgfortran-12-dev \
                                          libhwloc15 \
                                          libhwloc-dev \
                                          libhwloc-plugins \
                                          libibverbs1 \
                                          libibverbs-dev \
                                          libjs-jquery-ui \
                                          libjs-jquery \
                                          libmunge2 \
                                          libnl-3-200 \
                                          libnl-3-dev \
                                          libnl-route-3-200 \
                                          libnl-route-3-dev \
                                          libnuma-dev \
                                          libopenmpi3 \
                                          libopenmpi-dev \
                                          libpciaccess0 \
                                          libpmix2 \
                                          libpmix-dev \
                                          libprotobuf32 \
                                          libprotobuf-dev \
                                          libprotobuf-lite32 \
                                          libprotoc32 \
                                          libprotoc-dev \
                                          libpsm2-2 \
                                          libpsm-infinipath1 \
                                          librdmacm1 \
                                          libucx0 \
                                          libxnvctrl0 \
                                          ocl-icd-libopencl1 \
                                          openmpi-bin \
                                          openmpi-common \
                                          protobuf-compiler


RUN curl -fsSL "https://github.com/oneapi-src/oneTBB/archive/${TBB_VERSION}.tar.gz" -o "${TBB_VERSION}.tar.gz" && \
    tar -zxf ${TBB_VERSION}.tar.gz && \
    rm ${TBB_VERSION}.tar.gz && \
    cmake -DCMAKE_INSTALL_PREFIX=/opt/usher \
          -DTBB_DIR=${PWD}/oneTBB-2019_U9 \
          -DISAL_LIB=/opt/isa-l/lib/libisal.so \
          -DCMAKE_CXX_FLAGS=-I/opt/isa-l/include .. && \
    make -j $(nproc) && \
    make -j $(nproc) install && \
    mkdir -p /opt/tbb/lib && \
    cp /downloads/usher/build/tbb_cmake_build/tbb_cmake_build_subdir_release/libtbb_preview.so.2 /opt/tbb/lib


# production
FROM python:3.11.8-slim-bookworm AS production

LABEL maintainer="Michal Lazniewski <mlazniewski@pzh.gov.pl>"
LABEL maintainer="Michał Kadlof <mkadlof@pzh.gov.pl>"

RUN echo 'alias ls="ls --color --group-directories-first"' >> /root/.bashrc && \
    echo "PS1='\[\033[36m\][\[\033[m\]\[\033[34m\]\[\e[1m \u@\h \e[0m\] \[\033[32m\]\W\[\033[m\]\[\033[36m\]]\[\033[m\]# '" >> /root/.bashrc

ENV PATH="/opt/usher/bin:/opt/blat/bin:/opt/gofasta/bin:/opt/minimap2/bin:/opt/bwa/bin:/opt/FastQC:/opt/nextclade/bin:/opt/freebayes/bin:/opt/vcftools/bin:/opt/lofreq/bin:/opt/ivar/bin:/opt/bedtools2/bin:/opt/htslib/bin:/opt/bcftools/bin:/opt/samtools/bin:/opt/venv/bin:/usr/local/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
ENV LD_LIBRARY_PATH="/opt/htslib/lib:/opt/tbb/lib:/opt/isa-l/lib:${LD_LIBRARY_PATH}"

RUN apt update && \
    apt install -y --no-install-recommends bc  \
                                           openjdk-17-jre-headless \
                                           libdeflate-dev \
                                           libcurl3-gnutls\
                                           libcurl4 \
                                           mafft \
                                           libboost-program-options1.74.0 \
                                           libboost-iostreams1.74.0 \
                                           libboost-filesystem1.74.0 \
                                           libprotobuf32 \
                                           libgsl27 \
                                           procps

COPY --from=builder-base /opt /opt
COPY --from=builder-samtools /opt/samtools /opt/samtools
COPY --from=builder-bcftools /opt/bcftools /opt/bcftools
COPY --from=builder-htslib /opt/htslib /opt/htslib
COPY --from=builder-bedtools /opt/bedtools2 /opt/bedtools2
COPY --from=builder-ivar /opt/ivar /opt/ivar
COPY --from=builder-varscan /opt/varscan /opt/varscan
COPY --from=builder-lofreq /opt/lofreq /opt/lofreq
COPY --from=builder-vcftools /opt/vcftools /opt/vcftools
COPY --from=builder-freebayes /opt/freebayes /opt/freebayes
COPY --from=builder-snpeff /opt/snpEff /opt/snpEff
COPY --from=builder-python-stage1-pangolin /opt/venv/ /opt/venv
COPY --from=builder-nextclade /opt/nextclade /opt/nextclade
COPY --from=builder-fastqc /opt/FastQC /opt/FastQC
COPY --from=builder-trimmomatic /opt/trimmomatic /opt/trimmomatic
COPY --from=builder-bwa /opt/bwa /opt/bwa
COPY --from=builder-picard /opt/picard /opt/picard
COPY --from=builder-minimap2 /opt/minimap2 /opt/minimap2
COPY --from=builder-gofasta /opt/gofasta /opt/gofasta
COPY --from=builder-faToVcf /opt/blat /opt/blat
# Usher builds with several binaries, but pelines use only a one. We copy only that one.
COPY --from=builder-usher /opt/usher/bin/usher /opt/usher/bin/usher
COPY --from=builder-usher /opt/tbb /opt/tbb

ADD data/generic/generic_data.tar.gz /home/SARS-CoV2
ADD data/nextclade/SARS-CoV2.tar.gz /home

WORKDIR /home
