FROM ubuntu:16.04

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    cmake \
    ghostscript \
    git \
    gnuplot \
    groff \
    help2man \
    lsb-release \
    python \
    python-pip \
    r-base \
    rpm

RUN pip install cget

ENV SRC_DIR=/tmp/epacts-src
COPY . $SRC_DIR
WORKDIR $SRC_DIR

RUN rm -rf cget/
RUN cget install -DCMAKE_C_FLAGS="-fPIC" -DCMAKE_CXX_FLAGS="-fPIC" -f requirements.txt

RUN mkdir -p $SRC_DIR/build
WORKDIR $SRC_DIR/build
RUN rm -rf ./*
RUN cmake -DCMAKE_TOOLCHAIN_FILE=../cget/cget/cget.cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make install
