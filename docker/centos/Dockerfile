FROM centos:centos7 as base

ENV EPACTS_ROOT="/opt/epacts"

RUN yum -y install \
	centos-release-scl epel-release && \
    yum -y update && yum -y upgrade && \
    yum -y install \
	curl \
	groff \
	ghostscript \
	gnuplot \
	help2man \
	R R-studio \
	zlib

FROM base as builder

RUN yum -y install \
	autoconf \
	automake \
	gcc-c++ \
	libtool \
	python3 \
	python-pip \
	zlib-devel

WORKDIR /usr/src
COPY . .

RUN g++ --version && mkdir -p $EPACTS_ROOT && \
	autoreconf -f -i && ./configure --prefix $EPACTS_ROOT && \
	make && make install && rm -rf ${SRC_DIR}


FROM base as release

ENV PATH=$EPACTS_ROOT/bin:$PATH
COPY --from=builder $EPACTS_ROOT/ $EPACTS_ROOT/
ENTRYPOINT [ "epacts" ]
