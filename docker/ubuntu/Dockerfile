FROM ubuntu:18.04 as base

ENV EPACTS_ROOT="/opt/epacts"

RUN apt -y update && apt -y upgrade && \
	DEBIAN_FRONTEND=noninteractive apt -y install --no-install-suggests \
	curl \
	groff \
	ghostscript \
	gnuplot \
	help2man \
	r-base-core \
	zlib1g


FROM base as builder

RUN DEBIAN_FRONTEND=noninteractive apt -y install --no-install-suggests\
	autoconf \
	build-essential \
	g++-5 \
	libtool \
	lsb-release \
	python \
	python3-pip

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 5 && \
	 update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-5 5

WORKDIR /usr/src
COPY . .

RUN g++ --version && mkdir -p $EPACTS_ROOT && \
	autoreconf -f -i && ./configure --prefix $EPACTS_ROOT && \
	make && make install && rm -rf ${SRC_DIR}


FROM base as release

ENV PATH=$EPACTS_ROOT/bin:$PATH
COPY --from=builder $EPACTS_ROOT/ $EPACTS_ROOT/
ENTRYPOINT [ "epacts" ]
