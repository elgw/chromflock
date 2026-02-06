FROM ubuntu:24.04

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils \
    ssh \
    apt-transport-https ca-certificates \
    cmake \
    libz-dev \
    libreadline-dev \
    make \
    gcc \
    git

RUN git clone https://github.com/elgw/chromflock.git
WORKDIR /chromflock/
RUN git checkout dev
RUN make -C src/lua-5.3.5/ linux
RUN mkdir build
WORKDIR /chromflock/build/
RUN cmake ..
RUN make
RUN ./chromflock --help
