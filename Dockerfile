FROM ubuntu:22.04

COPY . /root/sb2l

RUN apt-get update

# Graphics
RUN apt-get install -y \
    libx11-xcb1 \
    libxcb1 \
    libxcb-xinerama0 \
    libxcb-cursor0 \
    libxkbcommon-x11-0 \
    libglib2.0-0 \
    libgl1 \
    libsm6 \
    libxext6 \
    libxrender1 \
    libx11-6 \
    libxau6 \
    libxdmcp6 \
    fontconfig \
    libfreetype6 \
    libfontconfig1 \
    ca-certificates

# Dependencies
RUN apt-get install -y \
    cmake \ 
    libgmp-dev \
    python2.7 \
    flex \
    bison \
    gcc \
    g++ \
    make \
    libfuse2 \
    libqhull-dev \
    fuse \
    git

# Build
WORKDIR /root/sb2l
RUN git submodule update --init --recursive
WORKDIR /root/sb2l/build
RUN cmake ..
RUN make -j
RUN make install
WORKDIR /root/sb2l/
RUN rm -rf build
WORKDIR /root/sb2l/examples/build
RUN cmake ..
RUN make

CMD ["bash"]