FROM ubuntu:22.04

COPY . /root/sb2l

# Update and install in one layer, so the package index can never be a stale
# one cached from an earlier build.
RUN apt-get update && apt-get install -y --no-install-recommends \
    `# graphics, for the VIBes viewer and the editor` \
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
    ca-certificates \
    `# build` \
    cmake \
    libgmp-dev \
    python3 \
    flex \
    bison \
    gcc \
    g++ \
    make \
    libfuse2 \
    fuse \
    git \
    libglfw3-dev \
    libgl1-mesa-dev \
    && rm -rf /var/lib/apt/lists/*

# Build the library and the editor, optimised. The build tree is kept: the
# editor lives in it, and rebuilding from the image would mean configuring
# everything again.
WORKDIR /root/sb2l
RUN git submodule update --init --recursive
WORKDIR /root/sb2l/build
RUN cmake -DCMAKE_BUILD_TYPE=Release -DSB2L_BUILD_GUI=ON ..
RUN make -j
RUN make install

# Examples, ready to run (the viewer must be started first, see the README).
WORKDIR /root/sb2l/examples/build
RUN cmake -DCMAKE_BUILD_TYPE=Release ..
RUN make -j

WORKDIR /root/sb2l
CMD ["bash"]
