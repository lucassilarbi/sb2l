# SB2L: Set-Based B-spline Library

This library allows you to compute uniform, non-uniform, rational, non-rational set-based B-splines. The curve parameter and control points can be real, interval, or affine.

# Using Docker (Recommended)

Install docker [here](https://docs.docker.com/engine/install/).

Compile the image

```bash
sudo docker build -t sb2l .
```

Enable display

```
xhost +local:docker
```

Run the docker

```bash
sudo docker run -it --privileged -e DISPLAY=$DISPLAY -v /tmp/.X11-unix:/tmp/.X11-unix --net=host docker.io/sb2l
```

Launch the visualization (Optional): In another terminal:

- Interact with the docker
    ```
    docker exec -it <docker_name> bash
    ```

- Find the [Vibes](https://github.com/ENSTABretagneRobotics/VIBES) AppImage

    ```bash
    cd /root/sb2l/VIBes-viewer
    ```

- Start the viewer

    ```
    ./VIBes-0.2.3-linux.AppImage
    ```

Launch the example

```bash
./example
```

# Building from source

tested on Ubuntu22.04

Install prerequisites:

```bash
apt-get install cmake libgmp-dev python2.7 flex bison gcc g++ make pkg-config libfuse2 libqhull-dev
```

Update git submodules:

```bash
git submodule update --init --recursive
```

Install SB2L:

```bash
mkdir build && cd build
cmake ..
make
make install
```

# Building examples

Some examples are available in the ```examples``` folder

```bash
cd examples
mkdir build && cd build
cmake ..
make
```

Launch example

```bash
./example
```

# Visualization

Used the [Vibes](https://github.com/ENSTABretagneRobotics/VIBES) library. First launch the ```VIBes-0.2.3-linux.AppImage``` and then an examples.

```bash
cd VIBes-viewer
./VIBes-0.2.3-linux.AppImage
```