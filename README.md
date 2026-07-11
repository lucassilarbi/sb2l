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
sudo apt-get install cmake libgmp-dev python3 flex bison gcc g++ make pkg-config libfuse2 libqhull-dev
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
sudo make install
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
./basic_example
```

# Visualization

Used the [Vibes](https://github.com/ENSTABretagneRobotics/VIBES) library. First launch the ```VIBes-0.2.3-linux.AppImage``` and then an examples.

```bash
cd VIBes-viewer
./VIBes-0.2.3-linux.AppImage
```

# Interactive GUI

An optional Dear ImGui editor (`gui/`) displays a B-spline of any type and lets you
edit it directly on its control points: drag a control point, drag inside a control box
or zonotope to move it, and drag a corner to resize it. Zonotope rendering uses a fast
O(m log m) boundary so evaluation stays interactive.

Extra prerequisites (already in the Docker image):

```bash
sudo apt-get install libglfw3-dev libgl1-mesa-dev
```

Build and run (Dear ImGui is vendored under `3rd/imgui`; the build is off by default):

```bash
mkdir build && cd build
cmake -DSB2L_BUILD_GUI=ON ..
make sb2l_gui
./gui/sb2l_gui
```

A display is required (use the Docker X11 setup above).
