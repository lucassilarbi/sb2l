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

Dear ImGui is vendored under `3rd/imgui`; the GUI is off by default and enabled with
`-DSB2L_BUILD_GUI=ON`.

## Prerequisites

The core prerequisites (see "Building from source") plus GLFW and OpenGL (already in the
Docker image):

```bash
sudo apt-get install libglfw3-dev libgl1-mesa-dev
```

## Compile

From the repository root, on a fresh clone:

```bash
git submodule update --init --recursive
mkdir -p build && cd build
cmake -DSB2L_BUILD_GUI=ON ..
make sb2l_gui -j
```

The executable is produced at `build/gui/sb2l_gui`.

## Launch

A display is required (use the Docker X11 setup above).

```bash
./gui/sb2l_gui
```

Controls: pick the curve type / degree / number of control points in the "Controls" window,
together with the two sets the B-spline is built on, chosen independently of each other:

* **curve parameter**: the set the parameter `u` is taken in (`R`, `IR` or `Z`), i.e. how the
  basis functions themselves are evaluated. Changing it rebuilds the basis.
* **control points**: the set the control points are taken in (points, boxes or zonotopes),
  i.e. which of `eval_point` / `eval_box` / `eval_zonotope` is run on that basis.

All 9 pairs are valid: real control points over an interval parameter give the tube of the
curve, control boxes over a real parameter give their images with no parameter enclosure, and
so on.

In the "Canvas" window drag a control handle to edit, drag inside a box or zonotope to move it
and drag a corner to resize it, scroll to zoom, and drag empty space to pan.
