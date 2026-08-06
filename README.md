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

The curve can live in the plane or **in space**: the "dimension" combo switches between
a 2D canvas (pan/zoom) and a 3D one (orbit camera). SB2 reads the dimension as the number
of coordinate rows of the control points, so switching never rebuilds the symbolic basis,
it only re-runs the evaluation. Every set is editable in 3D too -- points, boxes with
their 8 corner handles, and zonotope generator tips free in any direction of space.

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
* **dimension**: whether the curve lives in the plane or in space. It is the number of
  coordinate rows handed to `eval_*`, so it is free of any rebuild.

All 9 pairs are valid: real control points over an interval parameter give the tube of the
curve, control boxes over a real parameter give their images with no parameter enclosure, and
so on.

In the "Canvas" window drag a control handle to edit, drag inside a box or zonotope to move it
and drag a corner to resize it, scroll to zoom, and drag empty space to pan.

In 3D, drag empty space to orbit and shift-drag (or right-drag) to pan. A grabbed handle moves
in the view-parallel plane through it, so orbit to reach the third axis. Result zonotopes are
drawn as their *exact* screen silhouette: the camera is orthographic on purpose, so the
projection of a zonotope is again a zonotope and the same O(m log m) boundary walk applies --
no 3D hull is ever built. Control zonotopes get their true facet mesh, enumerated directly
(every facet of a 3D zonotope is a parallelogram spanned by a pair of generators).

`--shot <file.ppm>` / `--shot3 <file.ppm>` render one screenshot of the 2D / 3D canvas and exit.
The geometry behind the 3D rendering is checked by `zonotope3d_checks`
(`cmake -DSB2L_BUILD_TESTS=ON`): exact support planes, enclosure of all 2^m vertices,
conservative generator reduction, and the projected mesh staying inside the projected silhouette.
