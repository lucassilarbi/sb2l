# SB2L: Set-Based B-spline Library

This library allows you to compute uniform, non-uniform, rational, non-rational set-based B-splines. The curve parameter and control points can be real, interval, or affine.

# Building from source

supported and tested on Ubuntu22.04

Install prerequisites:

```apt-get install cmake libgmp-dev python2.7 flex bison gcc g++ make pkg-config libfuse2```

Update git submodules:

```git submodule update --init --recursive```

Install SB2L:

```
mkdir build && cd build
cmake ..
make
make install
```

# Examples

Some examples are available in the ```examples``` folder

```
cd examples
mkdir build && cd build
cmake ..
make
```

# Visualization

Used the [Vibes](https://github.com/ENSTABretagneRobotics/VIBES) library. First launch the ```VIBes-0.2.3-linux.AppImage``` and then an examples.

```
cd VIBes-viewer
./VIBes-0.2.3-linux.AppImage
```