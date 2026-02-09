FBGA 3D
=========

FBGA 3D is an implementation of Forward backward method with generic complex GGGV diagram for 3D spatial paths written is a C++. The library computes lap-time prediction and velocity profile optimization of vehicle dynamics consistent with the provided GGGV diagram. It is designed to be used in the context of autonomous driving and vehicle dynamics research. 

## Requirements

### General

- git
- cmake
- make

### Third party libraries

Some third party libraries are needed to compile the project.

need to run the following command to install the third party libraries:

```{shell}
./third_party.sh
```

Most of the library are needed just for the tests, so if you don't need them you can skip the installation.

  
## Building

To install the project you need to run the following commands:

```{shell}
./build.sh
```

Default will build the project in release mode, if you want to build in debug mode you can run:

```{shell}
./build.sh -debug
```

For further information you can run:

```{shell}
./build.sh -h
```

### Alternative build

If you want to install the project in a different directory you can run:

```{shell}
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE="$BUILD_TYPE" ..
make -j
```
