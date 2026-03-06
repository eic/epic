#!/bin/bash

### Compilation
# To configure, build, and install the geometry (to the `install` directory), use the following commands:
cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j20
cmake --install build

# To load the geometry, you can use the scripts in the `install` directory:   
source install/bin/thisepic.sh 
