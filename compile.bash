#!/bin/bash

cmake -B build -S . -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j20
cmake --install build
source install/setup.sh
