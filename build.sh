#!/bin/bash

# Default to debug mode if no argument provided
BUILD_TYPE=${1:-Debug}
BUILD_DIR="cmake-build-${BUILD_TYPE,,}"  # Convert to lowercase

echo "Building in ${BUILD_TYPE} mode..."

mkdir -p ${BUILD_DIR}
cd ${BUILD_DIR}
cmake -DCMAKE_BUILD_TYPE=${BUILD_TYPE} -G "CodeBlocks - Unix Makefiles" ..
make
cd ..

#if [ ! -d "./data" ]
#then
#    tar -xzvf data.tar.gz
#fi
#cd script
