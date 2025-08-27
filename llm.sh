#!/bin/bash

set -x

# Check if debug mode is enabled
DEBUG_MODE=${1:-false}

if [ "$DEBUG_MODE" = "debug" ] || [ "$DEBUG_MODE" = "gdb" ]; then
    echo "Building in Debug mode for gdb debugging..."
    ./build.sh Debug
    BUILD_DIR="cmake-build-debug"
else
    echo "Building in Release mode..."
    ./build.sh Release
    BUILD_DIR="cmake-build-release"
fi

/usr/bin/cmake --build ./${BUILD_DIR} --target demo_llm_run -- -j 6

run_file=./${BUILD_DIR}/src/demo_llm_run

if [ "$DEBUG_MODE" = "debug" ] || [ "$DEBUG_MODE" = "gdb" ]; then
    echo "Running with gdb debugger..."
    gdb --args ${run_file}
else
    echo "Running normally..."
    ${run_file}
fi 