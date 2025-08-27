cmake_minimum_required(VERSION 3.10)
project(zkCNN)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "-mcmodel=large -O3 -lpthread -pthread")
link_directories(3rd/mcl)

include_directories(src)
include_directories(3rd)
include_directories(3rd/mcl/include)

add_subdirectory(src)
add_subdirectory(3rd/mcl)