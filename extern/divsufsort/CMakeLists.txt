cmake_minimum_required(VERSION 3.16.3)
## Project information ##
project(libdivsufsort C)
set(PROJECT_VENDOR "Yuta Mori")
set(PROJECT_CONTACT "yuta.256@gmail.com")
set(PROJECT_URL "https://github.com/y-256/libdivsufsort")
set(PROJECT_DESCRIPTION "A lightweight suffix sorting library")

set(CMAKE_C_STANDARD 99)

add_library(
        divsufsort OBJECT
        include/config.h
        include/divsufsort.h
        include/divsufsort64.h
        include/divsufsort_private.h
        include/lfs.h
        src/divsufsort.c
        src/sssort.c
        src/trsort.c
        src/utils.c
)

add_library(
        divsufsort_static OBJECT
        include/config.h
        include/divsufsort.h
        include/divsufsort64.h
        include/divsufsort_private.h
        include/lfs.h
        src/divsufsort.c
        src/sssort.c
        src/trsort.c
        src/utils.c
)

target_include_directories(
        divsufsort PUBLIC "include"
)
target_include_directories(
        divsufsort_static PUBLIC "include"
)