cmake_minimum_required(VERSION 3.20)
project(fmd C)

set(CMAKE_C_STANDARD 99)

set(BUILD_OPENMP OFF)
set(DEBUG ON)

if(DEBUG)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fPIC")
endif()

if(BUILD_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        add_definitions(-DLS_OMP)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

add_executable(fmd main.c include/fmd_common.h include/fmd_table.h include/fmd_bitstream.h src/fmd_bitstream.c)
