cmake_minimum_required(VERSION 3.20)
project(fmd C CXX)
set(PROJECT_VENDOR "Unsal Ozturk")
set(PROJECT_CONTACT "uensal.oeztuerk@gmail.com")
set(PROJECT_URL "https://github.com/uensalo/fmd")
set(PROJECT_DESCRIPTION "FMD: Sequence domain full-text graph indexing")
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_C_STANDARD 99)

set(BUILD_OPENMP             OFF)
set(BUILD_MARCH_NATIVE_FLAG  OFF)
set(BUILD_BUILTIN_POPCOUNT   OFF)
set(BUILD_DEBUG              OFF)
set(BUILD_SDSL               OFF)
set(BUILD_ORACLE             OFF)
set(BUILD_TESTS              OFF)

set(PROJECT_ROOT_DIR    ${CMAKE_SOURCE_DIR})
set(PROJECT_SOURCE_DIR  ${PROJECT_ROOT_DIR}/src)
set(PROJECT_HEADER_DIR  ${PROJECT_ROOT_DIR}/include)
set(PROJECT_EXTERN_DIR  ${PROJECT_ROOT_DIR}/extern)
set(PROJECT_TESTS_DIR   ${PROJECT_ROOT_DIR}/tests)
set(PROJECT_WRAPPER_DIR ${PROJECT_ROOT_DIR}/wrapper)
################################################################################
# library file definitions
################################################################################
set(C_FILES
    ${PROJECT_SOURCE_DIR}/fmd_bitstream.c
    ${PROJECT_SOURCE_DIR}/fmd_table.c
    ${PROJECT_SOURCE_DIR}/fmd_tree.c
    ${PROJECT_SOURCE_DIR}/fmd_vector.c
    ${PROJECT_SOURCE_DIR}/fmd_string.c
    ${PROJECT_SOURCE_DIR}/fmd_graph.c
    ${PROJECT_SOURCE_DIR}/fmd_min_heap.c
    ${PROJECT_SOURCE_DIR}/fmd_fmi.c
    ${PROJECT_SOURCE_DIR}/fmd_fmd.c
    ${PROJECT_SOURCE_DIR}/fmd_interval_merge_tree.c
    ${PROJECT_SOURCE_DIR}/fmd_oimt.c
    ${PROJECT_SOURCE_DIR}/fmd_constraint_set.c
    ${PROJECT_SOURCE_DIR}/fmd_annealing.c
    ${PROJECT_SOURCE_DIR}/fmd_encoded_graph.c
)
set(H_FILES
    ${PROJECT_HEADER_DIR}/fmd_common.h
    ${PROJECT_HEADER_DIR}/fmd_table.h
    ${PROJECT_HEADER_DIR}/fmd_bitstream.h
    ${PROJECT_HEADER_DIR}/fmd_tree.h
    ${PROJECT_HEADER_DIR}/fmd_vector.h
    ${PROJECT_HEADER_DIR}/fmd_string.h
    ${PROJECT_HEADER_DIR}/fmd_graph.h
    ${PROJECT_HEADER_DIR}/fmd_min_heap.h
    ${PROJECT_HEADER_DIR}/fmd_fmi.h
    ${PROJECT_HEADER_DIR}/fmd_fmd.h
    ${PROJECT_HEADER_DIR}/fmd_interval_merge_tree.h
    ${PROJECT_HEADER_DIR}/fmd_oimt.h
    ${PROJECT_HEADER_DIR}/fmd_constraint_set.h
    ${PROJECT_HEADER_DIR}/fmd_annealing.h
    ${PROJECT_HEADER_DIR}/fmd_encoded_graph.h
    ${PROJECT_HEADER_DIR}/rgfa_parser.h
    ${PROJECT_HEADER_DIR}/permutation_parser.h
    ${PROJECT_HEADER_DIR}/fmdg_parser.h
)

################################################################################
# option specific behaviour
################################################################################

if(BUILD_DEBUG)
    add_definitions(-DFMD_DEBUG)
    add_link_options(-fsanitize=address)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    set(BUILD_OPENMP OFF)
    set(CMAKE_BUILD_TYPE "Debug")
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fPIC")
    if(BUILD_MARCH_NATIVE_FLAG)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    endif()
endif()

if(BUILD_SDSL)
    list(APPEND C_FILES ${PROJECT_WRAPPER_DIR}/sdsl_wrapper.cpp)
    list(APPEND H_FILES ${PROJECT_WRAPPER_DIR}/sdsl_wrapper.h)
endif()

if(BUILD_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        add_definitions(-DFMD_OMP)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        #set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

if(BUILD_BUILTIN_POPCOUNT)
    add_definitions(-DFMD_BUILTIN_POPCOUNT)
endif()

# build dependencies first
add_subdirectory(${PROJECT_EXTERN_DIR}/)

#add_library(fmd_lib SHARED ${C_FILES} ${H_FILES})
add_library(fmd_lib_static STATIC ${C_FILES} ${H_FILES})

# Specify where the project should look for headers
#target_include_directories(
#    fmd_lib PUBLIC ${PROJECT_HEADER_DIR}
#    fmd_lib PUBLIC ${CMAKE_BINARY_DIR}/extern/libdivsufsort/include
#)
target_include_directories(
    fmd_lib_static PUBLIC ${PROJECT_HEADER_DIR}
    fmd_lib_static PUBLIC ${CMAKE_BINARY_DIR}/extern/libdivsufsort/include
)

################################################################################
# dependencies
################################################################################

# divsufsort
#target_link_libraries(fmd_lib PUBLIC divsufsort)
target_link_libraries(fmd_lib_static PUBLIC divsufsortf)
#target_link_libraries(fmd_lib PUBLIC divsufsort64)
target_link_libraries(fmd_lib_static PUBLIC divsufsort64f)

# math
#target_link_libraries(fmd_lib PUBLIC m)
target_link_libraries(fmd_lib_static PUBLIC m)

#sdsl
if(BUILD_SDSL)
    #target_include_directories(
    #    fmd_lib PUBLIC ${PROJECT_EXTERN_DIR}/sdsl/include
    #)
    target_include_directories(
            fmd_lib_static PUBLIC ${PROJECT_EXTERN_DIR}/sdsl/include
    )
    #target_link_libraries(fmd_lib PUBLIC sdsl)
    target_link_libraries(fmd_lib_static PUBLIC sdsl)
    add_definitions(-DFMD_SDSL)
endif()

#oracle
if(BUILD_ORACLE)
    add_definitions(-DFMD_ORACLE)
endif()

################################################################################
# executables
################################################################################
add_executable(fmd fmd.c)
target_link_libraries(fmd fmd_lib_static)

if(BUILD_TESTS)
    file(GLOB unit_tests "${PROJECT_TESTS_DIR}/*.c")
    foreach(unit_test ${unit_tests})
        get_filename_component(TEST_NAME ${unit_test} NAME_WE [CACHE])
        add_executable(${TEST_NAME} ${unit_test})
        target_link_libraries(${TEST_NAME} fmd_lib_static)
    endforeach()
endif()
