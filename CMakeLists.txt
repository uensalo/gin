cmake_minimum_required(VERSION 3.20)
project(fmd C CXX)

set(CMAKE_C_STANDARD 99)

set(BUILD_OPENMP OFF)
set(BUILD_MARCH_NATIVE_FLAG ON)
set(BUILD_BUILTIN_POPCOUNT ON)

set(BUILD_TESTS OFF)
set(BUILD_DEBUG OFF)

set(BUILD_SDSL ON)

if(BUILD_DEBUG)
    add_link_options(-fsanitize=address)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    set(BUILD_OPENMP OFF)
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fPIC")
    if(BUILD_MARCH_NATIVE_FLAG)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -march=native")
    endif()
endif()

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
    ${PROJECT_SOURCE_DIR}/fmd_interval_merge_tree.c
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
    ${PROJECT_HEADER_DIR}/fmd_interval_merge_tree.h
    ${PROJECT_HEADER_DIR}/fmd_constraint_set.h
    ${PROJECT_HEADER_DIR}/fmd_annealing.h
    ${PROJECT_HEADER_DIR}/fmd_encoded_graph.h
    ${PROJECT_HEADER_DIR}/rgfa_parser.h
    ${PROJECT_HEADER_DIR}/permutation_parser.h
    ${PROJECT_HEADER_DIR}/fmdg_parser.h
)

if(BUILD_SDSL)
    list(APPEND C_FILES ${PROJECT_SOURCE_DIR}/fmd_fmd_sdsl.c)
    list(APPEND C_FILES ${PROJECT_HEADER_DIR}/fmd_fmd_sdsl.h)
    list(APPEND C_FILES ${PROJECT_WRAPPER_DIR}/sdsl_wrapper.cpp)
    list(APPEND C_FILES ${PROJECT_WRAPPER_DIR}/sdsl_wrapper.h)
else()
    list(APPEND C_FILES ${PROJECT_SOURCE_DIR}/fmd_fmd.c)
    list(APPEND C_FILES ${PROJECT_HEADER_DIR}/fmd_fmd.h)
endif()

add_library(fmd_lib        SHARED ${C_FILES} ${H_FILES})
add_library(fmd_lib_static STATIC ${C_FILES} ${H_FILES})
target_include_directories(
        fmd_lib        PUBLIC ${PROJECT_HEADER_DIR}/
)
target_include_directories(
        fmd_lib_static PUBLIC ${PROJECT_HEADER_DIR}/
)
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

################################################################################
# dependencies
################################################################################
add_subdirectory(${PROJECT_EXTERN_DIR}/)

# divsufsort
# Find divsufsort headers
find_path(DIVSUFSORT_INCLUDE_DIR divsufsort.h)
include_directories(${DIVSUFSORT_INCLUDE_DIR})

# Find divsufsort library
find_library(DIVSUFSORT_LIBRARY divsufsort64)
if(NOT DIVSUFSORT_LIBRARY)
    message(FATAL_ERROR "Divsufsort library not found")
endif()
#target_link_libraries(fmd_lib PUBLIC ${DIVSUFSORT_LIBRARY})
target_link_libraries(fmd_lib_static PUBLIC divsufsort64)
target_link_libraries(fmd_lib_static PUBLIC divsufsort)

target_link_libraries(fmd_lib PUBLIC divsufsort64)
target_link_libraries(fmd_lib PUBLIC divsufsort)

# math
target_link_libraries(fmd_lib PUBLIC m)
target_link_libraries(fmd_lib_static PUBLIC m)

#sdsl
if(BUILD_SDSL)
    target_link_libraries(fmd_lib PUBLIC sdsl)
    target_link_libraries(fmd_lib_static PUBLIC sdsl)
    add_definitions(-DFMD_SDSL)
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
