cmake_minimum_required(VERSION 3.20)
project(fmd C)

set(CMAKE_C_STANDARD 99)

set(BUILD_OPENMP OFF)
set(DEBUG ON)

if(DEBUG)
    add_link_options(-fsanitize=address)
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
    set(BUILD_OPENMP OFF)
else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -fPIC")
endif()

set(PROJECT_ROOT_DIR   ${CMAKE_SOURCE_DIR})
set(PROJECT_SOURCE_DIR ${PROJECT_ROOT_DIR}/src)
set(PROJECT_HEADER_DIR ${PROJECT_ROOT_DIR}/include)
set(PROJECT_EXTERN_DIR ${PROJECT_ROOT_DIR}/extern)
set(PROJECT_TESTS_DIR  ${PROJECT_ROOT_DIR}/tests)
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
    ${PROJECT_SOURCE_DIR}/fmd_fmd.c
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
    ${PROJECT_HEADER_DIR}/fmd_fmd.h
)
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
################################################################################
# dependencies
################################################################################
add_subdirectory(${PROJECT_EXTERN_DIR}/)

# divsufsort
target_include_directories(
        divsufsort  PUBLIC ${PROJECT_EXTERN_DIR}/divsufsort/include/
)
target_link_libraries(fmd_lib        PUBLIC divsufsort)
target_link_libraries(fmd_lib_static PUBLIC divsufsort_static)
################################################################################
# executables
################################################################################
add_executable(fmd main.c)
target_link_libraries(fmd fmd_lib_static)

file(GLOB unit_tests "${PROJECT_TESTS_DIR}/*.c")
foreach(unit_test ${unit_tests})
    get_filename_component(TEST_NAME ${unit_test} NAME_WE [CACHE])
    add_executable(${TEST_NAME} ${unit_test})
    target_link_libraries(${TEST_NAME} fmd_lib_static)
endforeach()
