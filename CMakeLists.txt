cmake_minimum_required(VERSION 3.20)
project(gin C CXX)
set(PROJECT_VENDOR "Unsal Ozturk")
set(PROJECT_CONTACT "uensal.oeztuerk@gmail.com")
set(PROJECT_URL "https://github.com/uensalo/gin")
set(PROJECT_DESCRIPTION "GIN-TONIC: Non-hierarchical full-text indexing for graphs")
set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_C_STANDARD 99)

set(BUILD_OPENMP             OFF)
set(BUILD_MARCH_NATIVE_FLAG  OFF)
set(BUILD_BUILTIN_POPCOUNT   OFF)
set(BUILD_DEBUG              OFF)
set(BUILD_DNA_FMI            ON)
set(BUILD_SDSL               OFF)
set(BUILD_ORACLE             OFF)

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
    ${PROJECT_SOURCE_DIR}/gin_bitstream.c
    ${PROJECT_SOURCE_DIR}/gin_table.c
    ${PROJECT_SOURCE_DIR}/gin_tree.c
    ${PROJECT_SOURCE_DIR}/gin_vector.c
    ${PROJECT_SOURCE_DIR}/gin_string.c
    ${PROJECT_SOURCE_DIR}/gin_graph.c
    ${PROJECT_SOURCE_DIR}/gin_min_heap.c
    ${PROJECT_SOURCE_DIR}/gin_fmi.c
    ${PROJECT_SOURCE_DIR}/gin_gin.c
    ${PROJECT_SOURCE_DIR}/gin_interval_merge_tree.c
    ${PROJECT_SOURCE_DIR}/gin_oimt.c
    ${PROJECT_SOURCE_DIR}/gin_constraint_set.c
    ${PROJECT_SOURCE_DIR}/gin_annealing.c
    ${PROJECT_SOURCE_DIR}/gin_encoded_graph.c
)
set(H_FILES
    ${PROJECT_HEADER_DIR}/gin_common.h
    ${PROJECT_HEADER_DIR}/gin_table.h
    ${PROJECT_HEADER_DIR}/gin_bitstream.h
    ${PROJECT_HEADER_DIR}/gin_tree.h
    ${PROJECT_HEADER_DIR}/gin_vector.h
    ${PROJECT_HEADER_DIR}/gin_string.h
    ${PROJECT_HEADER_DIR}/gin_graph.h
    ${PROJECT_HEADER_DIR}/gin_min_heap.h
    ${PROJECT_HEADER_DIR}/gin_fmi.h
    ${PROJECT_HEADER_DIR}/gin_gin.h
    ${PROJECT_HEADER_DIR}/gin_interval_merge_tree.h
    ${PROJECT_HEADER_DIR}/gin_oimt.h
    ${PROJECT_HEADER_DIR}/gin_constraint_set.h
    ${PROJECT_HEADER_DIR}/gin_annealing.h
    ${PROJECT_HEADER_DIR}/gin_encoded_graph.h
    ${PROJECT_HEADER_DIR}/rgfa_parser.h
    ${PROJECT_HEADER_DIR}/permutation_parser.h
    ${PROJECT_HEADER_DIR}/ging_parser.h
)

################################################################################
# option specific behaviour
################################################################################

if(BUILD_DEBUG)
    add_definitions(-DGIN_DEBUG)
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

if(BUILD_DNA_FMI)
    list(APPEND C_FILES ${PROJECT_SOURCE_DIR}/gin_dna_fmi.c)
    list(APPEND H_FILES ${PROJECT_HEADER_DIR}/gin_dna_fmi.h)
elseif(BUILD_SDSL)
    list(APPEND C_FILES ${PROJECT_WRAPPER_DIR}/sdsl_wrapper.cpp)
    list(APPEND H_FILES ${PROJECT_WRAPPER_DIR}/sdsl_wrapper.h)
endif()

if(BUILD_OPENMP)
    find_package(OpenMP)
    if (OPENMP_FOUND)
        add_definitions(-DGIN_OMP)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

if(BUILD_BUILTIN_POPCOUNT)
    add_definitions(-DGIN_BUILTIN_POPCOUNT)
endif()

# build dependencies first
add_subdirectory(${PROJECT_EXTERN_DIR}/)

#add_library(gin_lib SHARED ${C_FILES} ${H_FILES})
add_library(gin_lib_static STATIC ${C_FILES} ${H_FILES})

# Specify where the project should look for headers
#target_include_directories(
#    gin_lib PUBLIC ${PROJECT_HEADER_DIR}
#    gin_lib PUBLIC ${CMAKE_BINARY_DIR}/extern/libdivsufsort/include
#)
target_include_directories(
    gin_lib_static PUBLIC ${PROJECT_HEADER_DIR}
    gin_lib_static PUBLIC ${CMAKE_BINARY_DIR}/extern/libdivsufsort-gin/include
)

################################################################################
# dependencies
################################################################################

# divsufsort
#target_link_libraries(gin_lib PUBLIC divsufsort)
target_link_libraries(gin_lib_static PUBLIC divsufsort-gin)
#target_link_libraries(gin_lib PUBLIC divsufsort64)
target_link_libraries(gin_lib_static PUBLIC divsufsort64-gin)

# math
#target_link_libraries(gin_lib PUBLIC m)
target_link_libraries(gin_lib_static PUBLIC m)

#sdsl
if(BUILD_DNA_FMI)
    add_definitions(-DGIN_DNA_FMI)
elseif(BUILD_SDSL)
    #target_include_directories(
    #    gin_lib PUBLIC ${PROJECT_EXTERN_DIR}/sdsl/include
    #)
    target_include_directories(
            gin_lib_static PUBLIC ${PROJECT_EXTERN_DIR}/sdsl/include
    )
    #target_link_libraries(gin_lib PUBLIC sdsl)
    target_link_libraries(gin_lib_static PUBLIC sdsl)
    add_definitions(-DGIN_SDSL)
endif()

#oracle
if(BUILD_ORACLE)
    add_definitions(-DGIN_ORACLE)
endif()

################################################################################
# executables
################################################################################
add_executable(gin gin.c)
target_link_libraries(gin gin_lib_static)
