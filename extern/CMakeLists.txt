#libdivsufsort
set(LIBDIVSUFSORT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/libdivsufsort")
if(NOT EXISTS "${LIBDIVSUFSORT_PATH}/.git")
    message(STATUS "libdivsufsort not found or not initialized. Initializing and updating...")
    execute_process(
            COMMAND git submodule update --init extern/libdivsufsort --recursive
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
            RESULT_VARIABLE GIT_SUBMOD_INIT_RESULT
            )
    if(NOT GIT_SUBMOD_INIT_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --init extern/libdivsufsort failed for libdivsufsort with ${GIT_SUBMOD_INIT_RESULT}")
    endif()
endif()
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")
add_subdirectory(libdivsufsort)

#sdsl
if(BUILD_SDSL)
    set(SDSL_PATH "${CMAKE_CURRENT_SOURCE_DIR}/sdsl")
    if(NOT EXISTS "${SDSL_PATH}/.git")
        message(STATUS "sdsl not found or not initialized. Initializing and updating...")
        execute_process(
                COMMAND git submodule update --init extern/sdsl --recursive
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/..
                RESULT_VARIABLE GIT_SUBMOD_INIT_RESULT
        )
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed for sdsl with ${GIT_SUBMOD_RESULT}")
        endif()
    endif()
    add_subdirectory(sdsl)
endif()


