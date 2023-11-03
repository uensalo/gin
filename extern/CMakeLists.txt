find_package(Git QUIET)
if(NOT GIT_FOUND)
    message(FATAL_ERROR "Git not found, please make sure it is installed and in the PATH.")
endif()

function(init_submodule path)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${path}/.git")
        message(STATUS "${path} submodule not found or not initialized. Initializing and updating...")
        execute_process(
                COMMAND ${GIT_EXECUTABLE} submodule update --init -- ${path}
                WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                RESULT_VARIABLE GIT_SUBMOD_RESULT
        )
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init failed for ${path} with ${GIT_SUBMOD_RESULT}")
        endif()
    else()
        message(STATUS "${path} submodule already initialized.")
    endif()
endfunction()

option(UPDATE_SUBMODULES "Update submodules to the latest commit" OFF)
if(UPDATE_SUBMODULES)
    message(STATUS "Updating submodules to the latest commit...")
    execute_process(
            COMMAND ${GIT_EXECUTABLE} submodule update --init --remote
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
            RESULT_VARIABLE GIT_SUBMOD_REM_RESULT
    )
    if(NOT GIT_SUBMOD_REM_RESULT EQUAL "0")
        message(FATAL_ERROR "git submodule update --remote failed with ${GIT_SUBMOD_REM_RESULT}")
    endif()
endif()

init_submodule(libdivsufsort-gin)
init_submodule(sdsl)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort examples")
set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bit mode")
add_subdirectory(libdivsufsort-gin)

if(BUILD_SDSL)
    add_subdirectory(sdsl)
endif()
