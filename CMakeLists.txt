cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12.0)


project (OpenVolumeMesh
    VERSION 2.1.0
    LANGUAGES CXX
    HOMEPAGE_URL https://openvolumemesh.org/
    )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
    set(OVM_STANDALONE_BUILD TRUE)
else()
    set(OVM_STANDALONE_BUILD FALSE)
endif()


#message("OVM_STANDALONE_BUILD: ${OVM_STANDALONE_BUILD}")
set(OVM_ENABLE_UNITTESTS ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh unit tests in 'all' target")
set(OVM_ENABLE_EXAMPLES ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh examples in 'all' target")
set(OVM_BUILD_DOCUMENTATION ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh documentation in 'all' target")

if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    if (WIN32)
        set (ACG_PROJECT_LIBDIR "lib")
        set (ACG_PROJECT_BINDIR ".")
    else ()
        set (ACG_PROJECT_LIBDIR "lib")
        set (ACG_PROJECT_BINDIR "bin")
    endif ()

    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/${ACG_PROJECT_BINDIR}")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/${ACG_PROJECT_LIBDIR}")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/${ACG_PROJECT_LIBDIR}")
endif()

set (CMAKE_DEBUG_POSTFIX "d")

enable_testing()
add_subdirectory (src)

if (OVM_ENABLE_EXAMPLES)
    add_subdirectory (examples)
endif()

if (OVM_BUILD_DOCUMENTATION)
    add_subdirectory (documentation)
endif()

#export(PACKAGE OpenVolumeMesh)

list(REMOVE_AT CMAKE_MODULE_PATH -1)
