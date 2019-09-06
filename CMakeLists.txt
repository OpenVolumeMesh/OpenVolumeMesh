cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
cmake_policy(VERSION 3.12.0)


project (OpenVolumeMesh
    VERSION 2.1.0
    LANGUAGES CXX
    HOMEPAGE_URL https://openvolumemesh.org/
    )

#list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

if( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
    set(OVM_STANDALONE_BUILD TRUE)
else()
    set(OVM_STANDALONE_BUILD FALSE)
endif()

set(OVM_ENABLE_APPLICATIONS ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh applications")
set(OVM_ENABLE_UNITTESTS ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh unit tests")
set(OVM_ENABLE_EXAMPLES ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh examples")
set(OVM_BUILD_DOCUMENTATION ${OVM_STANDALONE_BUILD} CACHE BOOL "Build OpenVolumeMesh documentation")

set(OVM_CXX_VERSION 14 CACHE STRING "C++ version to use")
set_property(CACHE OVM_CXX_VERSION PROPERTY STRINGS 14 17 20)

if (OVM_STANDALONE_BUILD)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/bin")
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/Build/lib")
    set(OVM_TARGET_PREFIX "")
    enable_testing()
else()
    # OVM_TARGET_PREFIX is prepended to targets with generic names,
    # i.e., "doc", "unittests"
    set(OVM_TARGET_PREFIX "OVM_")
endif()

set (CMAKE_DEBUG_POSTFIX "d")

add_subdirectory (src)

if (OVM_ENABLE_EXAMPLES)
    add_subdirectory(examples)
endif()

if (OVM_BUILD_DOCUMENTATION)
   add_subdirectory(documentation)
endif()

#export(PACKAGE OpenVolumeMesh)

#list(REMOVE_AT CMAKE_MODULE_PATH -1)
