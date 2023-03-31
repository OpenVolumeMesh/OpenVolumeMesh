SET(SOURCE_FILES
    unittests.cc
    unittests_attribs.cc
    unittests_basics.cc
    unittests_files.cc
    unittests_local_topology.cc
    unittests_common.cc
    unittests_handles.cc
    unittests_swap_entities.cc
    unittests_iterators.cc
    unittests_mesh_copies.cc
    unittests_smart_tagger.cc
    unittests_properties.cc)

if (NOT TARGET OpenVolumeMesh::OpenVolumeMesh)
    find_package(OpenVolumeMesh REQUIRED)
endif()


include(FetchContent)

if (NOT DEFINED FETCHCONTENT_UPDATES_DISCONNECTED)
    set (DEFINED FETCHCONTENT_UPDATES_DISCONNECTED ON)
endif()

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.12.1
)

if(WIN32)
    # avoid linking errors, cf https://stackoverflow.com/questions/12540970/how-to-make-gtest-build-mdd-instead-of-mtd-by-default-using-cmake
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
endif()
FetchContent_MakeAvailable(googletest)



set(BUILD_GMOCK OFF)
include(GoogleTest)
set(INSTALL_GTEST OFF)

set_target_properties(gtest PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY} )

enable_testing()

set(TARGET_NAME "${OVM_TARGET_PREFIX}unittests")
add_executable(${TARGET_NAME} ${SOURCE_FILES})
target_link_libraries(${TARGET_NAME}
    OpenVolumeMesh::OpenVolumeMesh
    GTest::gtest GTest::gtest_main
    )

gtest_add_tests(TARGET ${TARGET_NAME}
    EXTRA_ARGS "--gtest_color=yes;--gtest_output=xml"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Unittests"
    )

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/TestFiles/"
    DESTINATION "${CMAKE_BINARY_DIR}/Unittests"
    )

#set_target_properties(${TARGET_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Unittests )
