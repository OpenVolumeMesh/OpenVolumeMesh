SET(SOURCE_FILES
    unittests_smartptr.cc
    unittests_basics.cc
    unittests_files.cc
    unittests_common.cc
    unittests.cc
    unittests_iterators.cc
    unittests_properties.cc)

find_package(OpenVolumeMesh REQUIRED)

enable_testing()
find_package(GTest)

if(GTEST_FOUND)
    include(GoogleTest)

    add_executable(unittests ${SOURCE_FILES})
    target_link_libraries(unittests
        OpenVolumeMesh::OpenVolumeMesh
        GTest::GTest
        GTest::Main
        )

    gtest_add_tests(TARGET unittests
        EXTRA_ARGS "--gtest_color=yes;--gtest_output=xml"
        WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/Unittests"
        )

    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/TestFiles/"
        DESTINATION "${CMAKE_BINARY_DIR}/Unittests"
        )

    set_target_properties(unittests PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/Unittests)

    #[[=
    if ( WIN32 )
            add_custom_command (TARGET unittests POST_BUILD
                                COMMAND ${CMAKE_COMMAND} -E
                                    copy_if_different
                                    ${CMAKE_CURRENT_BINARY_DIR}/${CMAKE_CFG_INTDIR}/unittests.exe
                                    ${CMAKE_BINARY_DIR}/Unittests/unittests.exe)
    endif()
    =]]
endif()
