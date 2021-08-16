add_executable(ovmb_converter ovmb_converter.cc)
target_link_libraries(ovmb_converter OpenVolumeMesh::OpenVolumeMesh)

find_package(Boost 1.74.0 QUIET)
if(Boost_FOUND)

  add_executable(file_converter ovm_converter.cc)
    target_link_libraries(file_converter OpenVolumeMesh::OpenVolumeMesh Boost::boost)

else()
  message(STATUS "Boost not found. Disabling file_converter.")
endif()
