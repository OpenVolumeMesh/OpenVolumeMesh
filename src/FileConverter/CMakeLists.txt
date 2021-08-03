add_executable(ovmb_converter ovmb_converter.cc)
target_link_libraries(ovmb_converter OpenVolumeMesh::OpenVolumeMesh)
set_target_properties(ovmb_converter PROPERTIES CXX_STANDARD 20)

find_package(Boost 1.43.0 QUIET)
if(Boost_FOUND)

  add_executable(file_converter ovm_converter.cc)
    target_link_libraries(file_converter OpenVolumeMesh::OpenVolumeMesh Boost::boost)

else()
  message(STATUS "Boost not found. Disabling file_converter.")
endif()
