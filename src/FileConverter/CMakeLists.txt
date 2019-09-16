find_package(Boost 1.43.0)
if(Boost_FOUND)

  add_executable(file_converter ovm_converter.cc)
  set_target_properties (file_converter PROPERTIES
      CXX_STANDARD ${OVM_CXX_STANDARD}
      CXX_STANDARD_REQUIRED YES
      CXX_EXTENSIONS NO
    )
    target_link_libraries(file_converter OpenVolumeMesh::OpenVolumeMesh Boost::boost)

else()
  message(STATUS "Boost not found. Disabling file_converter.")
endif()
