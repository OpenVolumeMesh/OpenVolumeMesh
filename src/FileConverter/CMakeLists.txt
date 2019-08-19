include (ACGCommon)

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
)

find_package(Boost 1.43.0)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})

  # Add target for first example
  add_executable(file_converter EXCLUDE_FROM_ALL ovm_converter.cc)
  set_target_properties (file_converter PROPERTIES
      CXX_STANDARD ${OVM_CXX_VERSION}
      CXX_STANDARD_REQUIRED YES
      CXX_EXTENSIONS NO
    )
  target_link_libraries(file_converter OpenVolumeMesh)

  if(NOT WIN32)
      set_target_properties(file_converter PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/FileConverter)
  endif()

else()

  message(STATUS "Boost not found. Disabling file_converter.")

endif()
