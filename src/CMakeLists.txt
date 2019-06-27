include (ACGCommon)

set(OVM_ENABLE_DEPRECATED_APIS ON CACHE BOOL "Support deprecated OpenVolumeMesh APIs")

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/Config/DeprecationConfig.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/Config/DeprecationConfig.hh"
)

include_directories (
  ..
  ${CMAKE_CURRENT_SOURCE_DIR}
)

# source code directories
set (directories 
  .
  OpenVolumeMesh/Attribs
  OpenVolumeMesh/Core
  OpenVolumeMesh/Config
  OpenVolumeMesh/FileManager
  OpenVolumeMesh/Geometry
  OpenVolumeMesh/Mesh
)

# collect all header and source files
acg_append_files (headers "*.hh" ${directories})
acg_append_files (sources "*.cc" ${directories})

# Disable Library installation when not building OpenVolumeMesh on its own but as part of another project!
if ( NOT ${PROJECT_NAME} MATCHES "OpenVolumeMesh")
  set(ACG_NO_LIBRARY_INSTALL true)
endif()

if (WIN32)
    # OpenVolumeMesh has no dll exports so we have to build a static library on windows
    acg_add_library (OpenVolumeMesh STATIC ${sources} ${headers})
else ()
    acg_add_library (OpenVolumeMesh SHAREDANDSTATIC ${sources} ${headers})
    set_target_properties (OpenVolumeMesh PROPERTIES VERSION ${OPENVOLUMEMESH_VERSION_MAJOR}.${OPENVOLUMEMESH_VERSION_MINOR}
                                          SOVERSION ${OPENVOLUMEMESH_VERSION_MAJOR}.${OPENVOLUMEMESH_VERSION_MINOR} )
endif ()
foreach(TARGET_NAME OpenVolumeMesh OpenVolumeMeshStatic)
    if (TARGET ${TARGET_NAME})
        target_include_directories(${TARGET_NAME} PUBLIC "${CMAKE_CURRENT_BINARY_DIR}/Config/")
    endif()
endforeach(TARGET_NAME)

# Only install if the project name matches OpenVolumeMesh.
if (NOT APPLE AND ${PROJECT_NAME} MATCHES "OpenVolumeMesh")

# Install Header Files)
install(DIRECTORY . 
        DESTINATION include
        FILES_MATCHING 
        PATTERN "*.hh"
        PATTERN "Unittests" EXCLUDE
        PATTERN "FileConverter" EXCLUDE
        PATTERN "CVS" EXCLUDE
        PATTERN ".svn" EXCLUDE
        PATTERN "tmp" EXCLUDE
        PATTERN "Templates" EXCLUDE
        PATTERN "Debian*" EXCLUDE)
endif ()

# install config files
install( FILES ${CMAKE_CURRENT_BINARY_DIR}/Config/DeprecationConfig.hh
           DESTINATION include/OpenVolumeMesh/Config )

# Only build unittests and file converter
# if not built as external library
if(${PROJECT_NAME} MATCHES "OpenVolumeMesh")
    # Add unittests target
    add_subdirectory(Unittests)
    add_subdirectory(FileConverter)
endif()
