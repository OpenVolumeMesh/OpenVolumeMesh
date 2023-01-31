

SET(SOURCE_FILES
    OpenVolumeMesh/Attribs/InterfaceAttrib.cc
    OpenVolumeMesh/Attribs/OpenVolumeMeshStatus.cc
    OpenVolumeMesh/Attribs/StatusAttrib.cc
    OpenVolumeMesh/FileManager/FileManager.cc
    OpenVolumeMesh/FileManager/TypeNames.cc
    OpenVolumeMesh/FileManager/Serializers.cc
    OpenVolumeMesh/Core/Handles.cc
    OpenVolumeMesh/Core/ResourceManager.cc
    OpenVolumeMesh/Core/BaseEntities.cc
    OpenVolumeMesh/Core/TopologyKernel.cc
    OpenVolumeMesh/Core/Iterators.cc
    OpenVolumeMesh/Core/detail/internal_type_name.cc
    OpenVolumeMesh/Core/Properties/PropertyStorageBase.cc
    OpenVolumeMesh/IO/enums.cc
    OpenVolumeMesh/IO/PropertyCodecs.cc
    OpenVolumeMesh/IO/detail/BinaryIStream.cc
    OpenVolumeMesh/IO/detail/BinaryFileReader.cc
    OpenVolumeMesh/IO/detail/BinaryFileWriter.cc
    OpenVolumeMesh/IO/detail/GeometryWriter.cc
    OpenVolumeMesh/IO/detail/GeometryReader.cc
    OpenVolumeMesh/IO/detail/Decoder.cc
    OpenVolumeMesh/IO/detail/Encoder.cc
    OpenVolumeMesh/IO/detail/ovmb_format.cc
    OpenVolumeMesh/IO/detail/ovmb_codec.cc
    OpenVolumeMesh/IO/detail/WriteBuffer.cc
    OpenVolumeMesh/Mesh/TetrahedralMeshIterators.cc
    OpenVolumeMesh/Mesh/HexahedralMeshIterators.cc
    OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.cc
    OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.cc
    OpenVolumeMesh/Unstable/Topology/TetTopology.cc
    OpenVolumeMesh/Unstable/Topology/TriangleTopology.cc
)

# No Qt in any of our headers:
file(GLOB_RECURSE HEADER_FILES OpenVolumeMesh/*.hh)
set_property(SOURCE ${HEADER_FILES} PROPERTY SKIP_AUTOMOC ON)
set_property(SOURCE ${HEADER_FILES} PROPERTY SKIP_AUTOUIC ON)
set_property(SOURCE ${HEADER_FILES} PROPERTY SKIP_AUTOGEN ON)

set(OVM_ENABLE_DEPRECATED_APIS OFF CACHE BOOL "Support deprecated OpenVolumeMesh APIs")

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/Config/DeprecationConfig.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config/DeprecationConfig.hh"
)

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/Config/Version.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config/Version.hh"
)


# If we are in a subdirectory of another project that has a VCICommon module
# (e.g. OpenFlipper), we will use vci_add_library; otherwise behave "normally".
include(VCICommon OPTIONAL RESULT_VARIABLE IN_VCI_BUILD)
if(IN_VCI_BUILD)
    message("OVM is being built as part of an VCI project, using vci_add_library!")
    vci_add_library(OpenVolumeMesh ${SOURCE_FILES})
else(IN_VCI_BUILD)
    add_library(OpenVolumeMesh ${SOURCE_FILES})
endif(IN_VCI_BUILD)

add_library(OpenVolumeMesh::OpenVolumeMesh ALIAS OpenVolumeMesh)

include(GenerateExportHeader)
generate_export_header(OpenVolumeMesh
    BASE_NAME OVM
    EXPORT_FILE_NAME "OpenVolumeMesh/Config/Export.hh"
    DEPRECATED_MACRO_NAME "CMAKE_OVM_DEPRECATED"
    NO_DEPRECATED_MACRO_NAME "CMAKE_OVM_NO_DEPRECATED"
)

set_target_properties (OpenVolumeMesh PROPERTIES
    VERSION ${OpenVolumeMesh_VERSION_MAJOR}.${OpenVolumeMesh_VERSION_MINOR}
    SOVERSION ${OpenVolumeMesh_VERSION_MAJOR}.${OpenVolumeMesh_VERSION_MINOR}
    POSITION_INDEPENDENT_CODE ON
    AUTOMOC OFF
    AUTOUIC OFF
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN TRUE
    )


target_compile_options(OpenVolumeMesh PUBLIC
    "$<$<CXX_COMPILER_ID:MSVC>:/bigobj>"
    "$<$<CXX_COMPILER_ID:MSVC>:/Zc:__cplusplus>" # Enable updated __cplusplus macro
    PRIVATE
    "$<$<CXX_COMPILER_ID:GNU>:-Wall>"
    "$<$<CXX_COMPILER_ID:Clang>:-Wall;-Wno-c++98-compat;-Wno-c++98-compat-pedantic;-Wno-weak-vtables;-Wno-padded;-Wno-sign-conversion;-Wno-missing-prototypes>"
    "$<$<CXX_COMPILER_ID:MSVC>:/W1>"
    # disable MSVC inline warnings that are enabled with /Ob2, which is in the default cmake release cxxflags:
    "$<$<CXX_COMPILER_ID:MSVC>:/wd4710>"
    "$<$<CXX_COMPILER_ID:MSVC>:/wd4711>"
    )

target_include_directories(OpenVolumeMesh
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

include(CMakePackageConfigHelpers)
include(GNUInstallDirs)


set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/OpenVolumeMesh)

configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMeshConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMeshConfig.cmake"
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
    )

write_basic_package_version_file(
    OpenVolumeMeshConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMinorVersion
    )

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMeshConfigVersion.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMeshConfig.cmake"
    DESTINATION ${INSTALL_CONFIGDIR})


install(DIRECTORY
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config"
    DESTINATION include/OpenVolumeMesh)

# Install Header Files
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh
    DESTINATION include
    FILES_MATCHING PATTERN "*.hh"
    )

install(TARGETS OpenVolumeMesh
    EXPORT OpenVolumeMeshTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(EXPORT OpenVolumeMeshTargets
    FILE OpenVolumeMeshTargets.cmake
    NAMESPACE OpenVolumeMesh::
    DESTINATION ${INSTALL_CONFIGDIR}
)

export(EXPORT OpenVolumeMeshTargets
    NAMESPACE OpenVolumeMesh::)
#export(PACKAGE OpenVolumeMesh)

if (OVM_ENABLE_APPLICATIONS)
    add_subdirectory(FileConverter)
endif()

if (OVM_ENABLE_UNITTESTS)
    add_subdirectory(Unittests)
endif()
