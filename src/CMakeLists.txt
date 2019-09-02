SET(HEADER_FILES
    OpenVolumeMesh/System/FunctionalInclude.hh
    OpenVolumeMesh/System/MemoryInclude.hh
    OpenVolumeMesh/System/Deprecation.hh
    OpenVolumeMesh/Attribs/TexCoordAttrib.hh
    OpenVolumeMesh/Attribs/OpenVolumeMeshStatus.hh
    OpenVolumeMesh/Attribs/ColorAttribT_impl.hh
    OpenVolumeMesh/Attribs/ColorAttrib.hh
    OpenVolumeMesh/Attribs/InterfaceAttrib.hh
    OpenVolumeMesh/Attribs/NormalAttrib.hh
    OpenVolumeMesh/Attribs/NormalAttribT_impl.hh
    OpenVolumeMesh/Attribs/TexCoordAttribT_impl.hh
    OpenVolumeMesh/Attribs/StatusAttribT_impl.hh
    OpenVolumeMesh/Attribs/StatusAttrib.hh
    OpenVolumeMesh/Geometry/Vector11T.hh
    OpenVolumeMesh/Geometry/VectorT_inc.hh
    OpenVolumeMesh/Geometry/VectorT.hh
    OpenVolumeMesh/FileManager/FileManagerT_impl.hh
    OpenVolumeMesh/FileManager/FileManager.hh
    OpenVolumeMesh/Core/BaseEntities.hh
    OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh
    OpenVolumeMesh/Core/BaseProperty.hh
    OpenVolumeMesh/Core/ResourceManager.hh
    OpenVolumeMesh/Core/PropertyHandles.hh
    OpenVolumeMesh/Core/TypeName.hh
    OpenVolumeMesh/Core/PropertyPtrT_impl.hh
    OpenVolumeMesh/Core/PropertyDefinesT_impl.hh
    OpenVolumeMesh/Core/GeometryKernel.hh
    OpenVolumeMesh/Core/TopologyKernel.hh
    OpenVolumeMesh/Core/ResourceManagerT_impl.hh
    OpenVolumeMesh/Core/PropertyDefines.hh
    OpenVolumeMesh/Core/Serializers.hh
    OpenVolumeMesh/Core/OpenVolumeMeshBaseProperty.hh
    OpenVolumeMesh/Core/ForwardDeclarations.hh
    OpenVolumeMesh/Core/PropertyPtr.hh
    OpenVolumeMesh/Core/Entities.hh
    OpenVolumeMesh/Core/SerializersT_impl.hh
    OpenVolumeMesh/Core/OpenVolumeMeshProperty.hh
    OpenVolumeMesh/Core/Iterators.hh
    OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.hh
    OpenVolumeMesh/Mesh/HexahedralMesh.hh
    OpenVolumeMesh/Mesh/HexahedralMeshIterators.hh
    OpenVolumeMesh/Mesh/TetrahedralMesh.hh
    OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.hh
    OpenVolumeMesh/Mesh/TetrahedralMeshIterators.hh
    OpenVolumeMesh/Mesh/PolyhedralMesh.hh
    OpenVolumeMesh/Mesh/TetrahedralGeometryKernel.hh
    Unittests/unittests_common.hh)

SET(SOURCE_FILES
    OpenVolumeMesh/Attribs/InterfaceAttrib.cc
    OpenVolumeMesh/Attribs/OpenVolumeMeshStatus.cc
    OpenVolumeMesh/Attribs/StatusAttrib.cc
    OpenVolumeMesh/Geometry/Vector.cc
    OpenVolumeMesh/FileManager/FileManager.cc
    OpenVolumeMesh/Core/PropertyDefines.cc
    OpenVolumeMesh/Core/OpenVolumeMeshHandle.cc
    OpenVolumeMesh/Core/ResourceManager.cc
    OpenVolumeMesh/Core/BaseEntities.cc
    OpenVolumeMesh/Core/BaseProperty.cc
    OpenVolumeMesh/Core/TopologyKernel.cc
    OpenVolumeMesh/Core/Serializers.cc
    OpenVolumeMesh/Core/Iterators.cc
    OpenVolumeMesh/Core/OpenVolumeMeshBaseProperty.cc
    OpenVolumeMesh/Mesh/TetrahedralMeshIterators.cc
    OpenVolumeMesh/Mesh/HexahedralMeshIterators.cc
    OpenVolumeMesh/Mesh/TetrahedralMeshTopologyKernel.cc
    OpenVolumeMesh/Mesh/HexahedralMeshTopologyKernel.cc
)

set(OVM_ENABLE_DEPRECATED_APIS ON CACHE BOOL "Support deprecated OpenVolumeMesh APIs")

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/Config/DeprecationConfig.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config/DeprecationConfig.hh"
)

configure_file (
    "${CMAKE_CURRENT_SOURCE_DIR}/OpenVolumeMesh/Config/Version.hh.in"
    "${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config/Version.hh"
)


# If we are in a subdirectory of another project that has a ACGCommon module
# (e.g. OpenFlipper), we will use acg_add_library; otherwise behave "normally".
include(ACGCommon OPTIONAL RESULT_VARIABLE IN_ACG_BUILD)
if(IN_ACG_BUILD)
    message("OVM is being built as part of an ACG project, using acg_add_library!")
    if (WIN32)
        acg_add_library(OpenVolumeMesh STATIC ${SOURCE_FILES})
    else()
        acg_add_library(OpenVolumeMesh SHARED ${SOURCE_FILES})
    endif()
else(IN_ACG_BUILD)
    if (WIN32)
        # Windows .dll build is not supported yet.
        # TODO: turn this on once we have OVM_EXPORT macros in all required places,
        # cf CXX_VISIBILITY_PRESET:
        add_library(OpenVolumeMesh STATIC ${SOURCE_FILES})
    else()
        option(BUILD_SHARED_LIBS "Build shared library (*.dll, *.so, *.dylib) instead of static library (*.a, *.lib)" ON)
        add_library(OpenVolumeMesh ${SOURCE_FILES})
    endif()
endif(IN_ACG_BUILD)

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
    CXX_STANDARD ${OVM_CXX_VERSION}
    CXX_STANDARD_REQUIRED YES
    CXX_EXTENSIONS NO
    # TODO: turn this on once we have OVM_EXPORT macros in all required places:
    #CXX_VISIBILITY_PRESET hidden
    POSITION_INDEPENDENT_CODE ON
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

target_compile_definitions(OpenVolumeMesh
    PUBLIC
        "INCLUDE_TEMPLATES"
        "$<$<CXX_COMPILER_ID:MSVC>:NOMINMAX>"
        "$<$<CXX_COMPILER_ID:MSVC>:_USE_MATH_DEFINES>"
    )

target_include_directories(OpenVolumeMesh
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/OpenVolumeMesh/Config>
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
    COMPATIBILITY SameMajorVersion
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
export(PACKAGE OpenVolumeMesh)

if (OVM_ENABLE_APPLICATIONS)
    add_subdirectory(FileConverter)
else()
    add_subdirectory(FileConverter EXCLUDE_FROM_ALL)
endif()

if (OVM_ENABLE_UNITTESTS)
    add_subdirectory(Unittests)
else()
    add_subdirectory(Unittests EXCLUDE_FROM_ALL)
endif()

