find_package(Doxygen REQUIRED dot)
set(TARGET_NAME "${OVM_TARGET_PREFIX}doc")

if (DOXYGEN_FOUND)
    set(DOXYGEN_PROJECT_LOGO "${CMAKE_CURRENT_SOURCE_DIR}/html/logo_vci.png")
    set(DOXYGEN_IMAGE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/images")
    set(DOXYGEN_WARN_IF_UNDOCUMENTED "FALSE")
    set(DOXYGEN_CITE_BIB_FILES "${CMAKE_CURRENT_SOURCE_DIR}/references.bib")

    if (${ACG_PROJECT_DATADIR})
        set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build/${ACG_PROJECT_DATADIR}/Doc/html")
    else()
        set(DOXYGEN_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/Build/share/doc/OpenVolumeMesh")
    endif()

    doxygen_add_docs(${TARGET_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/..")

endif()
