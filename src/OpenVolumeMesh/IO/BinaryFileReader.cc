#include "BinaryFileReader.hh"

#include <array>

namespace OpenVolumeMesh::IO {


const char* to_string(ReadCompatibility rc)
{
    static const std::array<const char*, 6> strings {
        "Ok",
        "MeshVertexDimensionIncompatible",
        "MeshTopologyIncompatible",
        "MeshHandleIncompatible",
        "FileVersionUnsupported",
        "InvalidFile"
    };
    size_t idx = static_cast<size_t>(rc);
    if (idx >= strings.size())
        return nullptr;
    return strings[idx];
}
const char* to_string(ReadState rs) {
    static const std::array<const char*, 17> strings {
        "Init",
        "HeaderRead",
        "ReadingChunks",
        "Finished",
        "Error",
        "ErrorInvalidMagic",
        "ErrorEndNotReached",
        "ErrorIncompatible",
        "ErrorChunkTooBig",
        "ErrorMissingData",
        "ErrorUnsupportedChunkType",
        "ErrorUnsupportedChunkVersion",
        "ErrorInvalidTopoType",
        "ErrorHandleRange",
        "ErrorInvalidEncoding",
        "ErrorEmptyList",
        "ErrorInvalidChunkSize",
    };

    size_t idx = static_cast<size_t>(rs);
    if (idx >= strings.size())
        return nullptr;
    return strings[idx];
}

} // namespace OpenVolumeMesh::IO
