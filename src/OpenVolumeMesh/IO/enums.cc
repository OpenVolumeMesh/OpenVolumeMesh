#include <array>
#include <OpenVolumeMesh/IO/enums.hh>
#include <cstddef>

namespace OpenVolumeMesh::IO {

const char* to_string(ReadCompatibility rc)
{
    static const std::array<const char*, 8> strings {
        "Ok",
        "CannotOpenFile",
        "BadStream",
        "MeshVertexDimensionIncompatible",
        "MeshTopologyIncompatible",
        "MeshHandleIncompatible",
        "FileVersionUnsupported",
        "InvalidFile"
    };
    size_t idx = static_cast<size_t>(rc);
    if (idx >= strings.size())
        return nullptr;
    // cppcheck-suppress CastIntegerToAddressAtReturn ; false positive
    return strings[idx];
}
const char* to_string(ReadState rs) {
    static const std::array<const char*, 20> strings {
        "Ok",
        "CannotOpenFile",
        "BadStream",
        "Init",
        "HeaderRead",
        "ReadingChunks",
        "Finished",
        "Error",
        "ErrorInvalidFile",
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
    // cppcheck-suppress CastIntegerToAddressAtReturn ; false positive
    return strings[idx];
}

const char* to_string(ReadResult result) {
    static const std::array<const char*, 6> strings {
        "Ok",
        "OtherError",
        "InvalidFile",
        "CannotOpenFile",
        "BadStream",
        "IncompatibleMesh",
    };

    size_t idx = static_cast<size_t>(result);
    if (idx >= strings.size())
        return nullptr;
    // cppcheck-suppress CastIntegerToAddressAtReturn ; false positive
    return strings[idx];
}
const char* to_string(WriteResult result) {
    static const std::array<const char*, 4> strings {
        "Ok",
        "Error",
        "CannotOpenFile",
        "BadStream",
    };

    size_t idx = static_cast<size_t>(result);
    if (idx >= strings.size())
        return nullptr;
    // cppcheck-suppress CastIntegerToAddressAtReturn ; false positive
    return strings[idx];
}

} // namespace OpenVolumeMesh::IO
