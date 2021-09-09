#pragma once
#include <OpenVolumeMesh/Config/Export.hh>

namespace OpenVolumeMesh::IO {

enum class ReadCompatibility {
    Ok,
    CannotOpenFile,
    BadStream,
    MeshVertexDimensionIncompatible,
    MeshTopologyIncompatible,
    MeshHandleIncompatible,
    FileVersionUnsupported,
    InvalidFile
};
OVM_EXPORT const char* to_string(ReadCompatibility);

enum class ReadState {
    Ok,
    CannotOpenFile,
    BadStream,
    Init,
    HeaderRead,
    ReadingChunks,
    Finished,
    Error,
    ErrorInvalidFile,
    ErrorEndNotReached,
    ErrorIncompatible,
    ErrorChunkTooBig,
    ErrorMissingData,
    ErrorUnsupportedChunkType,
    ErrorUnsupportedChunkVersion,
    ErrorInvalidTopoType,
    ErrorHandleRange,
    ErrorInvalidEncoding,
    ErrorEmptyList,
    ErrorInvalidChunkSize,
};
OVM_EXPORT const char* to_string(ReadState);


enum class ReadResult {
    Ok,
    OtherError,
    InvalidFile,
    CannotOpenFile,
    BadStream,
    IncompatibleMesh,
    };
OVM_EXPORT const char* to_string(ReadResult);

enum class WriteResult {
    Ok,
    Error,
    CannotOpenFile,
    BadStream,
};
OVM_EXPORT const char* to_string(WriteResult);

} // namespace OpenVolumeMesh::IO
