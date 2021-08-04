#pragma once

#include <OpenVolumeMesh/IO/BinaryIO.hh>

#include <vector>
#include <limits>
#include <limits>
#include <any>
#include <functional>


namespace OpenVolumeMesh::IO {

enum class ReadCompatibility {
    Ok,
    MeshVertexDimensionIncompatible,
    MeshTopologyIncompatible,
    MeshHandleIncompatible,
    FileVersionUnsupported,
    InvalidFile
};
const char* to_string(ReadCompatibility);

enum class ReadState {
    Init,
    HeaderRead,
    ReadingChunks,
    Finished,
    Error,
    ErrorInvalidMagic,
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
const char* to_string(ReadState);

template<typename MeshT>
class BinaryFileReader
{
public:
    BinaryFileReader(std::istream &_s, MeshT &_mesh)
        : stream_(_s)
        , mesh_(_mesh)
    {}

    ReadCompatibility compatibility();
    bool read_file();
    ReadState state() const {return state_;}
    void enable_topology_check(bool enabled) { topology_check_ = enabled;}
    void enable_bottom_up_incidences(bool enabled) { bottom_up_incidences_ = enabled;}
private:
    void read_header();
    void read_chunk();
    void readVerticesChunk(BufferReader &reader);
    void readEdgesChunk(BufferReader &reader);
    void readFacesChunk(BufferReader &reader);
    void readCellsChunk(BufferReader &reader);
    void readPropDirChunk(BufferReader &reader);
    void readPropChunk(BufferReader &reader);
    template<typename HandleT, typename ReadFunc, typename AddFunc>
    void readFacesOrCells(BufferReader &reader,
                          TopoChunkHeader const &header,
                          uint8_t fixed_valence,
                          IntEncoding valence_enc,
                          uint64_t n,
                          ReadFunc read_handle,
                          AddFunc add_entity);
    bool validate_span(uint64_t total, uint64_t read, uint64_t base, uint64_t count);

    std::vector<uint32_t> read_valences(BufferReader &reader, IntEncoding enc, size_t count);

private:
    BinaryIStream stream_;
    MeshT &mesh_;
    bool topology_check_ = true;
    bool bottom_up_incidences_ = true;

    FileHeader file_header_;
    ReadState state_ = ReadState::Init;
    uint64_t n_verts_read_ = 0;
    uint64_t n_edges_read_ = 0;
    uint64_t n_faces_read_ = 0;
    uint64_t n_cells_read_ = 0;

    struct Property {
        PropertyEntity entity;
        std::function<bool(BufferReader&, size_t base, size_t count)> read_func = nullptr;
    };

    std::vector<Property> props_;

};

} // namespace OpenVolumeMesh::IO


#include "BinaryFileReaderT_impl.hh"
