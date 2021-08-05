#pragma once

#include <OpenVolumeMesh/IO/BinaryIO.hh>
#include <vector>
#include <string>

namespace OpenVolumeMesh::IO {

template<typename MeshT>
class BinaryFileWriter
{
public:
    BinaryFileWriter(std::ostream &_stream, MeshT &_mesh)
        : writer_(_stream)
        , mesh_(_mesh)
    {}
    bool write();
    /// returns bytes of required padding
    uint8_t write_chunk_header(ChunkType type, uint64_t payload_length);
    void write_vertices(uint64_t first, uint32_t count);
    void write_edges   (uint64_t first, uint32_t count);
    void write_faces   (uint64_t first, uint32_t count);
    void write_cells   (uint64_t first, uint32_t count);

private:
    StreamWriter writer_;
    MeshT const &mesh_;
};


} // namespace OpenVolumeMesh::IO

#include <OpenVolumeMesh/IO/BinaryFileWriterT_impl.hh>
