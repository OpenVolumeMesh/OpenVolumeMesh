#pragma once
#include "BinaryFileWriter.hh"

#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>
#include <sstream>

namespace OpenVolumeMesh {
class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO {

static void write_valences (StreamWriter &writer,
                            std::vector<uint32_t> const &valences)
{
    auto minmax = std::minmax_element(valences.begin(), valences.end());
    uint32_t minval = *minmax.first;
    uint32_t maxval = *minmax.second;
    if (minval == maxval && minval < 0xff) {
        writer.u8(minval);
        writer.reserved<3>();
        return;
    }
    writer.u8(0); // variable valences
    writer.reserved<3>();
    IntEncoding enc = suitable_int_encoding(maxval);
    writer.write(enc);
    writer.reserved<3>();

    uint64_t bytes_required = 4 + valences.size() * elem_size(enc);

    auto write_all = [&](auto write_one) {
        for (const auto val: valences) {
            write_one(writer, val);
        }
    };

    call_with_encoder(enc, write_all);
}

template<typename MeshT>
bool BinaryFileWriter<MeshT>::write()
{
    if (mesh_.needs_garbage_collection()) {
        throw std::runtime_error("run garbage collection first!");
    }
    TopoType topo_type = TopoType::Polyhedral;

    if (std::is_base_of<TetrahedralMeshTopologyKernel, MeshT>::value) {
        topo_type = TopoType::Tetrahedral;
    }
    if (std::is_base_of<HexahedralMeshTopologyKernel, MeshT>::value) {
        topo_type = TopoType::Hexahedral;
    }

    FileHeader header {
        .version = 0,
        .vertex_dim = 3,
        .topo_type = topo_type,
        .n_verts = mesh_.n_vertices(),
        .n_edges = mesh_.n_edges(),
        .n_faces = mesh_.n_faces(),
        .n_cells = mesh_.n_cells()};
    writer_.write(header);
    write_vertices(0, header.n_verts);
    write_edges(0, header.n_edges);
    write_faces(0, header.n_faces);
    write_cells(0, header.n_cells);
    return true;
}

template<typename MeshT>
uint8_t BinaryFileWriter<MeshT>::write_chunk_header(ChunkType type, uint64_t payload_length)
{
    if (payload_length >= std::numeric_limits<uint64_t>::max() - 7) {
        throw std::runtime_error("too much data");
    }
    uint64_t padded = (payload_length + 7) & ~7LL;
    ChunkHeader chunk_header = {
        .type = type,
        .version = 0,
        .padding_bytes = static_cast<uint8_t>(padded - payload_length),
        .compression = 0,
        .flags = ChunkFlags::Mandatory,
        .file_length = padded,
        .payload_length = payload_length,
    };
    writer_.write(chunk_header);
    return chunk_header.padding_bytes;
}


template<typename MeshT>
void BinaryFileWriter<MeshT>::write_vertices(uint64_t first, uint32_t count)
{
    if (count == 0) return;
    const auto N = MeshT::PointT::dim();

    VertexChunkHeader header = {
        .base = first,
        .count = count,
        .enc = VertexEncoding::Double};

    size_t payload_length = ovmb_size<VertexChunkHeader>
            + count * N * elem_size(header.enc);

    uint8_t pad = write_chunk_header(ChunkType::Vertices, payload_length);
    writer_.write(header);

    auto end = first + count;
    assert(end <= mesh_.n_vertices());

    auto write_all = [&](auto write_one)
    {
        for (uint64_t i = first; i < end; ++i) {
            const auto &pos = mesh_.vertex(VertexHandle::from_unsigned(i));
            for (size_t dim = 0; dim < N; ++dim) {
                write_one(writer_, pos[dim]);
            }
        }
    };
    call_with_encoder(header.enc, write_all);
    writer_.padding(pad);
}

template<typename MeshT>
void BinaryFileWriter<MeshT>::write_edges(uint64_t first, uint32_t count)
{
    if (count == 0) return;

    TopoChunkHeader header = {
        .base = first,
        .count = count,
        .enc = suitable_int_encoding(mesh_.n_vertices()),
        .handle_offset = 0};

    size_t payload_length = ovmb_size<TopoChunkHeader>
            + count * 2 * elem_size(header.enc);

    uint8_t pad = write_chunk_header(ChunkType::Edges, payload_length);
    writer_.write(header);

    auto end = first + count;
    assert(end <= mesh_.n_edges());

    auto write_all = [&](auto write_one)
    {
        for (uint64_t i = first; i < end; ++i) {
            auto heh = mesh_.halfedge_handle(EdgeHandle::from_unsigned(i), 0);
            write_one(writer_, mesh_.from_vertex_handle(heh).uidx());
            write_one(writer_, mesh_.  to_vertex_handle(heh).uidx());
        }
    };
    call_with_encoder(header.enc, write_all);

    writer_.padding(pad);
}

template<typename MeshT>
void BinaryFileWriter<MeshT>::write_faces(uint64_t first, uint32_t count)
{
    if (count == 0) return;

    std::ostringstream ss;
    auto buf_writer = StreamWriter(ss);

    TopoChunkHeader header = {
        .base = first,
        .count = count,
        .enc = suitable_int_encoding(mesh_.n_halfedges()),
        .handle_offset = 0};
    buf_writer.write(header);

    std::vector<uint32_t> valences;
    valences.reserve(count);
    auto end = first + count;
    assert(end <= mesh_.n_faces());

    for (uint64_t i = first; i < end; ++i) {
        auto fh = FaceHandle::from_unsigned(i);
        valences.push_back(mesh_.valence(fh));
    }
    write_valences(buf_writer, valences);

    auto write_all = [&](auto write_one)
    {
        for (uint64_t i = first; i < end; ++i) {
            auto fh = FaceHandle::from_unsigned(i);
            for (const auto heh: mesh_.face_halfedges(fh)) {
                write_one(buf_writer, heh.uidx());
            }
        }
    };
    call_with_encoder(header.enc, write_all);

    auto buf = ss.str();
    uint8_t pad = write_chunk_header(ChunkType::Faces, buf.size());
    writer_.write(buf.data(), buf.size());
    writer_.padding(pad);
}

template<typename MeshT>
void BinaryFileWriter<MeshT>::write_cells(uint64_t first, uint32_t count)
{
    if (count == 0) return;

    std::ostringstream ss;
    auto buf_writer = StreamWriter(ss);

    TopoChunkHeader header = {
        .base = first,
        .count = count,
        .enc = suitable_int_encoding(mesh_.n_halffaces()),
        .handle_offset = 0};
    buf_writer.write(header);

    std::vector<uint32_t> valences;
    valences.reserve(count);
    auto end = first + count;
    assert(end <= mesh_.n_cells());

    for (uint64_t i = first; i < end; ++i) {
        auto ch = CellHandle::from_unsigned(i);
        valences.push_back(mesh_.valence(ch));
    }

    write_valences(buf_writer, valences);

    auto write_all = [&](auto write_one)
    {
        for (uint64_t i = first; i < end; ++i) {
            auto ch = CellHandle::from_unsigned(i);
            for (const auto hfh: mesh_.cell_halffaces(ch)) {
                write_one(buf_writer, hfh.uidx());
            }
        }
    };
    call_with_encoder(header.enc, write_all);

    auto buf = ss.str();
    uint8_t pad = write_chunk_header(ChunkType::Cells, buf.size());
    writer_.write(buf.data(), buf.size());
    writer_.padding(pad);
}

} // namespace OpenVolumeMesh::IO
