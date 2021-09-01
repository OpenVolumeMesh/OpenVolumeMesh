#pragma once
#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <OpenVolumeMesh/Core/EntityUtils.hh>

#include <OpenVolumeMesh/IO/BinaryFileWriter.hh>
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_codec.hh>


#include <sstream>

namespace OpenVolumeMesh {
class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO::detail {

static void write_valences (WriteBuffer &_buffer,
                            std::vector<uint32_t> const &_valences)
{
    auto minmax = std::minmax_element(_valences.begin(), _valences.end());
    uint32_t minval = *minmax.first;
    uint32_t maxval = *minmax.second;

    Encoder encoder(_buffer);

    if (minval == maxval && minval < 0xff) {
        encoder.u8(minval);
        encoder.reserved<3>();
        return;
    }
    encoder.u8(0); // variable valences
    auto valence_enc = suitable_int_encoding(maxval);
    write(encoder, valence_enc);
    encoder.reserved<2>();

    size_t bytes_required = 4 + _valences.size() * elem_size(valence_enc);

    _buffer.need(bytes_required);

    auto write_all = [&](auto write_one) {
        for (const auto val: _valences) {
            write_one(encoder, val);
        }
    };

    call_with_encoder(valence_enc, write_all);
}

template<typename MeshT>
WriteResult BinaryFileWriterImpl<MeshT>::write_file()
{
    if (!ostream_.good()) {
        return WriteResult::BadStream;
    }
    using TopoType = TopoType;
    if (mesh_.needs_garbage_collection()) {
        throw std::runtime_error("run garbage collection first!");
    }
    TopoType topo_type = TopoType::Polyhedral;

    if (std::is_base_of<TetrahedralMeshTopologyKernel, MeshT>::value
            || (options_.detect_specialized_mesh && mesh_is_tetrahedral()))
    {
        topo_type = TopoType::Tetrahedral;
    }
    if (std::is_base_of<HexahedralMeshTopologyKernel, MeshT>::value
            || (options_.detect_specialized_mesh && mesh_is_hexahedral()))
    {
        topo_type = TopoType::Hexahedral;
    }

    FileHeader header;
    header.header_version = 0;
    header.file_version = 0;
    header.vertex_dim = 3;
    header.topo_type = topo_type;
    header.n_verts = mesh_.n_vertices();
    header.n_edges = mesh_.n_edges();
    header.n_faces = mesh_.n_faces();
    header.n_cells = mesh_.n_cells();

    header_buffer_.reset();
    Encoder encoder(header_buffer_);
    write(encoder, header);
    header_buffer_.write_to_stream(ostream_);

    write_propdir();
    write_vertices(0, static_cast<uint32_t>(header.n_verts));
    write_edges(0, static_cast<uint32_t>(header.n_edges));
    write_faces(0, static_cast<uint32_t>(header.n_faces));
    write_cells(0, static_cast<uint32_t>(header.n_cells));
    write_all_props();

    chunk_buffer_.reset();
    write_chunk(ChunkType::EndOfFile);

    std::cout << "DEBUG: chunk buffer size: " << chunk_buffer_.allocated_size() << std::endl;
    std::cout << "DEBUG: header buffer size: " << header_buffer_.allocated_size() << std::endl;
    if (ostream_.good()) {
        return WriteResult::Ok;
    } else {
        return WriteResult::Error;
    }
}

template<typename MeshT>
void BinaryFileWriterImpl<MeshT>::write_chunk(ChunkType type)
{
    auto payload_length = chunk_buffer_.size();
    size_t padded = (payload_length + 7) & ~7LL;

    ChunkHeader header;
    header.type = type;
    header.version = 0;
    header.padding_bytes = static_cast<uint8_t>(padded - payload_length);
    header.compression = 0;
    header.flags = ChunkFlags::Mandatory;
    header.file_length = padded;
    header.payload_length = payload_length;

    header_buffer_.reset();
    Encoder encoder(header_buffer_);
    write(encoder, header);
    header_buffer_.write_to_stream(ostream_);
    chunk_buffer_.write_to_stream(ostream_);
    encoder.padding(header.padding_bytes);
    header_buffer_.write_to_stream(ostream_);
}


template<typename MeshT>
void BinaryFileWriterImpl<MeshT>::write_vertices(uint32_t first, uint32_t count)
{
    if (count == 0) return;


    const auto N = MeshT::PointT::dim();

    VertexChunkHeader header;
    header.span.base = first;
    header.span.count = count;
    header.enc = VertexEncoding::Double;

    chunk_buffer_.reset();
    chunk_buffer_.need(
        ovmb_size<VertexChunkHeader>
        + count * N * elem_size(header.enc));

    Encoder encoder(chunk_buffer_);

    write(encoder, header);

    auto end = first + count;
    assert(end <= mesh_.n_vertices());

    auto write_all = [&](auto write_one)
    {
        for (uint32_t i = first; i < end; ++i) {
            const auto &pos = mesh_.vertex(VertexHandle::from_unsigned(i));
            for (size_t dim = 0; dim < N; ++dim) {
                write_one(encoder, pos[dim]);
            }
        }
    };
    call_with_encoder(header.enc, write_all);
    write_chunk(ChunkType::Vertices);
}

template<typename MeshT>
void BinaryFileWriterImpl<MeshT>::write_edges(uint32_t first, uint32_t count)
{
    if (count == 0) return;


    TopoChunkHeader header;
    header.span.base = first;
    header.span.count = count;
    header.enc = suitable_int_encoding(mesh_.n_vertices());
    header.handle_offset = 0;


    chunk_buffer_.reset();
    chunk_buffer_.need(
                ovmb_size<TopoChunkHeader>
                + count * 2 * elem_size(header.enc));

    Encoder encoder(chunk_buffer_);
    write(encoder, header);

    auto end = first + count;
    assert(end <= mesh_.n_edges());

    auto write_all = [&](auto write_one)
    {
        for (uint32_t i = first; i < end; ++i) {
            auto heh = mesh_.halfedge_handle(EdgeHandle::from_unsigned(i), 0);
            write_one(encoder, mesh_.from_vertex_handle(heh).uidx());
            write_one(encoder, mesh_.  to_vertex_handle(heh).uidx());
        }
    };
    call_with_encoder(header.enc, write_all);

    write_chunk(ChunkType::Edges);
}

template<typename MeshT>
void BinaryFileWriterImpl<MeshT>::write_faces(uint32_t first, uint32_t count)
{
    if (count == 0) return;

    chunk_buffer_.reset();
    Encoder encoder(chunk_buffer_);

    TopoChunkHeader header;
    header.span.base = first;
    header.span.count = count;
    header.enc = suitable_int_encoding(mesh_.n_halfedges());
    header.handle_offset = 0;

    write(encoder, header);

    std::vector<uint32_t> valences;
    valences.reserve(count);
    auto end = first + count;
    assert(end <= mesh_.n_faces());

    for (uint32_t i = first; i < end; ++i) {
        auto fh = FaceHandle::from_unsigned(i);
        valences.push_back(mesh_.valence(fh));
    }
    write_valences(chunk_buffer_, valences);

    auto write_all = [&](auto write_one)
    {
        for (uint32_t i = first; i < end; ++i) {
            auto fh = FaceHandle::from_unsigned(i);
            for (const auto heh: mesh_.face_halfedges(fh)) {
                write_one(encoder, heh.uidx());
            }
        }
    };
    call_with_encoder(header.enc, write_all);

    write_chunk(ChunkType::Faces);
}

template<typename MeshT>
void BinaryFileWriterImpl<MeshT>::write_cells(uint32_t first, uint32_t count)
{
    if (count == 0) return;

    chunk_buffer_.reset();
    Encoder encoder(chunk_buffer_);

    TopoChunkHeader header;
    header.span.base = first;
    header.span.count = count;
    header.enc = suitable_int_encoding(mesh_.n_halffaces());
    header.handle_offset = 0;

    write(encoder, header);

    std::vector<uint32_t> valences;
    valences.reserve(count);
    auto end = first + count;
    assert(end <= mesh_.n_cells());

    for (uint32_t i = first; i < end; ++i) {
        auto ch = CellHandle::from_unsigned(i);
        valences.push_back(mesh_.valence(ch));
    }

    write_valences(chunk_buffer_, valences);

    auto write_all = [&](auto write_one)
    {
        for (uint32_t i = first; i < end; ++i) {
            auto ch = CellHandle::from_unsigned(i);
            for (const auto hfh: mesh_.cell_halffaces(ch)) {
                write_one(encoder, hfh.uidx());
            }
        }
    };
    call_with_encoder(header.enc, write_all);

    write_chunk(ChunkType::Cells);
}

template<typename MeshT>
void BinaryFileWriterImpl<MeshT>::write_propdir()
{
    ResourceManager const &resman = mesh_;

    chunk_buffer_.reset();
    Encoder encoder(chunk_buffer_);

    WriteBuffer serialized_default;

    for_each_entity([&](auto entity_tag){
        const auto begin = resman.persistent_props_begin<decltype(entity_tag)>();
        const auto end   = resman.persistent_props_end<decltype(entity_tag)>();

        for (auto prop_it = begin; prop_it != end; ++prop_it)
        {
            PropertyStorageBase *prop_base  = *prop_it;
            auto prop_enc = prop_codecs_.get_encoder(prop_base->internal_type_name());
            if (!prop_enc) {
                std::cerr << "Could not find encoder for property '" << prop_base->name()
                          << "' of type '" << prop_base->internal_type_name()
                          << "', ignoring."
                          << std::endl;
                continue;
            }
            serialized_default.reset();
            prop_enc->serialize_default(prop_base, serialized_default);

            PropertyInfo prop_info;
            prop_info.entity_type = as_prop_entity(prop_base->entity_type());
            prop_info.serialized_default = serialized_default.vec();
            prop_info.name = prop_base->name();
            prop_info.data_type_name = prop_enc->ovmb_type_name();
            write(encoder, prop_info);

            props_.emplace_back(prop_base, prop_enc);
        }
    });

    if (chunk_buffer_.size() == 0)
        return;

    write_chunk(ChunkType::PropertyDirectory);
}

template<typename MeshT>
void BinaryFileWriterImpl<MeshT>::write_all_props()
{
    uint32_t idx = 0;
    for (auto prop: props_)
    {
        chunk_buffer_.reset();
        Encoder encoder(chunk_buffer_);

        PropChunkHeader chunk_header;
        chunk_header.span.base = 0;
        chunk_header.span.count = prop.prop->size();
        chunk_header.idx = idx++;
        write(encoder, chunk_header);

        prop.encoder->serialize(prop.prop, chunk_buffer_, chunk_header.span.base, chunk_header.span.base + chunk_header.span.count);

        write_chunk(ChunkType::Property);
    }
}

template<typename MeshT>
bool
BinaryFileWriterImpl<MeshT>::
mesh_is_tetrahedral() {
    // for now, just check face and cell valences, and require at least one cell.
    if (mesh_.n_cells() == 0) {
        return false;
    }
    for (const auto fh: mesh_.faces()) {
        if (mesh_.valence(fh) != 3) {
            return false;
        }
    }
    for (const auto ch: mesh_.cells()) {
        if (mesh_.valence(ch) != 4) {
            return false;
        }
    }
    return true;
    // TODO: check order of faces in cells
}
template<typename MeshT>
bool
BinaryFileWriterImpl<MeshT>::
mesh_is_hexahedral() {
    if (mesh_.n_cells() == 0) {
        return false;
    }
    for (const auto fh: mesh_.faces()) {
        if (mesh_.valence(fh) != 4) {
            return false;
        }
    }
    for (const auto ch: mesh_.cells()) {
        if (mesh_.valence(ch) != 6) {
            return false;
        }
    }
    return true;
    // TODO: check order of faces in cells
}

} // namespace OpenVolumeMesh::IO::detail
