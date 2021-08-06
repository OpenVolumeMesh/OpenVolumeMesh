#pragma once

#include <OpenVolumeMesh/IO/BinaryFileReader.hh>
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <istream>
#include <cassert>
#include <numeric>
#include <OpenVolumeMesh/Core/OpenVolumeMeshHandle.hh>

#include <iostream>

namespace OpenVolumeMesh {
class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO {


template<typename MeshT>
ReadCompatibility BinaryFileReader<MeshT>::compatibility()
{
    read_header();
    if (state_ == ReadState::Error) {
        return ReadCompatibility::InvalidFile;
    }
    if (file_header_.version > 0) {
        return ReadCompatibility::FileVersionUnsupported;
    }
    // TODO: allow converting vertex dimension?
    if (file_header_.vertex_dim != MeshT::PointT::dim()) {
        return ReadCompatibility::MeshVertexDimensionIncompatible;
    }
    if (std::is_base_of<TetrahedralMeshTopologyKernel, MeshT>::value) {
        if (file_header_.topo_type != TopoType::Tetrahedral) {
            return ReadCompatibility::MeshTopologyIncompatible;
        }
    }
    if (std::is_base_of<HexahedralMeshTopologyKernel, MeshT>::value) {
        if (file_header_.topo_type != TopoType::Hexahedral) {
            return ReadCompatibility::MeshTopologyIncompatible;
        }
    }
    size_t max_index = std::max(
                std::max(file_header_.n_verts, file_header_.n_edges),
                std::max(file_header_.n_faces, file_header_.n_cells));

    auto max_handle_idx = static_cast<size_t>(std::numeric_limits<int>::max());
    if (max_index > max_handle_idx)
    {
        return ReadCompatibility::MeshHandleIncompatible;
    }
    return ReadCompatibility::Ok;
}

template<typename MeshT>
bool BinaryFileReader<MeshT>::read_file()
{
    read_header();
    if (compatibility() != ReadCompatibility::Ok) {
        state_ = ReadState::ErrorIncompatible;
        return false;
    }
    if (state_ != ReadState::HeaderRead) {
        return false;
    }
    mesh_.clear();
    mesh_.enable_bottom_up_incidences(false);
    mesh_.reserve_vertices(file_header_.n_verts);
    mesh_.reserve_edges(file_header_.n_edges);
    mesh_.reserve_faces(file_header_.n_faces);
    mesh_.reserve_cells(file_header_.n_cells);

    state_ = ReadState::ReadingChunks;
    while (stream_.remaining_bytes() > 0) {
        read_chunk();
        if (state_ != ReadState::ReadingChunks) {
            return false;
        }
    }
    if (stream_.remaining_bytes() != 0) {
        state_ = ReadState::ErrorEndNotReached;
        return false;
    }
    if (file_header_.n_verts != mesh_.n_vertices()
            || file_header_.n_edges != mesh_.n_edges()
            || file_header_.n_faces != mesh_.n_faces()
            || file_header_.n_cells != mesh_.n_cells())
    {
        state_ = ReadState::ErrorMissingData;
        return false;
    }
    state_ = ReadState::Finished;
    if (bottom_up_incidences_)
        mesh_.enable_bottom_up_incidences(true);
    return true;
}

template<typename MeshT>
void BinaryFileReader<MeshT>::read_header()
{
    if (state_ == ReadState::Init) {
        //stream_.seek(0);
        auto valid_magic = stream_.make_reader(ovmb_size<FileHeader>).read(file_header_);
        if (valid_magic) {
            state_ = ReadState::HeaderRead;
        } else {
            state_ = ReadState::ErrorInvalidMagic;
        }
    }
}


template<typename MeshT>
void BinaryFileReader<MeshT>::read_chunk()
{
    // TODO: enforce chunk alignment!

    ChunkHeader header;
    stream_.make_reader(ovmb_size<ChunkHeader>).read(header);
    if (header.file_length > stream_.remaining_bytes()) {
        state_ = ReadState::ErrorChunkTooBig;
        return;
    }
    assert(header.compression == 0);
    assert(header.version == 0);
    auto chunk_reader = stream_.make_reader(header.payload_length);
    if (header.version != 0) {
        if (header.isMandatory()) {
            state_ = ReadState::ErrorUnsupportedChunkVersion;
            return;
        } else {
            chunk_reader.skip();
        }
    } else {

        switch(header.type)
        {
        case ChunkType::PropertyDirectory:
            readPropDirChunk(chunk_reader);
            break;
        case ChunkType::Property:
            readPropChunk(chunk_reader);
            break;
        case ChunkType::Vertices:
            readVerticesChunk(chunk_reader);
            break;
        case ChunkType::Edges:
            readEdgesChunk(chunk_reader);
            break;
        case ChunkType::Faces:
            readFacesChunk(chunk_reader);
            break;
        case ChunkType::Cells:
            readCellsChunk(chunk_reader);
            break;
        default:
            if (header.isMandatory()) {
                state_ = ReadState::ErrorUnsupportedChunkType;
            } else {
                chunk_reader.skip();
            }
            break;
        }
    }
    if (state_ != ReadState::ReadingChunks)
        return;
    assert(chunk_reader.finished());
    // TODO: this is ugly
    stream_.make_reader(header.padding_bytes).padding(header.padding_bytes);
}

template<typename MeshT>
void BinaryFileReader<MeshT>::readVerticesChunk(BufferReader &reader)
{
    VertexChunkHeader header;
    reader.read(header);

    if (!is_valid(header.enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }
    if (!validate_span(file_header_.n_verts, n_verts_read_, header.base, header.count))
        return;


    auto read_all = [&](auto read_one)
    {
        const auto N = MeshT::PointT::dim();
        if (reader.remaining_bytes() != header.count * N * elem_size(header.enc)) {
            std::cerr << "vert chunk size" << std::endl;
            std::cerr << "remaining: " << reader.remaining_bytes() << std::endl;
            std::cerr << "expected: " << header.count * N * elem_size(header.enc) << std::endl;
            state_ = ReadState::ErrorInvalidChunkSize;
            return;
        }
        typename MeshT::PointT point;
        for (size_t i = 0; i < header.count; ++i) {
            for (size_t dim = 0; dim < N; ++dim) {
                point[dim] = read_one(reader);
            }
            mesh_.add_vertex(point);
        }
    };

    call_with_decoder(header.enc, read_all);

    if (state_ == ReadState::ReadingChunks) {
        n_verts_read_ += header.count;
        assert(mesh_.n_vertices() == n_verts_read_);
    }
}

template<typename MeshT>
void BinaryFileReader<MeshT>::readEdgesChunk(BufferReader &reader)
{
    TopoChunkHeader header;
    reader.read(header);
    if (!is_valid(header.enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }
    if (!validate_span(file_header_.n_edges, n_edges_read_, header.base, header.count))
        return;


    auto read_all = [&](auto read_one)
    {
        if (reader.remaining_bytes() != header.count * 2 * elem_size(header.enc)) {
            std::cerr << "edge chunk size " << std::endl;
            state_ = ReadState::ErrorInvalidChunkSize;
            return;
        }
        const auto n_verts = mesh_.n_vertices();
        for (size_t i = 0; i < header.count; ++i) {
            uint64_t src = read_one(reader) + header.base;
            uint64_t dst = read_one(reader) + header.base;
            if (src >= n_verts || dst >= n_verts) {
                state_ = ReadState::ErrorHandleRange;
            } else {
                auto vh_src = VertexHandle::from_unsigned(src);
                auto vh_dst = VertexHandle::from_unsigned(dst);
                mesh_.add_edge(vh_src, vh_dst, true);
            }
        }
    };

    call_with_decoder(header.enc, read_all);

    if (state_ == ReadState::ReadingChunks) {
        n_edges_read_ += header.count;
        assert(mesh_.n_edges() == n_edges_read_);
    }
}

template<typename MeshT>
void BinaryFileReader<MeshT>::readFacesChunk(BufferReader &reader)
{
    TopoChunkHeader header;
    reader.read(header);

    uint32_t fixed_valence = reader.u8();
    IntEncoding valence_enc; reader.read(valence_enc);
    reader.reserved<2>();

    if (!is_valid(header.enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }
    if (!validate_span(file_header_.n_faces, n_faces_read_, header.base, header.count))
        return;
    if (file_header_.topo_type == TopoType::Tetrahedral && fixed_valence != 3) {
        state_ = ReadState::ErrorInvalidTopoType;
        return;
    }
    if (file_header_.topo_type == TopoType::Hexahedral && fixed_valence != 4) {
        state_ = ReadState::ErrorInvalidTopoType;
        return;
    }
    auto read_all = [&](auto read_one)
    {
        readFacesOrCells<HalfEdgeHandle>(
                    reader, header, fixed_valence, valence_enc, mesh_.n_halfedges(),
                    read_one, [&](auto handles) {mesh_.add_face(handles, topology_check_);});
    };

    call_with_decoder(header.enc, read_all);

    if (state_ == ReadState::ReadingChunks) {
        n_faces_read_ += header.count;
        assert(mesh_.n_faces() == n_faces_read_);
    }
}

template<typename MeshT>
void BinaryFileReader<MeshT>::readCellsChunk(BufferReader &reader)
{
    TopoChunkHeader header;
    reader.read(header);

    uint32_t fixed_valence = reader.u8();
    IntEncoding valence_enc; reader.read(valence_enc);
    reader.reserved<2>();

    if (!is_valid(header.enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }

    if (!validate_span(file_header_.n_cells, n_cells_read_, header.base, header.count))
        return;

    if (file_header_.topo_type == TopoType::Tetrahedral && fixed_valence != 4) {
        state_ = ReadState::ErrorInvalidTopoType;
        return;
    }

    if (file_header_.topo_type == TopoType::Hexahedral && fixed_valence != 6) {
        state_ = ReadState::ErrorInvalidTopoType;
        return;
    }
    auto read_all = [&](auto read_one)
    {
        readFacesOrCells<HalfFaceHandle>(
                    reader, header, fixed_valence, valence_enc, mesh_.n_halffaces(),
                    read_one, [&](auto handles) {mesh_.add_cell(handles, topology_check_);});
    };

    call_with_decoder(header.enc, read_all);

    if (state_ == ReadState::ReadingChunks) {
        n_cells_read_ += header.count;
        assert(mesh_.n_cells() == n_cells_read_);
    }
}

template<typename MeshT>
bool BinaryFileReader<MeshT>::validate_span(uint64_t total, uint64_t read, uint64_t base, uint64_t count)
{
    if (base != read) {
        state_ = ReadState::Error;
        return false;
    }
    if (total - read < count) {
        state_ = ReadState::Error;
        return false;
    }
    return true;
}

template<typename MeshT>
std::vector<uint32_t> BinaryFileReader<MeshT>::read_valences(BufferReader &reader, IntEncoding enc, size_t count)
{
    if (!is_valid(enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return {};
    }

    std::vector<uint32_t> valences;
    valences.reserve(count);
    auto read_all = [&](auto read_one)
    {
        if (reader.remaining_bytes() < count * elem_size(enc)) {
            state_ = ReadState::Error;
            return;
        }
        for (size_t i = 0; i < count; ++i) {
            valences.push_back(read_one(reader));
        }
    };

    call_with_decoder(enc, read_all);

    return valences;
}

template<typename MeshT>
template<typename HandleT, typename ReadFunc, typename AddFunc>
void BinaryFileReader<MeshT>::readFacesOrCells(
        BufferReader &reader,
        TopoChunkHeader const &header,
        uint8_t fixed_valence,
        IntEncoding valence_enc,

        uint64_t n,
        ReadFunc read_handle,
        AddFunc add_entity)
{
    std::vector<HandleT> handles;
    if (header.count == 0) {
        state_ = ReadState::ErrorEmptyList;
        return;
    }

    auto add_all = [&](auto get_valence)
    {
        for (uint64_t i = 0; i < header.count; ++i) {
            auto val = get_valence(i);
            handles.resize(val);
            for (uint64_t h = 0; h < val; ++h) {
                size_t idx = read_handle(reader) + header.base;
                if (idx < n) {
                    handles[h] = HandleT::from_unsigned(idx);
                } else {
                    state_ = ReadState::ErrorHandleRange;
                    return;
                }
            }
            add_entity(handles);
        }

    };

    auto esize = elem_size(header.enc);

    if (fixed_valence == 0) {
        auto valences = read_valences(reader, valence_enc, header.count);
        auto max_valence = *std::max_element(valences.begin(), valences.end());
        auto total_handles = std::accumulate(valences.begin(), valences.end(), 0ULL);
        if (reader.remaining_bytes() != total_handles * esize) {
            state_ = ReadState::Error;
            return;
        }
        handles.reserve(max_valence);
        add_all([&](size_t idx){return valences[idx];});
    } else {
        if (reader.remaining_bytes() != header.count * fixed_valence * esize) {
            state_ = ReadState::Error;
            return;
        }
        handles.reserve(fixed_valence);
        add_all([val=fixed_valence](size_t){return val;});
    }
}

template<typename MeshT>
void BinaryFileReader<MeshT>::readPropDirChunk(BufferReader &reader)
{
    reader.need(4);
    uint32_t count = reader.u32();

    if (props_.size() != 0) {
        // we can only have one property directory!
        state_ = ReadState::Error; // TODO more specific error
        return;
    }

    if (count == 0) {
        state_ = ReadState::ErrorEmptyList;
        return;
    }
    reader.need(count * 8);
    props_.reserve(count);
    while (stream_.remaining_bytes() > 0)
    {
        PropertyInfo prop_info;
        reader.read(prop_info);
        auto dec_it = property_dec.find(prop_info.data_type_name);
        if (dec_it == property_dec.end()) {
            std::cerr << "Could not find decoder for type " << prop_info.data_type_name
                      << ", ignoring."
                      << std::endl;
            props_.push_back({}); // important to have the right indices
            continue;
        }
        PropertyDecoderBase *decoder = dec_it->second.get();
        decoder->request_prop(mesh_, as_entity_type(prop_info.entity_type), prop_info.name, prop_info.serialized_default);
        props_.emplace_back(Property{prop_info.entity_type, decoder});
    }
}

template<typename MeshT>
void BinaryFileReader<MeshT>::readPropChunk(BufferReader &reader)
{
    PropChunkHeader header;
    reader.read(header);
    if (header.idx >= props_.size()) {
        state_ = ReadState::Error; // TODO more specific error
        return;
    }
    auto opt_prop = props_[header.idx];
    if (!opt_prop.has_value()) {
        reader.skip();
        return;
    }
    Property prop = *opt_prop;
    uint64_t n = 0;
    switch (prop.entity) {
        case PropertyEntity::Vertex:   n = n_verts_read_; break;
        case PropertyEntity::Edge:     n = n_edges_read_; break;
        case PropertyEntity::Face:     n = n_faces_read_; break;
        case PropertyEntity::Cell:     n = n_cells_read_; break;
        case PropertyEntity::HalfEdge: n = 2 * n_edges_read_; break;
        case PropertyEntity::HalfFace: n = 2 * n_faces_read_; break;
        case PropertyEntity::Mesh:     n = 1; break;
    }
    if (header.base >= n || n - header.base < header.count) {
        state_ = ReadState::ErrorHandleRange;
        return;
    }
    prop.decoder->deserialize(reader, header.base, header.base + header.count);
}




} // namespace OpenVolumeMesh::IO
