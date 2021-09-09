#include <OpenVolumeMesh/IO/detail/BinaryFileReader.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReader_impl.hh>
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_codec.hh>
#include <OpenVolumeMesh/Core/Handles.hh>

#include <istream>
#include <cassert>
#include <numeric>
#include <iostream>

namespace OpenVolumeMesh {
class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO::detail {

template<typename HandleT, typename ReadFunc, typename AddFunc>
void BinaryFileReader::readFacesOrCells(
        Decoder &reader,
        TopoChunkHeader const &header,
        uint8_t fixed_valence,
        IntEncoding valence_enc,

        uint64_t n,
        ReadFunc read_handle,
        AddFunc add_entity)
{
    std::vector<HandleT> handles;
    if (header.span.count == 0) {
        state_ = ReadState::ErrorEmptyList;
        return;
    }

    auto add_all = [&](auto get_valence)
    {
        for (uint64_t i = 0; i < header.span.count; ++i) {
            auto val = get_valence(i);
            handles.resize(val);
            for (uint64_t h = 0; h < val; ++h) {
                size_t idx = read_handle(reader) + header.span.base;
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
        auto valences = read_valences(reader, valence_enc, header.span.count);
        auto max_valence = *std::max_element(valences.begin(), valences.end());
        auto total_handles = std::accumulate(valences.begin(), valences.end(), 0ULL);
        if (reader.remaining_bytes() != total_handles * esize) {
            state_ = ReadState::Error;
            return;
        }
        handles.reserve(max_valence);
        add_all([&](size_t idx){return valences[idx];});
    } else {
        if (reader.remaining_bytes() != header.span.count * fixed_valence * esize) {
            state_ = ReadState::Error;
            return;
        }
        handles.reserve(fixed_valence);
        add_all([val=fixed_valence](size_t){return val;});
    }
}

bool BinaryFileReader::read_header()
{
    if (state_ == ReadState::Init) {
        //stream_.seek(0);
        auto decoder = stream_.make_decoder(ovmb_size<FileHeader>);
        auto ok = read(decoder, file_header_);
        if (ok) {
            state_ = ReadState::HeaderRead;
            return true;
        } else {
            state_ = ReadState::ErrorInvalidFile;
            return false;
        }
    }
    return true;
}


bool BinaryFileReader::validate_span(uint64_t total, uint64_t read, ArraySpan const&span)
{
    if (span.base != read) {
        state_ = ReadState::Error;
        return false;
    }
    if (total - read < span.count) {
        state_ = ReadState::Error;
        return false;
    }
    return true;
}

std::vector<uint32_t> BinaryFileReader::read_valences(Decoder &reader, IntEncoding enc, size_t count)
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

void BinaryFileReader::readPropChunk(Decoder &reader)
{
    PropChunkHeader header;
    read(reader, header);
    if (header.idx >= props_.size()) {
        state_ = ReadState::Error; // TODO more specific error
        return;
    }
    auto &prop = props_[header.idx];
    if (!prop.decoder) {
        reader.skip();
        return;
    }
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
    if (header.span.base >= n || n - header.span.base < header.span.count) {
        state_ = ReadState::ErrorHandleRange;
        return;
    }
    prop.decoder->deserialize(prop.prop.get(),
            reader,
            static_cast<size_t>(header.span.base),
            static_cast<size_t>(header.span.base + header.span.count));
}

void BinaryFileReader::readEdgesChunk(Decoder &reader)
{
    TopoChunkHeader header;
    read(reader, header);
    if (!is_valid(header.enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }
    if (!validate_span(file_header_.n_edges, n_edges_read_, header.span))
        return;


    auto read_all = [&](auto read_one)
    {
        if (reader.remaining_bytes() != header.span.count * 2 * elem_size(header.enc)) {
            std::cerr << "edge chunk size " << std::endl;
            state_ = ReadState::ErrorInvalidChunkSize;
            return;
        }
        for (size_t i = 0; i < header.span.count; ++i) {
            uint64_t src = read_one(reader) + header.span.base;
            uint64_t dst = read_one(reader) + header.span.base;
            if (src >= n_verts_read_ || dst >= n_verts_read_) {
                state_ = ReadState::ErrorHandleRange;
            } else {
                auto vh_src = VertexHandle::from_unsigned(src);
                auto vh_dst = VertexHandle::from_unsigned(dst);
                mesh_->add_edge(vh_src, vh_dst, true);
            }
        }
    };

    call_with_decoder(header.enc, read_all);

    if (state_ == ReadState::ReadingChunks) {
        n_edges_read_ += header.span.count;
    }
}

void BinaryFileReader::readFacesChunk(Decoder &reader)
{
    TopoChunkHeader header;
    read(reader, header);

    uint32_t fixed_valence = reader.u8();
    IntEncoding valence_enc; read(reader, valence_enc);
    reader.reserved<2>();

    if (!is_valid(header.enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }
    if (!validate_span(file_header_.n_faces, n_faces_read_, header.span))
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
        auto add_face =  [&](auto halfedges)
        {
            return mesh_->add_face(std::move(halfedges), options_.topology_check);
        };
        readFacesOrCells<HalfEdgeHandle>(
                    reader, header, fixed_valence, valence_enc, n_edges_read_ * 2,
                    read_one, add_face);
    };

    call_with_decoder(header.enc, read_all);

    if (state_ == ReadState::ReadingChunks) {
        n_faces_read_ += header.span.count;
    }
}

void BinaryFileReader::readCellsChunk(Decoder &reader)
{
    TopoChunkHeader header;
    read(reader, header);

    uint32_t fixed_valence = reader.u8();
    IntEncoding valence_enc; read(reader, valence_enc);
    reader.reserved<2>();

    if (!is_valid(header.enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }

    if (!validate_span(file_header_.n_cells, n_cells_read_, header.span))
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
        auto add_cell =  [&](auto halffaces)
        {
            return mesh_->add_cell(std::move(halffaces), options_.topology_check);
        };
        readFacesOrCells<HalfFaceHandle>(
                    reader, header, fixed_valence, valence_enc, n_faces_read_ * 2,
                    read_one, add_cell);
    };

    call_with_decoder(header.enc, read_all);

    if (state_ == ReadState::ReadingChunks) {
        n_cells_read_ += header.span.count;
    }
}

void
BinaryFileReader::
readVerticesChunk(Decoder &reader)
{
    VertexChunkHeader header;
    read(reader, header);

    if (!is_valid(header.enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }
    if (!validate_span(file_header_.n_verts, n_verts_read_, header.span))
        return;

    auto pos_size = elem_size(file_header_.vertex_encoding) * file_header_.vertex_dim;
    if (reader.remaining_bytes() != header.span.count  * pos_size) {
#if 0
        std::cerr << "vert chunk size" << std::endl;
        std::cerr << "remaining: " << reader.remaining_bytes() << std::endl;
        std::cerr << "expected: " << header.span.count *  pos_size << std::endl;
#endif
        state_ = ReadState::ErrorInvalidChunkSize;
        return;
    }

    geometry_reader_->read(reader, file_header_.vertex_encoding, header.span.base, header.span.count);

    if (state_ == ReadState::ReadingChunks) {
        n_verts_read_ += header.span.count;
    }
}

std::optional<TopoType> BinaryFileReader::topo_type()
{
    if (!read_header())
        return {};
    return file_header_.topo_type;
}

std::optional<uint8_t> BinaryFileReader::vertex_dim()
{
    if (!read_header())
        return {};
    return file_header_.vertex_dim;
}

std::optional<VertexEncoding> BinaryFileReader::vertex_encoding()
{
    if (!read_header())
        return {};
    return file_header_.vertex_encoding;
}


ReadResult BinaryFileReader::internal_read_file(TopologyKernel &out)
{
    mesh_ = &out;
    if (file_header_.n_verts > max_handle_idx
            || file_header_.n_edges > max_handle_idx
            || file_header_.n_faces > max_handle_idx
            || file_header_.n_cells > max_handle_idx)
    {
        return ReadResult::IncompatibleMesh;
    }

    out.clear();
    out.enable_bottom_up_incidences(false);

    out.add_n_vertices(static_cast<int>(file_header_.n_verts));
    out.reserve_edges(static_cast<int>(file_header_.n_edges));
    out.reserve_faces(static_cast<int>(file_header_.n_faces));
    out.reserve_cells(static_cast<int>(file_header_.n_cells));

    state_ = ReadState::ReadingChunks;
    while (stream_.remaining_bytes() > 0) {
        read_chunk();
        if (state_ != ReadState::ReadingChunks) {
            return ReadResult::InvalidFile;
        }
    }
    if (stream_.remaining_bytes() != 0) {
        state_ = ReadState::ErrorEndNotReached;
        return ReadResult::InvalidFile;
    }
    if (!reached_eof_chunk) {
        state_ = ReadState::ErrorEndNotReached;
    }
    if (file_header_.n_verts != out.n_vertices()
            || file_header_.n_edges != out.n_edges()
            || file_header_.n_faces != out.n_faces()
            || file_header_.n_cells != out.n_cells())
    {
        state_ = ReadState::ErrorMissingData;
        return ReadResult::InvalidFile;
    }
    state_ = ReadState::Ok;
    if (options_.bottom_up_incidences) {
        out.enable_bottom_up_incidences(true);
    }
    return ReadResult::Ok;
}


void
BinaryFileReader::
read_chunk()
{
    // TODO: enforce chunk alignment!

    ChunkHeader header;
    auto decoder = stream_.make_decoder(ovmb_size<ChunkHeader>);
    read(decoder, header);
    if (header.file_length > stream_.remaining_bytes()) {
        state_ = ReadState::ErrorChunkTooBig;
        return;
    }
    assert(header.compression == 0);
    assert(header.version == 0);
    auto chunk_reader = stream_.make_decoder(header.payload_length);
    if (header.version != 0) {
        if (header.isMandatory()) {
            state_ = ReadState::ErrorUnsupportedChunkVersion;
            return;
        } else {
            chunk_reader.skip();
        }
    } else {
        using ChunkType = ChunkType;
        switch(header.type)
        {
        case ChunkType::EndOfFile:
            if (header.payload_length != 0) {
                state_ = ReadState::Error;
            }
            if (reached_eof_chunk) {
                state_ = ReadState::Error;
            }
            reached_eof_chunk = true;
            break;
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
    stream_.make_decoder(header.padding_bytes).padding(header.padding_bytes);
}


void
BinaryFileReader::
readPropDirChunk(Decoder &reader)
{
    if (props_.size() != 0) {
        // we can only have one property directory!
        state_ = ReadState::Error; // TODO more specific error
        return;
    }
    while (reader.remaining_bytes() > 0)
    {
        PropertyInfo prop_info;
        read(reader, prop_info);
        auto prop_decoder = prop_codecs_.get_decoder(prop_info.data_type_name);
        if (!prop_decoder) {
            std::cerr << "Could not find decoder for type " << prop_info.data_type_name
                      << ", ignoring."
                      << std::endl;
            props_.emplace_back(); // important to have the right indices
            continue;
        }
        auto prop = prop_decoder->request_property(*mesh_, as_entity_type(prop_info.entity_type), prop_info.name, prop_info.serialized_default);
        props_.emplace_back(prop_info.entity_type, std::move(prop), prop_decoder);
    }
}




} // namespace OpenVolumeMesh::IO::detail
