#include <OpenVolumeMesh/IO/detail/BinaryFileReader.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReader_impl.hh>
#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_codec.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>
#include <OpenVolumeMesh/Core/Handles.hh>

#include <istream>
#include <cassert>
#include <numeric>
#include <iostream>
#include <string>

namespace OpenVolumeMesh {
class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO::detail {

bool BinaryFileReader::read_header()
{
    try {
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
    } catch (parse_error &e) {
        state_ = ReadState::ErrorInvalidFile;
        error_msg_ = std::string("parse_error: ") + e.what();
        return false;
    }
}


bool BinaryFileReader::validate_span(uint64_t total, uint64_t read, ArraySpan const&span)
{
    if (span.first != read) {
        state_ = ReadState::Error;
        error_msg_ = "Invalid span start, must resume where the last chunk of the same topo type left off";
        return false;
    }
    if (total - read < span.count) {
        state_ = ReadState::Error;
        error_msg_ = "Invalid span, end exceeds total entity count";
        return false;
    }
    return true;
}

template<typename T, typename FuncMakeT>
bool
BinaryFileReader::read_n_ints(Decoder &reader,
                              IntEncoding enc,
                              std::vector<T> &out_vec,
                              size_t count,
                              FuncMakeT make_t)
{
    if (!is_valid(enc)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return false;
    }
    auto needed = count * elem_size(enc);
    if (reader.remaining_bytes() < needed)
    {
        state_ = ReadState::ErrorInvalidFile;
        error_msg_ = "read_n_ints: not enough data in chunk, need " +
                std::to_string(needed) + ", have " +
                std::to_string(reader.remaining_bytes());
        return false;
    }

    out_vec.clear();
    out_vec.reserve(count);

    auto read_all = [&](auto read_one)
    {
        for (size_t i = 0; i < count; ++i) {
            out_vec.push_back(make_t(read_one(reader)));
        }
    };

    call_with_decoder(enc, read_all);

    return true;
}

void BinaryFileReader::read_topo_chunk(Decoder &reader)
{
    TopoChunkHeader header;
    read(reader, header);

    if (header.span.count == 0) {
        state_ = ReadState::ErrorEmptyList;
        error_msg_ = "TOPO chunk contains no data";
        return;
    }

    if (!is_valid(header.handle_encoding)) {
        state_ = ReadState::ErrorInvalidEncoding;
        error_msg_ = "TOPO chunk: invalid handle encoding";
        return;
    }

    if (!is_valid(header.entity)) {
        state_ = ReadState::ErrorInvalidFile;
        error_msg_ = "TOPO chunk: Invalid topology entity " + std::to_string(static_cast<size_t>(header.entity));
        return;
    }

    if (header.valence != 0 && header.valence_encoding != IntEncoding::None) {
        state_ = ReadState::ErrorInvalidFile;
        error_msg_ = "TOPO edge chunk: valence encoding must be None for fixed valences";
        return;
    }


    ValenceVec valences;

    uint64_t total_handles = header.valence * header.span.count;

    if (header.valence == 0)
    {
        if (header.valence_encoding == IntEncoding::None) {
            state_ = ReadState::ErrorInvalidFile;
            error_msg_ = "TOPO edge chunk: valence encoding must not be None for variable valences";
            return;
        }
        auto success = read_n_ints(reader,
                                   header.valence_encoding,
                                   valences,
                                   header.span.count,
                                   [](uint32_t x){return x;});
        if (!success) {
            return;
        }
        assert(header.span.count == valences.size());
        total_handles = std::accumulate(valences.begin(), valences.end(), 0ULL);
    }
    auto handle_size = elem_size(header.handle_encoding);
    auto expected_bytes = total_handles * handle_size;
    if (reader.remaining_bytes() != expected_bytes) {
        state_ = ReadState::Error;
        error_msg_ = "TOPO chunk: number of remaining bytes incorrect, expecting "
+   std::to_string(expected_bytes) + ", have " + std::to_string(reader.remaining_bytes());
        return;
    }

    switch (header.entity)  {
    case OpenVolumeMesh::IO::detail::TopoEntity::Edge:
        read_edges(reader, header);
        return;
    case OpenVolumeMesh::IO::detail::TopoEntity::Face:
        read_faces(reader, header, valences);
        return;
    case OpenVolumeMesh::IO::detail::TopoEntity::Cell:
        read_cells(reader, header, valences);
        return;
    }
}

void BinaryFileReader::read_edges(Decoder &reader, const TopoChunkHeader &header)
{
    if (!validate_span(file_header_.n_edges, n_edges_read_, header.span)) {
        return;
    }
    if (header.valence != 2) {
        state_ = ReadState::ErrorInvalidFile;
        error_msg_ = "TOPO edge chunk: valence must be 2";
        return;
    }

    auto read_all = [&](auto read_one)
    {
        for (size_t i = 0; i < header.span.count; ++i) {
            uint64_t src = read_one(reader) + header.span.first;
            uint64_t dst = read_one(reader) + header.span.first;
            if (src >= n_verts_read_ || dst >= n_verts_read_) {
                state_ = ReadState::ErrorHandleRange;
                return;
            } else {
                auto vh_src = VertexHandle::from_unsigned(src);
                auto vh_dst = VertexHandle::from_unsigned(dst);
                mesh_->add_edge(vh_src, vh_dst, true);
            }
        }
    };

    call_with_decoder(header.handle_encoding, read_all);

    if (state_ == ReadState::ReadingChunks) {
        n_edges_read_ += header.span.count;
    }
}

void BinaryFileReader::read_faces(Decoder &reader, const TopoChunkHeader &header, const ValenceVec &_valences)
{
    if (!validate_span(file_header_.n_faces, n_faces_read_, header.span))
        return;

    if (file_header_.topo_type == TopoType::Tetrahedral && header.valence != 3) {
        state_ = ReadState::ErrorInvalidTopoType;
        error_msg_ = "TOPO chunk: Faces of tetrahedral meshes must have a fixed valence of 3";
        return;
    }
    if (file_header_.topo_type == TopoType::Hexahedral && header.valence != 4) {
        state_ = ReadState::ErrorInvalidTopoType;
        error_msg_ = "TOPO chunk: Faces of hexahedral meshes must have a fixed valence of 4";
        return;
    }
    assert(header.valence != 0 || _valences.size() == header.span.count);

    auto read_heh = [&](uint32_t x)
    {
        auto idx = x + header.handle_offset;
        if (idx >= 2 * n_edges_read_) {
            throw parse_error("Invalid Halfedge handle");
        }
        return HEH::from_unsigned(idx);
    };
    for (uint64_t i = 0; i < header.span.count; ++i)
    {
        uint32_t valence = header.valence == 0 ? _valences[i] : header.valence;
        std::vector<HEH> halfedges;
        halfedges.reserve(valence);
        auto success = read_n_ints(reader,
                                   header.handle_encoding,
                                   halfedges,
                                   valence,
                                   read_heh);
        if (!success) break;
        mesh_->add_face(std::move(halfedges), options_.topology_check);
    };

    if (state_ == ReadState::ReadingChunks) {
        n_faces_read_ += header.span.count;
    }
}

void BinaryFileReader::read_cells(Decoder &reader, const TopoChunkHeader &header, const ValenceVec &_valences)
{
    if (!validate_span(file_header_.n_cells, n_cells_read_, header.span))
        return;

    if (file_header_.topo_type == TopoType::Tetrahedral && header.valence != 4) {
        state_ = ReadState::ErrorInvalidTopoType;
        error_msg_ = "TOPO chunk: Cells of tetrahedral meshes must have a fixed valence of 4";
        return;
    }

    if (file_header_.topo_type == TopoType::Hexahedral && header.valence != 6) {
        state_ = ReadState::ErrorInvalidTopoType;
        error_msg_ = "TOPO chunk: Cells of hexahedral meshes must have a fixed valence of 6";
        return;
    }

    assert(header.valence != 0 || _valences.size() == header.span.count);

    auto read_hfh = [&](uint32_t x)
    {
        auto idx = x + header.handle_offset;
        if (idx >= 2 * n_faces_read_) {
            throw parse_error("Invalid Halfface handle");
        }
        return HFH::from_unsigned(idx);
    };

    for (uint64_t i = 0; i < header.span.count; ++i)
    {
        uint32_t valence = header.valence == 0 ? _valences[i] : header.valence;
        std::vector<HFH> halffaces;
        halffaces.reserve(valence);
        auto success = read_n_ints(reader,
                                   header.handle_encoding,
                                   halffaces,
                                   valence,
                                   read_hfh);
        if (!success) break;
        mesh_->add_cell(std::move(halffaces), options_.topology_check);
    };

    if (state_ == ReadState::ReadingChunks) {
        n_cells_read_ += header.span.count;
    }

}

void BinaryFileReader::read_prop_chunk(Decoder &reader)
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
    if (header.span.first >= n || n - header.span.first < header.span.count) {
        state_ = ReadState::ErrorHandleRange;
        return;
    }
    prop.decoder->deserialize(prop.prop.get(),
            reader,
            static_cast<size_t>(header.span.first),
            static_cast<size_t>(header.span.first + header.span.count));
}


void
BinaryFileReader::
read_vertices_chunk(Decoder &reader)
{
    VertexChunkHeader header;
    read(reader, header);

    if (!is_valid(header.vertex_encoding)) {
        state_ = ReadState::ErrorInvalidEncoding;
        return;
    }
    if (!validate_span(file_header_.n_verts, n_verts_read_, header.span))
        return;

    auto pos_size = elem_size(header.vertex_encoding) * file_header_.vertex_dim;
    if (reader.remaining_bytes() != header.span.count  * pos_size) {
#if 0
        std::cerr << "vert chunk size" << std::endl;
        std::cerr << "remaining: " << reader.remaining_bytes() << std::endl;
        std::cerr << "expected: " << header.span.count *  pos_size << std::endl;
#endif
        state_ = ReadState::ErrorInvalidChunkSize;
        return;
    }

    geometry_reader_->read(reader, header.vertex_encoding, header.span.first, header.span.count);

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
            read_propdir_chunk(chunk_reader);
            break;
        case ChunkType::Property:
            read_prop_chunk(chunk_reader);
            break;
        case ChunkType::Vertices:
            read_vertices_chunk(chunk_reader);
            break;
        case ChunkType::Topo:
            read_topo_chunk(chunk_reader);
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
    if (!chunk_reader.finished()) {
        state_ = ReadState::ErrorInvalidFile;
        error_msg_ = "Extra data at end of chunk, remaining bytes:" +
                std::to_string(chunk_reader.remaining_bytes());
        return;

    }
    assert(chunk_reader.finished());
    // TODO: this is ugly
    stream_.make_decoder(header.padding_bytes).padding(header.padding_bytes);
}


void
BinaryFileReader::
read_propdir_chunk(Decoder &reader)
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

        EntityType entity_type;
        try {
            entity_type = as_entity_type(prop_info.entity_type);
        } catch (std::runtime_error &e) {
            throw parse_error(e.what());
        }

        auto prop = prop_decoder->request_property(*mesh_, entity_type, prop_info.name, prop_info.serialized_default);
        props_.emplace_back(prop_info.entity_type, std::move(prop), prop_decoder);
    }
}





} // namespace OpenVolumeMesh::IO::detail
