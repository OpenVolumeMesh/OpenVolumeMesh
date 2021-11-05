#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <OpenVolumeMesh/Core/EntityUtils.hh>

#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReader_impl.hh>


#include <sstream>

namespace OpenVolumeMesh {
class TetrahedralMeshTopologyKernel;
class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO::detail {

WriteResult BinaryFileWriter::write_file()
{
    try {
        return do_write_file();
    }  catch (write_error &e) {
        error_msg_ = std::string("write_error: ") + e.what();
        return WriteResult::Error;
    }  catch (std::exception &e) {
        error_msg_ = std::string("exception: ") + e.what();
        return WriteResult::Error;
    }
}

WriteResult BinaryFileWriter::do_write_file()
{
    if (!ostream_.good()) {
        return WriteResult::BadStream;
    }
    using TopoType = TopoType;
    if (mesh_.needs_garbage_collection()) {
        throw write_error("run garbage collection first!");
    }

    // this should be handled before, in make_ovmb_writer
    assert (options_.topology_type != WriteOptions::TopologyType::AutoDetect);

    TopoType topo_type = TopoType::Polyhedral;
    if (options_.topology_type == WriteOptions::TopologyType::Tetrahedral) {
        topo_type = TopoType::Tetrahedral;
    } else  if (options_.topology_type == WriteOptions::TopologyType::Hexahedral) {
        topo_type = TopoType::Hexahedral;
    }

    FileHeader header;
    header.header_version = 1;
    header.file_version = 1;
    header.vertex_dim = geometry_writer_->dim();
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
    if (geometry_writer_->vertex_encoding() != VertexEncoding::None) {
        write_vertices({0, static_cast<uint32_t>(header.n_verts)});
    }
    write_edges({0, static_cast<uint32_t>(header.n_edges)});
    write_faces({0, static_cast<uint32_t>(header.n_faces)});
    write_cells({0, static_cast<uint32_t>(header.n_cells)});
    write_all_props();

    chunk_buffer_.reset();
    write_chunk(ChunkType::EndOfFile);

    if (ostream_.good()) {
        return WriteResult::Ok;
    } else {
        return WriteResult::Error;
    }
}

void BinaryFileWriter::write_chunk(ChunkType type)
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

void BinaryFileWriter::write_vertices(ArraySpan const&_span)
{
    if (_span.count == 0) return;

    VertexChunkHeader header;
    header.span = _span;
    header.vertex_encoding = VertexEncoding::Double;

    chunk_buffer_.reset();
    chunk_buffer_.need(
        ovmb_size<VertexChunkHeader>
        + _span.count * geometry_writer_->elem_size());

    Encoder encoder(chunk_buffer_);
    write(encoder, header);

    geometry_writer_->write(chunk_buffer_, _span);
    write_chunk(ChunkType::Vertices);
}


/// Write header and valence information of a topo chunk.
template<typename ValenceValueOrFunc>
static void start_topo_chunk(WriteBuffer &_buffer,
                             ArraySpan const &_span,
                             TopoEntity _topo_entity,
                             IntEncoding _handle_encoding,
                             ValenceValueOrFunc _valence)
{
    if (_span.count == 0) return;
    auto end = _span.first + _span.count;

    TopoChunkHeader header;
    header.span = _span;
    header.entity = _topo_entity;
    header.handle_encoding = _handle_encoding;
    header.handle_offset = 0;


    uint64_t valence_sum = 0;
    if constexpr (std::is_integral_v<ValenceValueOrFunc>) {
        header.valence = _valence;
        header.valence_encoding = IntEncoding::None;
    } else {
        uint32_t min_valence = std::numeric_limits<uint32_t>::max();
        uint32_t max_valence = std::numeric_limits<uint32_t>::min();

        for (uint64_t idx = _span.first; idx < end; ++idx) {
            auto val = _valence(idx);
            if (val < min_valence) { min_valence = val;}
            if (val > max_valence) { max_valence = val;}
            valence_sum += val;
        }
        if (min_valence == max_valence && min_valence <= std::numeric_limits<uint8_t>::max()) {
            header.valence = min_valence;
            header.valence_encoding = IntEncoding::None;
        } else {
            header.valence = 0;
            header.valence_encoding = suitable_int_encoding(max_valence);
        }
    }

    _buffer.reset();
    Encoder encoder(_buffer);
    write(encoder, header);

    if (header.valence != 0) {
        // fixed valence
        _buffer.need(
                    ovmb_size<TopoChunkHeader>
                    + _span.count * header.valence * elem_size(header.handle_encoding));
    } else {
        if constexpr (std::is_integral_v<ValenceValueOrFunc>) {
            assert(false);
        } else {
            // variable valence
            _buffer.need(
                        ovmb_size<TopoChunkHeader>
                        + _span.count * elem_size(header.valence_encoding)
                        +  valence_sum * elem_size(header.handle_encoding));

            auto write_all = [&](auto write_one) {
                for (uint64_t idx = _span.first; idx < end; ++idx) {
                    write_one(encoder, _valence(idx));
                }
            };
            call_with_encoder(header.valence_encoding, write_all);
        }
    }
}

void BinaryFileWriter::write_edges(ArraySpan const&_span)
{
    if (_span.count == 0) return;
    auto end = _span.first + _span.count;
    assert(end <= mesh_.n_edges());

    auto handle_encoding = suitable_int_encoding(mesh_.n_vertices());

    start_topo_chunk(chunk_buffer_,
                     _span,
                     TopoEntity::Edge,
                     handle_encoding,
                     2);

    Encoder encoder(chunk_buffer_);
    auto write_all = [&](auto write_one) {
        for (uint64_t i = _span.first; i < end; ++i) {
            auto heh = mesh_.halfedge_handle(EdgeHandle::from_unsigned(i), 0);
            write_one(encoder, mesh_.from_vertex_handle(heh).uidx());
            write_one(encoder, mesh_.  to_vertex_handle(heh).uidx());
        }
    };

    call_with_encoder(handle_encoding, write_all);

    write_chunk(ChunkType::Topo);
}

void BinaryFileWriter::write_faces(ArraySpan const&_span)
{
    if (_span.count == 0) return;
    auto end = _span.first + _span.count;
    assert(end <= mesh_.n_faces());

    auto handle_encoding = suitable_int_encoding(mesh_.n_halfedges());

    auto get_valence = [&](uint64_t idx){return mesh_.valence(FH::from_unsigned(idx));};

    start_topo_chunk(chunk_buffer_,
                     _span,
                     TopoEntity::Face,
                     handle_encoding,
                     get_valence);


    Encoder encoder(chunk_buffer_);
    auto write_all = [&](auto write_one)
    {
        for (uint64_t i = _span.first; i < end; ++i) {
            auto fh = FaceHandle::from_unsigned(i);
            for (const auto heh: mesh_.face_halfedges(fh)) {
                write_one(encoder, heh.uidx());
            }
        }
    };

    call_with_encoder(handle_encoding, write_all);
    write_chunk(ChunkType::Topo);
}

void BinaryFileWriter::write_cells(ArraySpan const&_span)
{
    if (_span.count == 0) return;
    auto end = _span.first + _span.count;
    assert(end <= mesh_.n_cells());

    auto handle_encoding = suitable_int_encoding(mesh_.n_halffaces());

    auto get_valence = [&](uint64_t idx){return mesh_.valence(CH::from_unsigned(idx));};

    start_topo_chunk(chunk_buffer_,
                     _span,
                     TopoEntity::Cell,
                     handle_encoding,
                     get_valence);


    Encoder encoder(chunk_buffer_);
    auto write_all = [&](auto write_one)
    {
        for (uint64_t i = _span.first; i < end; ++i) {
            auto ch = CellHandle::from_unsigned(i);
            for (const auto hfh: mesh_.cell_halffaces(ch)) {
                write_one(encoder, hfh.uidx());
            }
        }
    };

    call_with_encoder(handle_encoding, write_all);
    write_chunk(ChunkType::Topo);
}

void BinaryFileWriter::write_propdir()
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

void BinaryFileWriter::write_all_props()
{
    uint32_t idx = 0;
    for (auto prop: props_)
    {
        chunk_buffer_.reset();
        Encoder encoder(chunk_buffer_);

        PropChunkHeader chunk_header;
        chunk_header.span.first = 0;
        chunk_header.span.count = prop.prop->size();
        chunk_header.idx = idx++;
        write(encoder, chunk_header);

        prop.encoder->serialize(prop.prop, chunk_buffer_, chunk_header.span.first, chunk_header.span.first + chunk_header.span.count);

        write_chunk(ChunkType::Property);
    }
}

} // namespace OpenVolumeMesh::IO::detail
