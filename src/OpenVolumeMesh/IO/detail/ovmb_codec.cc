#include <OpenVolumeMesh/IO/detail/ovmb_codec.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>

#include <type_traits>
#include <cassert>

namespace OpenVolumeMesh::IO::detail {

template<typename Enum>
void read_enum(Decoder& decoder, Enum& out)
{
    std::underlying_type_t<Enum> v;
    decoder.read(v);
    out = static_cast<Enum>(v);
    if (!is_valid(out)) {
        throw parse_error("Invalid enum value");
    }
}

template<typename Enum>
void write_enum(Encoder& encoder, Enum val)
{
    assert(is_valid(val));
    auto v = static_cast<std::underlying_type_t<Enum>>(val);
    encoder.write(v);
}


template void write_enum(Encoder&, IntEncoding);
template void write_enum(Encoder&, PropertyEntity);
template void write_enum(Encoder&, TopoType);
template void write_enum(Encoder&, VertexEncoding);
template void write_enum(Encoder&, ChunkFlags);
template void write_enum(Encoder&, ChunkType);

template void read_enum(Decoder&, IntEncoding& out);
template void read_enum(Decoder&, PropertyEntity& out);
template void read_enum(Decoder&, TopoType& out);
template void read_enum(Decoder&, VertexEncoding& out);
template void read_enum(Decoder&, ChunkFlags& out);
template void read_enum(Decoder&, ChunkType& out);


void write(Encoder &encoder, const FileHeader & header) {
    encoder.write(ovmb_magic);

    encoder.u8(header.file_version);
    encoder.u8(header.header_version);
    encoder.u8(header.vertex_dim);
    write_enum(encoder, header.topo_type);

    encoder.reserved<4>();

    encoder.u64(header.n_verts);
    encoder.u64(header.n_edges);
    encoder.u64(header.n_faces);
    encoder.u64(header.n_cells);
}

bool read(Decoder &decoder, FileHeader &_file_header)
{
    decoder.need(ovmb_size<FileHeader>);
    std::array<uint8_t, 8> magic_buf;
    decoder.read(magic_buf);

    if (magic_buf != ovmb_magic) {
        return false;
    }
    _file_header.file_version = decoder.u8();
    _file_header.header_version = decoder.u8();
    if (_file_header.header_version != 1)
        return false;
    _file_header.vertex_dim = decoder.u8();
    read(decoder, _file_header.topo_type);
    decoder.reserved<4>();
    _file_header.n_verts = decoder.u64();
    _file_header.n_edges = decoder.u64();
    _file_header.n_faces = decoder.u64();
    _file_header.n_cells = decoder.u64();

    return true;
}

void write(Encoder &encoder, const ArraySpan &span) {
    encoder.u64(span.first);
    encoder.u32(span.count);
}

void read(Decoder &decoder, ArraySpan &_out) {
    decoder.need(ovmb_size<ArraySpan>);
    _out.first = decoder.u64();
    _out.count = decoder.u32();
}

void write(Encoder &encoder, const ChunkHeader &header) {
    write(encoder, header.type);
    encoder.u8(header.version);
    encoder.u8(header.padding_bytes);
    encoder.u8(header.compression);
    encoder.u8(static_cast<uint8_t>(header.flags));
    encoder.u64(header.file_length);
}

void read(Decoder &decoder, ChunkHeader &header) {
    decoder.need(ovmb_size<ChunkHeader>);
    read(decoder, header.type);
    header.version = decoder.u8();
    header.padding_bytes = decoder.u8();
    header.compression = decoder.u8();
    read(decoder, header.flags);
    header.file_length = decoder.u64();
    if (header.padding_bytes > header.file_length) {
        throw parse_error("Cannot have more padding than total length");
    }
    header.payload_length = header.file_length - header.padding_bytes;
}

void write(Encoder &encoder, const PropChunkHeader & val)
{
    write(encoder, val.span);
    encoder.u32(val.idx);
}

void read(Decoder &decoder, PropChunkHeader &_out)
{
    decoder.need(ovmb_size<PropChunkHeader>);
    read(decoder, _out.span);
    _out.idx = decoder.u32();
}

void write(Encoder &encoder, const VertexChunkHeader & val)
{
    write(encoder, val.span);
    write(encoder, val.vertex_encoding);
    encoder.reserved<3>();
}

void read(Decoder &decoder, VertexChunkHeader &_out)
{
    decoder.need(ovmb_size<VertexChunkHeader>);
    read(decoder, _out.span);
    read(decoder, _out.vertex_encoding);
    decoder.reserved<3>();
}

void write(Encoder &encoder, const TopoChunkHeader & val)
{
    write(encoder, val.span);
    write(encoder, val.entity);
    encoder.u8(val.valence);
    write(encoder, val.valence_encoding);
    write(encoder, val.handle_encoding);
    encoder.u64(val.handle_offset);
}

void read(Decoder &decoder, TopoChunkHeader &_out)
{
    decoder.need(ovmb_size<TopoChunkHeader>);
    read(decoder, _out.span);
    read(decoder, _out.entity);
    _out.valence = decoder.u8();
    read(decoder, _out.valence_encoding);
    read(decoder, _out.handle_encoding);
    _out.handle_offset = decoder.u64();
}

void write(Encoder &encoder, PropertyInfo const &val) {
    write(encoder, val.entity_type);
    encoder.writeVec<uint32_t>(val.name);
    encoder.writeVec<uint32_t>(val.data_type_name);
    encoder.writeVec<uint32_t>(val.serialized_default);
}

void read(Decoder &decoder, PropertyInfo &_out)
{
    decoder.need(2+3*4);
    read(decoder, _out.entity_type);
    decoder.readVec<uint32_t>(_out.name);
    decoder.readVec<uint32_t>(_out.data_type_name);
    decoder.readVec<uint32_t>(_out.serialized_default);
}

} // namespace OpenVolumeMesh::IO::detail
