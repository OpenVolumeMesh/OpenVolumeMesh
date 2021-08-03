#include "BinaryIO.hh"

#include <cassert>
#include <cstring>
#include <limits>

namespace OpenVolumeMesh::IO {

static const std::array<uint8_t, 8> ovmb_magic {'O', 'V', 'M', 'B', '\n', '\r', '\n', 0xff};

uint8_t BufferReader::u8()
{
    std::array<uint8_t, 1> buf;
    read(buf);
    return buf[0];

}
uint16_t BufferReader::u16()
{
    std::array<uint8_t, 2> buf;
    read(buf);
    return (buf[0] + (uint16_t(buf[1]) << 8));

}
uint32_t BufferReader::u32()
{
    std::array<uint8_t, 4> buf;
    read(buf);
    return (buf[0]
            + (uint32_t(buf[1]) << 8)
            + (uint32_t(buf[2]) << 16)
            + (uint32_t(buf[3]) << 24));

}

uint64_t BufferReader::u64()
{
    std::array<uint8_t, 8> buf;
    read(buf);
    return (buf[0]
            + (uint64_t(buf[1]) << 8)
            + (uint64_t(buf[2]) << 16)
            + (uint64_t(buf[3]) << 24)
            + (uint64_t(buf[4]) << 32)
            + (uint64_t(buf[5]) << 40)
            + (uint64_t(buf[6]) << 48)
            + (uint64_t(buf[7]) << 56)
            );
}

double BufferReader::dbl()
{
    std::array<uint8_t, 8> buf;
    read(buf);
    double ret;
    std::memcpy(&ret, buf.data(), sizeof(double));
    return ret;
}

float BufferReader::flt()
{
    std::array<uint8_t, 4> buf;
    read(buf);
    float ret;
    std::memcpy(&ret, buf.data(), sizeof(float));
    return ret;
}


void BufferReader::read(uint8_t *s, size_t n)
{
    if (n <= remaining_bytes()) {
        std::memcpy(s, &data_[pos_], n);
        pos_ += n;
    } else {
        std::runtime_error("end of buffer!");
    }
}


BufferReader::~BufferReader() {
    assert (pos_ == size());
}

uint64_t BufferReader::remaining_bytes() const {
    assert (pos_ < size());
    return size() - pos_;
}

void BufferReader::seek(size_t off) {
    assert(off <= size());
    pos_ = off;
}

void BufferReader::padding(uint8_t n)
{
    std::array<uint8_t, 256> buf;
    read(buf.data(), n);
    for (int i = 0; i < n; ++i) {
        if (buf[i] != 0) {
            throw std::runtime_error("padding not 0");
        }
    }
}



void BufferReader::read(EntityType &_entity_type)
{
    auto v = u8();
    if (v > 6) {
        throw std::runtime_error("unsupported entity type");
    }
    _entity_type = static_cast<EntityType>(v);
}

void StreamWriter::write(EntityType val)
{ u8(static_cast<uint8_t>(val)); }


void BufferReader::read(ChunkType &_chunk_type)
{ _chunk_type = static_cast<ChunkType>(u32()); }

void StreamWriter::write(ChunkType val)
{ u32(static_cast<uint32_t>(val)); }


void BufferReader::read(VertexEncoding &_enc)
{ _enc = static_cast<VertexEncoding>(u8()); }

void StreamWriter::write(VertexEncoding val)
{ u8(static_cast<uint8_t>(val)); }


void BufferReader::read(IntEncoding &_enc)
{ _enc = static_cast<IntEncoding>(u8()); }

void StreamWriter::write(IntEncoding val)
{ u8(static_cast<uint8_t>(val)); }


void BufferReader::read(VertexChunkHeader &_out)
{
    _out.base = u64();
    _out.count = u32();
    read(_out.enc);
    reserved<3>();
}

void StreamWriter::write(const VertexChunkHeader & val)
{
    u64(val.base);
    u32(val.count);
    write(val.enc);
    reserved<3>();
}


void BufferReader::read(TopoChunkHeader &_out)
{
    _out.base = u64();
    _out.count = u32();
    read(_out.enc);
    reserved<3>();
    _out.handle_offset = u64();
}


void StreamWriter::write(const TopoChunkHeader & val)
{
    u64(val.base);
    u32(val.count);
    write(val.enc);
    reserved<3>();
    u64(val.handle_offset);
}

void BufferReader::read(ChunkHeader &header) {
    read(header.type);
    header.version = u8();
    header.padding_bytes = u8();
    header.compression = u8();
    header.flags = static_cast<ChunkFlags>(u8());
    header.file_length = u64();
    assert (header.padding_bytes < header.file_length); // TODO error handling
    header.payload_length = header.file_length - header.padding_bytes;
}

void StreamWriter::write(const ChunkHeader &header) {
    write(header.type);
    u8(header.version);
    u8(header.padding_bytes);
    u8(header.compression);
    u8(static_cast<uint8_t>(header.flags));
    u64(header.file_length);
}

bool BufferReader::read(FileHeader &_file_header)
{
    std::array<uint8_t, 8> magic_buf;
    read(magic_buf);

    if (magic_buf != ovmb_magic) {
        return false;
    }
    _file_header.version = u8();
    _file_header.vertex_dim = u8();
    _file_header.topo_type = static_cast<TopoType>(u8());
    reserved<5>();
    _file_header.n_verts = u64();
    _file_header.n_edges = u64();
    _file_header.n_faces = u64();
    _file_header.n_cells = u64();

    return true;
}

void StreamWriter::write(const FileHeader & header) {
    write(ovmb_magic);
    u8(header.version);
    u8(header.vertex_dim);
    u8(static_cast<uint8_t>(header.topo_type));
    reserved<5>();
    u64(header.n_verts);
    u64(header.n_edges);
    u64(header.n_faces);
    u64(header.n_cells);
}

void StreamWriter::u8(uint8_t val)
{
    s_.write(reinterpret_cast<const char*>(&val), 1);
}

void StreamWriter::u16(uint16_t val)
{
    std::array<uint8_t, 2> data;
    data[0] =  val       & 0xff;
    data[1] = (val >> 8) & 0xff;
    write(data);
}

void StreamWriter::u32(uint32_t val)
{
    std::array<uint8_t, 4> data;
    data[0] =  val        & 0xff;
    data[1] = (val >>  8) & 0xff;
    data[2] = (val >> 16) & 0xff;
    data[3] = (val >> 24) & 0xff;
    write(data);
}

void StreamWriter::u64(uint64_t val)
{
    std::array<uint8_t, 8> data;
    data[0] =  val        & 0xff;
    data[1] = (val >>  8) & 0xff;
    data[2] = (val >> 16) & 0xff;
    data[3] = (val >> 24) & 0xff;
    data[4] = (val >> 32) & 0xff;
    data[5] = (val >> 40) & 0xff;
    data[6] = (val >> 48) & 0xff;
    data[7] = (val >> 56) & 0xff;
    write(data);
}

void StreamWriter::dbl(double val)
{
    std::array<uint8_t, sizeof(double)> data;
    std::memcpy(data.data(), &val, data.size());
    write(data);
}

void StreamWriter::flt(float val)
{
    std::array<uint8_t, sizeof(float)> data;
    std::memcpy(data.data(), &val, data.size());
    write(data);
}

const std::array<char, 256> zero_buf = {0};

void StreamWriter::padding(size_t n)
{
    if (n > zero_buf.size()) {
        throw std::runtime_error("too much padding!");
    }
    s_.write(zero_buf.data(), n);
}

void StreamWriter::write(const char *s, size_t n)
{
    s_.write(s, n);
}



IntEncoding suitable_int_encoding(uint32_t max_value)
{
   if (max_value <= std::numeric_limits<uint8_t>::max())  return IntEncoding::U8;
   if (max_value <= std::numeric_limits<uint16_t>::max())  return IntEncoding::U16;
   return IntEncoding::U32;
}

static uint64_t stream_length(std::istream &_s) {
    auto orig = _s.tellg();
    _s.seekg(0, std::ios_base::end);
    auto res =_s.tellg();
    _s.seekg(orig, std::ios_base::beg);
    return res;
}
BinaryIStream::BinaryIStream(std::istream &_s)
    : BinaryIStream(_s, stream_length(_s))
{
}

BinaryIStream::BinaryIStream(std::istream &_s, uint64_t _size)
    : s_(_s)
    , size_(_size)
{
}

} // namespace OpenVolumeMesh::IO
