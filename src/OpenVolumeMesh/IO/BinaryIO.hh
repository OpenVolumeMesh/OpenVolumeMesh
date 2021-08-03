#pragma once

#include <istream>
#include <ostream>
#include <array>
#include <cstdint>
#include <memory>

namespace OpenVolumeMesh::IO {

template<typename T>
static inline size_t ovmb_size;

enum class EntityType : uint8_t {
    Vertex = 0,
    Edge = 1,
    Face = 2,
    Cell = 3,
    HalfEdge = 4,
    HalfFace = 5,
    Mesh = 6
};
template<> inline size_t ovmb_size<EntityType> = 1;

enum class TopoType : uint8_t {
    Polyhedral  = 0,
    Tetrahedral = 1,
    Hexahedral  = 2,
};
template<> inline size_t ovmb_size<TopoType> = 1;

struct FileHeader {
    uint8_t version;
    uint8_t vertex_dim;
    TopoType topo_type;
    // TODO: uint64_t file_size; // to detect truncation
    uint64_t n_verts;
    uint64_t n_edges;
    uint64_t n_faces;
    uint64_t n_cells;
};
template<> inline size_t ovmb_size<FileHeader> = 8 + 8 + 4*8;

static constexpr uint32_t fourcc(uint8_t a, uint8_t b, uint8_t c, uint8_t d) {
    return   static_cast<uint32_t>(a)
           |(static_cast<uint32_t>(b) << 8)
           |(static_cast<uint32_t>(c) << 16)
           |(static_cast<uint32_t>(d) << 24);
}
#define FOURCC(abcd) fourcc(abcd[0], abcd[1], abcd[2], abcd[3])

enum class ChunkType : uint32_t{
    Any               = 0, // not used in files!
    Vertices          = FOURCC("VERT"),
    Edges             = FOURCC("EDGE"),
    Faces             = FOURCC("FACE"),
    Cells             = FOURCC("CELL"),
    Property          = FOURCC("PROP")
};
#undef FOURCC

template<> inline size_t ovmb_size<ChunkType> = 4;

enum class ChunkFlags : uint8_t {
    Mandatory = 0x1,
};
template<> inline size_t ovmb_size<ChunkFlags> = 1;

struct ChunkHeader {
    ChunkType type;
    uint8_t version;
    uint8_t padding_bytes;
    uint8_t compression;
    ChunkFlags flags;
    uint64_t file_length; // excluding this header
    uint64_t payload_length; // computed from length - padding (or from compression header)

    bool isFlagSet(ChunkFlags flag) const {
        return (static_cast<uint8_t>(flags)
                & static_cast<uint8_t>(flag))
                == static_cast<uint8_t>(flag);

    }
    bool isMandatory() const {
        return isFlagSet(ChunkFlags::Mandatory);
    }
};
template<> inline size_t ovmb_size<ChunkHeader> =
        ovmb_size<ChunkType>
        + 4 * sizeof(uint8_t)
        + sizeof(uint64_t);


enum class VertexEncoding : uint8_t {
    Float = 1,
    Double = 2,
};
template<> inline size_t ovmb_size<VertexEncoding> = 1;
inline bool is_valid(VertexEncoding enc) {
    return enc == VertexEncoding::Float || enc == VertexEncoding::Double;

}
template<typename F>
inline void call_with_encoder(VertexEncoding enc, F const&f)
{
    switch(enc)
    {
    case VertexEncoding::Float:  f([](auto &w, float  v){w.flt(v);}); break;
    case VertexEncoding::Double: f([](auto &w, double v){w.dbl(v);}); break;
    }
}
template<typename F>
inline void call_with_decoder(VertexEncoding enc, F const&f)
{
    switch(enc)
    {
    case VertexEncoding::Float:  f([](auto &r){return r.flt();}); break;
    case VertexEncoding::Double: f([](auto &r){return r.dbl();}); break;
    }
}

enum class IntEncoding : uint8_t {
    U8 = 1,
    U16 = 2,
    U32 = 4,
};
template<> inline size_t ovmb_size<IntEncoding> = 1;
inline bool is_valid(IntEncoding enc) {
    return enc == IntEncoding::U8 || enc == IntEncoding::U16 || enc == IntEncoding::U32;
}
IntEncoding suitable_int_encoding(uint32_t max_value);

template<typename F>
inline void call_with_encoder(IntEncoding enc,F const&f)
{
    switch(enc)
    {
    case IntEncoding::U8:  f([](auto &w, uint8_t  v){w.u8 (v);}); break;
    case IntEncoding::U16: f([](auto &w, uint16_t v){w.u16(v);}); break;
    case IntEncoding::U32: f([](auto &w, uint32_t v){w.u32(v);}); break;
    }
}
template<typename F>
inline void call_with_decoder(IntEncoding enc, F const&f)
{
    switch(enc)
    {
    case IntEncoding::U8:  f([](auto &r){return r.u8(); }); break;
    case IntEncoding::U16: f([](auto &r){return r.u16();}); break;
    case IntEncoding::U32: f([](auto &r){return r.u32();}); break;
    }
}

inline uint8_t elem_size(IntEncoding enc) {
    switch(enc) {
    case IntEncoding::U8: return 1;
    case IntEncoding::U16: return 2;
    case IntEncoding::U32: return 4;
    default: return 0;
    }
}
inline uint8_t elem_size(VertexEncoding enc) {
    switch(enc) {
    case VertexEncoding::Float: return sizeof(float);
    case VertexEncoding::Double: return sizeof(double);
    default: return 0;
    }
}

struct VertexChunkHeader {
    uint64_t base;
    uint32_t count;
    VertexEncoding enc;
};
template<> inline size_t ovmb_size<VertexChunkHeader> = 8 + 4 + 4;

struct TopoChunkHeader {
    uint64_t base;
    uint32_t count;
    IntEncoding enc;
    uint64_t handle_offset;
};
template<> inline size_t ovmb_size<TopoChunkHeader> = 8 + 4 + 4 + 8;

class StreamWriter {
public:
    StreamWriter(std::ostream &_s)
        : s_(_s)
    {}
public:
    // basic types
    void u8(uint8_t);
    void u16(uint16_t);
    void u32(uint32_t);
    void u64(uint64_t);
    void dbl(double);
    void flt(float);

    // helpers
    void padding(size_t n);
    template<uint8_t N> void reserved();
    template<size_t N> void write(std::array<uint8_t, N> const &arr);
    void write(const char *s, size_t n);

    // ovmb types
    void write(EntityType);
    void write(ChunkType);
    void write(VertexEncoding);
    void write(IntEncoding);
    // ovmb structures
    void write(FileHeader const&);
    void write(ChunkHeader const&);
    void write(VertexChunkHeader const&);
    void write(TopoChunkHeader const&);

private:
    std::ostream &s_;
};

class BufferReader {
public:
    BufferReader(std::unique_ptr<uint8_t[]> _data, size_t _size)
        : data_(std::move(_data))
        , size_(_size)
        , cur_(data_.get())
        , end_(cur_ + size_)
    {}
    ~BufferReader();
public:
// file position handling:
    inline size_t size() const {return size_;};
    size_t remaining_bytes() const;
    inline size_t pos() const {return cur_ - data_.get();}
    void seek(size_t off);
    inline void skip() {seek(size());};

    void need(size_t n);
// basic types
    uint8_t  u8();
    uint16_t u16();
    uint32_t u32();
    uint64_t u64();
    double   dbl();
    float    flt();

// helpers
    void     padding(uint8_t n);
    template<uint8_t N> void reserved();
    template<size_t N> void read(std::array<uint8_t, N> &arr);
    void read(uint8_t *s, size_t n);

// ovmb types
    void read(EntityType &);
    void read(ChunkType &);
    void read(VertexEncoding &);
    void read(IntEncoding &);

// ovmb structures
    /// return true iff the magic number was correct
    bool read(FileHeader &);
    void read(ChunkHeader &);
    void read(VertexChunkHeader &);
    void read(TopoChunkHeader &);

    std::unique_ptr<uint8_t[]> data_;
    size_t size_;
    uint8_t *cur_;
    uint8_t *end_;
};

template<size_t N>
void BufferReader::read(std::array<uint8_t, N> &arr)
{
    read(arr.data(), arr.size());
}
template<size_t N>
void StreamWriter::write(const std::array<uint8_t, N> &arr)
{
    s_.write(reinterpret_cast<const char*>(arr.data()), arr.size());
}

extern const std::array<char, 256> zero_buf;
template<uint8_t N>
void StreamWriter::reserved()
{
    s_.write(zero_buf.data(), N);
}
template<uint8_t N>
void BufferReader::reserved()
{
    std::array<uint8_t, N> buf;
    read(buf);
    for (const auto ch: buf) {
        if (ch != 0) {
            throw std::runtime_error("reserved entry != 0");
        }
    }
}

class BinaryIStream {
public:
    explicit BinaryIStream(std::istream &_s);
    BinaryIStream(std::istream &_s,
                 uint64_t _size);
    size_t remaining_bytes() const {
        return size_ - pos_;
    }
    /// sub-readers share their istream; do not interleave use!
    BufferReader make_reader(size_t n)
    {
        if (remaining_bytes() < n) {
            throw std::runtime_error("make_reader: not enough bytes left.");
        }
        auto buf = std::make_unique<uint8_t[]>(n);
        s_.read(reinterpret_cast<char*>(buf.get()), n);
        pos_ += n;
        return BufferReader(std::move(buf), n);
    }
private:
    std::istream &s_;
    uint64_t size_;
    uint64_t pos_ = 0;
};



} // namespace OpenVolumeMesh::IO
