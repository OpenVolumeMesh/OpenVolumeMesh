#pragma once

#include <istream>
#include <ostream>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>
#include <limits>


namespace OpenVolumeMesh {
enum class EntityType;
}

namespace OpenVolumeMesh::IO {

template<typename T>
extern size_t ovmb_size;

enum class PropertyEntity : uint8_t {
    Vertex   = 0,
    Edge     = 1,
    Face     = 2,
    Cell     = 3,
    HalfEdge = 4,
    HalfFace = 5,
    Mesh     = 6
};

EntityType as_entity_type(PropertyEntity pe);
template<> inline size_t ovmb_size<PropertyEntity> = 1;

#if 0
enum class PropertyDataType : uint8_t {
    Bool    =  0,
    U8      =  1,
    U16     =  2,
    U32     =  3,
    U64     =  4,
    S8      =  5,
    S16     =  6,
    S32     =  7,
    S64     =  8,
    Float   =  9,
    Float2  = 10,
    Float3  = 11,
    Float4  = 12,
    Double  = 13,
    Double2 = 14,
    Double3 = 15,
    Double4 = 16,
    String  = 17,
    // complex type serialization:
    CistaPP = 0xff,
};
#endif

struct PropertyInfo {
    PropertyEntity entity_type;
    //PropertyDataType data_type;
    std::string name;
    std::string data_type_name; // for non-standard data types, empty otherwise
    std::vector<uint8_t> serialized_default;
};

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
    PropertyDirectory = FOURCC("DIRP"),
    Property          = FOURCC("PROP"),
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

struct PropChunkHeader {
    uint64_t base;
    uint32_t count;
    uint32_t idx; // which property from property dir?
};
template<> inline size_t ovmb_size<PropChunkHeader> = 8 + 4 + 4;

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
    inline void write(uint8_t  v) {u8(v);}
    inline void write(uint16_t v) {u16(v);}
    inline void write(uint32_t v) {u32(v);}
    inline void write(uint64_t v) {u64(v);}
    inline void write(float v)    {flt(v);}
    inline void write(double v)   {dbl(v);}

    // helpers
    void padding(size_t n);
    template<uint8_t N> void reserved();
    template<size_t N> void write(std::array<uint8_t, N> const &arr);
    void write(const uint8_t *s, size_t n);
    void write(const char *s, size_t n);

    template<typename LengthT, typename Vec>
    void writeVec(Vec const &vec);

    // ovmb types
    //void write(PropertyDataType);
    void write(PropertyEntity);
    void write(ChunkType);
    void write(VertexEncoding);
    void write(IntEncoding);
    // ovmb structures
    void write(FileHeader const&);
    void write(ChunkHeader const&);
    void write(PropChunkHeader const&);
    void write(VertexChunkHeader const&);
    void write(TopoChunkHeader const&);
    void write(PropertyInfo const&);

private:
    std::ostream &s_;
};


template<size_t N>
void StreamWriter::write(const std::array<uint8_t, N> &arr)
{
    s_.write(reinterpret_cast<const char*>(arr.data()), arr.size());
}

extern const std::array<uint8_t, 256> zero_buf;
template<uint8_t N>
void StreamWriter::reserved()
{
    write(zero_buf.data(), N);
}

class BufferReader {
public:
    BufferReader(std::unique_ptr<uint8_t[]> _data, size_t _size)
        : data_(std::move(_data))
        , size_(_size)
        , cur_(data_.get())
        , end_(cur_ + size_)
    {}
public:
// file position handling:
    bool finished() const {return cur_ == end_;}
    inline size_t size() const {return size_;};
    size_t remaining_bytes() const;
    inline size_t pos() const {return cur_ - data_.get();}
    void seek(size_t off);
    inline void skip() {seek(size());};

    void need(size_t n);

// basic types, WARNING: length unchecked!
    uint8_t  u8();
    uint16_t u16();
    uint32_t u32();
    uint64_t u64();
    double   dbl();
    float    flt();
    inline void read(uint8_t  &v) {v = u8();}
    inline void read(uint16_t &v) {v = u16();}
    inline void read(uint32_t &v) {v = u32();}
    inline void read(uint64_t &v) {v = u64();}
    inline void read(double &v)   {v = dbl();}
    inline void read(float &v)    {v = flt();}

// basic types, WARNING: length unchecked!
    void     padding(uint8_t n);
    template<uint8_t N> void reserved();
    template<size_t N> void read(std::array<uint8_t, N> &arr);
    void read(uint8_t *s, size_t n);
    void read(char *s, size_t n);


// basic types, WARNING: length unchecked!
    template<typename LengthT, typename Vec>
    void readVec(Vec &vec);
    //void read(PropertyDataType &);
    void read(PropertyEntity &);
    void read(ChunkType &);
    void read(VertexEncoding &);
    void read(IntEncoding &);

// ovmb structures, length checked!
    /// return true iff the magic number was correct
    bool read(FileHeader &);
    void read(ChunkHeader &);
    void read(PropChunkHeader &);
    void read(VertexChunkHeader &);
    void read(TopoChunkHeader &);
    void read(PropertyInfo &);

    std::unique_ptr<uint8_t[]> data_;
    size_t size_;
    uint8_t *cur_;
    uint8_t *end_;
};
template<typename LengthT, typename Vec>
void BufferReader::readVec(Vec &vec)
{
    need(sizeof(LengthT));
    LengthT len = 0;
    read(len);
    need(len);
    vec.resize(len);
    read(vec.data(), len);
}

template<typename LengthT, typename Vec>
void StreamWriter::writeVec(Vec const &vec)
{
    if (vec.size() > std::numeric_limits<LengthT>::max()) {
        throw std::runtime_error("vector too long for length data type");
        // TODO: maybe just truncate?
    }
    auto len = static_cast<LengthT>(vec.size());
    write(len);
    write(vec.data(), len);
}

template<size_t N>
void BufferReader::read(std::array<uint8_t, N> &arr)
{
    read(arr.data(), arr.size());
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
