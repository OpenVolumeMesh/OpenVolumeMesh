#pragma once

#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <istream>
#include <ostream>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>
#include <limits>
#include <cassert>

namespace OpenVolumeMesh::IO::detail {

template<typename T>
struct is_ovmb_enum : std::false_type {};

template<typename T>
inline constexpr bool is_ovmb_enum_v = is_ovmb_enum<T>::value;


template<typename T>
extern size_t ovmb_size;

enum class IntEncoding : uint8_t {
    None = 0,
    U8 = 1,
    U16 = 2,
    U32 = 4,
};
template<> inline size_t ovmb_size<IntEncoding> = 1;
template<> struct OVM_EXPORT is_ovmb_enum<IntEncoding> : std::true_type {};

inline bool is_valid(IntEncoding enc)
{
    return enc == IntEncoding::None
        || enc == IntEncoding::U8
        || enc == IntEncoding::U16
        || enc == IntEncoding::U32;
}
OVM_EXPORT IntEncoding suitable_int_encoding(uint32_t max_value);

enum class PropertyEntity : uint8_t {
    Vertex   = 0,
    Edge     = 1,
    Face     = 2,
    Cell     = 3,
    HalfEdge = 4,
    HalfFace = 5,
    Mesh     = 6
};
template<> inline size_t ovmb_size<PropertyEntity> = 1;
template<> struct OVM_EXPORT is_ovmb_enum<PropertyEntity> : std::true_type {};

inline bool is_valid(PropertyEntity ent) {
    return static_cast<uint8_t>(ent) <= 6;
}

OVM_EXPORT EntityType as_entity_type(PropertyEntity pe);
OVM_EXPORT PropertyEntity as_prop_entity(EntityType pe);

enum class TopoEntity : uint8_t {
    Edge     = 1,
    Face     = 2,
    Cell     = 3,
};
template<> inline size_t ovmb_size<TopoEntity> = 1;
template<> struct OVM_EXPORT is_ovmb_enum<TopoEntity> : std::true_type {};

inline bool is_valid(TopoEntity ent) {
    return static_cast<uint8_t>(ent) >= 1
        && static_cast<uint8_t>(ent) <= 3;
}

enum class TopoType : uint8_t {
    Polyhedral  = 0,
    Tetrahedral = 1,
    Hexahedral  = 2,
};
template<> inline size_t ovmb_size<TopoType> = 1;
template<> struct OVM_EXPORT is_ovmb_enum<TopoType> : std::true_type {};

inline bool is_valid(TopoType v) {
    return static_cast<uint8_t>(v) <= 2;
}

enum class VertexEncoding : uint8_t {
    None = 0, // topology-only mesh
    Float = 1,
    Double = 2,
};
template<> inline size_t ovmb_size<VertexEncoding> = 1;
template<> struct OVM_EXPORT is_ovmb_enum<VertexEncoding> : std::true_type {};

inline bool is_valid(VertexEncoding enc) {
    return enc == VertexEncoding::None
            || enc == VertexEncoding::Float
            || enc == VertexEncoding::Double;
}

struct PropertyInfo {
    PropertyEntity entity_type;
    std::string name;
    std::string data_type_name;
    std::vector<uint8_t> serialized_default;
};


OVM_EXPORT extern const std::array<uint8_t, 8> ovmb_magic;
struct FileHeader {
    uint8_t file_version = 0;
    uint8_t header_version = 0;
    uint8_t vertex_dim = 0;
    TopoType topo_type = TopoType::Polyhedral;
    // 4 bytes reserved
    uint64_t n_verts = 0;
    uint64_t n_edges = 0;
    uint64_t n_faces = 0;
    uint64_t n_cells = 0;
};
template<> inline size_t ovmb_size<FileHeader> =
        sizeof(ovmb_magic)
        + sizeof(FileHeader::file_version)
        + sizeof(FileHeader::header_version)
        + sizeof(FileHeader::vertex_dim)
        + ovmb_size<TopoType>
        + 4 // reserved
        + 4 * sizeof(uint64_t);

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
    Topo              = FOURCC("TOPO"),
    PropertyDirectory = FOURCC("DIRP"),
    Property          = FOURCC("PROP"),
    EndOfFile         = FOURCC("EOF "),
};
#undef FOURCC

template<> inline size_t ovmb_size<ChunkType> = 4;
template<> struct OVM_EXPORT is_ovmb_enum<ChunkType> : std::true_type {};

inline bool is_valid(ChunkType) {
    // Unkown chunk types will be handed in the chunk reader.
    return true;
}

enum class ChunkFlags : uint8_t {
    Mandatory = 0x1,
};
template<> inline size_t ovmb_size<ChunkFlags> = 1;
template<> struct OVM_EXPORT is_ovmb_enum<ChunkFlags> : std::true_type {};

inline bool is_valid(ChunkFlags flags) {
    return static_cast<uint8_t>(flags) <= 1;
}

struct OVM_EXPORT ChunkHeader {
    ChunkType type;
    uint8_t version;
    uint8_t padding_bytes;
    uint8_t compression;
    ChunkFlags flags;
    uint64_t file_length; // file bytes, excluding ChunkHeader

    /// not in file, computed from length - padding (or from compression header):
    uint64_t payload_length;

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
        + sizeof(ChunkHeader::version)
        + sizeof(ChunkHeader::padding_bytes)
        + sizeof(ChunkHeader::compression)
        + ovmb_size<ChunkFlags>
        + sizeof(ChunkHeader::file_length);


template<typename F>
inline void call_with_encoder(VertexEncoding enc, F const&f)
{
    switch(enc)
    {
    case VertexEncoding::None:  assert(false); break;
    case VertexEncoding::Float:  f([](auto &w, float  v){w.flt(v);}); break;
    case VertexEncoding::Double: f([](auto &w, double v){w.dbl(v);}); break;
    }
}
template<typename F>
inline void call_with_decoder(VertexEncoding enc, F const&f)
{
    switch(enc)
    {
    case VertexEncoding::None:  assert(false); break;
    case VertexEncoding::Float:  f([](auto &r){return r.flt();}); break;
    case VertexEncoding::Double: f([](auto &r){return r.dbl();}); break;
    }
}


template<typename F>
inline void call_with_encoder(IntEncoding enc,F const&f)
{
    switch(enc)
    {
    case IntEncoding::U8:  f([](auto &w, uint8_t  v){w.u8 (v);}); break;
    case IntEncoding::U16: f([](auto &w, uint16_t v){w.u16(v);}); break;
    case IntEncoding::U32: f([](auto &w, uint32_t v){w.u32(v);}); break;
    case IntEncoding::None: assert(false);
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
    case IntEncoding::None: assert(false);
    }
}

OVM_EXPORT constexpr inline uint8_t elem_size(IntEncoding enc) {
    switch(enc) {
    case IntEncoding::U8: return 1;
    case IntEncoding::U16: return 2;
    case IntEncoding::U32: return 4;
    default: return 0;
    }
}
OVM_EXPORT constexpr inline uint8_t elem_size(VertexEncoding enc) {
    switch(enc) {
    case VertexEncoding::Float: return sizeof(float);
    case VertexEncoding::Double: return sizeof(double);
    default: return 0;
    }
}

struct ArraySpan {
    uint64_t first;  // index of first element in this chunk
    uint32_t count; // number of elements in chunk
};
template<> inline size_t ovmb_size<ArraySpan> = sizeof(ArraySpan::first) + sizeof(ArraySpan::count);

struct PropChunkHeader {
    ArraySpan span;
    uint32_t idx; // which property from property dir?
};
template<> inline size_t ovmb_size<PropChunkHeader> = ovmb_size<ArraySpan> + sizeof(PropChunkHeader::idx);

struct VertexChunkHeader {
    ArraySpan span;
    VertexEncoding vertex_encoding;
    // 3 bytes reserved
};
template<> inline size_t ovmb_size<VertexChunkHeader> = ovmb_size<ArraySpan> + 1 + 3;

struct TopoChunkHeader {
    ArraySpan span;
    TopoEntity entity;
    uint8_t valence; // 0 for variable valence
    IntEncoding valence_encoding; // 'None' if valence is not variable
    IntEncoding handle_encoding;
    uint64_t handle_offset; // a value to add to every contained handle, can be used for file compression
};
template<> inline size_t ovmb_size<TopoChunkHeader> =  (
        ovmb_size<ArraySpan> +
        ovmb_size<TopoEntity> +
        sizeof(TopoChunkHeader::valence) +
        ovmb_size<IntEncoding> +
        ovmb_size<IntEncoding> +
        + sizeof(TopoChunkHeader::handle_offset));


} // namespace OpenVolumeMesh::IO::detail
namespace OpenVolumeMesh::IO {

// provide types that appear in public interfaces in non-detail namespace:
using detail::TopoType;
using detail::VertexEncoding;

} // namespace OpenVolumeMesh::IO
