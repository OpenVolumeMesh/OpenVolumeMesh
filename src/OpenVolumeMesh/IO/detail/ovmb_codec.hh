#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>

#include <istream>
#include <ostream>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>
#include <limits>


namespace OpenVolumeMesh::IO::detail {

class Encoder;
class Decoder;

template<typename Enum>
void write_enum(Encoder& encoder, Enum);

template<typename Enum>
void read_enum(Decoder& decoder, Enum&);

extern template OVM_EXPORT void write_enum(Encoder&, IntEncoding);
extern template OVM_EXPORT void write_enum(Encoder&, PropertyEntity);
extern template OVM_EXPORT void write_enum(Encoder&, TopoType);
extern template OVM_EXPORT void write_enum(Encoder&, VertexEncoding);
extern template OVM_EXPORT void write_enum(Encoder&, ChunkFlags);
extern template OVM_EXPORT void write_enum(Encoder&, ChunkType);

extern template OVM_EXPORT void read_enum(Decoder&, IntEncoding& out);
extern template OVM_EXPORT void read_enum(Decoder&, PropertyEntity& out);
extern template OVM_EXPORT void read_enum(Decoder&, TopoType& out);
extern template OVM_EXPORT void read_enum(Decoder&, VertexEncoding& out);
extern template OVM_EXPORT void read_enum(Decoder&, ChunkFlags& out);
extern template OVM_EXPORT void read_enum(Decoder&, ChunkType& out);

template<typename Enum>
std::enable_if_t<is_ovmb_enum_v<Enum>>
write(Encoder& encoder, Enum val) {
    write_enum(encoder, val);
}

template<typename Enum>
std::enable_if_t<is_ovmb_enum_v<Enum>>
read(Decoder& decoder, Enum& out) {
    read_enum(decoder, out);
}


template<typename Enum>
std::enable_if_t<is_ovmb_enum_v<Enum>>
write(Encoder &encoder, Enum);


OVM_EXPORT void write(Encoder&, const FileHeader&);
OVM_EXPORT bool read (Decoder&, FileHeader&);

OVM_EXPORT void write(Encoder&, const ChunkHeader&);
OVM_EXPORT void read (Decoder&, ChunkHeader&);

OVM_EXPORT void write(Encoder &, const PropChunkHeader&);
OVM_EXPORT void read (Decoder &, PropChunkHeader&);

OVM_EXPORT void write(Encoder &, const VertexChunkHeader&);
OVM_EXPORT void read (Decoder &, VertexChunkHeader&);

OVM_EXPORT void write(Encoder &, const TopoChunkHeader&);
OVM_EXPORT void read (Decoder &, TopoChunkHeader&);

OVM_EXPORT void write(Encoder &, const PropertyInfo&);
OVM_EXPORT void read (Decoder &, PropertyInfo&);


} // namespace OpenVolumeMesh::IO::detail
