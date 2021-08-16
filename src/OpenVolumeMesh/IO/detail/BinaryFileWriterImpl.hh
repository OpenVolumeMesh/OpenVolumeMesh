#pragma once

#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/IO/enums.hh>

#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <vector>

#include <string>

namespace OpenVolumeMesh::IO::detail {

template<typename MeshT>
class BinaryFileWriterImpl
{
public:
    BinaryFileWriterImpl(std::ostream &_stream,
                         MeshT const &_mesh,
                         PropertyCodecs const &_prop_codecs = g_default_property_codecs)
        : ostream_(_stream)
        , mesh_(_mesh)
        , prop_codecs_(_prop_codecs)
    {
        // preallocate to avoid reallocations
        header_buffer_.need(64);
        chunk_buffer_.need(1024*1024*100);
    }
    WriteResult write_file();
private:
    /// returns bytes of required padding
    void write_chunk(ChunkType type);
    void write_vertices(uint32_t first, uint32_t count);
    void write_edges   (uint32_t first, uint32_t count);
    void write_faces   (uint32_t first, uint32_t count);
    void write_cells   (uint32_t first, uint32_t count);

    void write_propdir();
    void write_all_props();
private:
    std::ostream &ostream_;
    MeshT const &mesh_;
    PropertyCodecs const &prop_codecs_;

    WriteBuffer chunk_buffer_;
    WriteBuffer header_buffer_;

    struct Property {
        Property(PropertyStorageBase *_prop, const PropertyEncoderBase *_encoder)
            : prop(_prop)
            , encoder(_encoder)
            {}
        PropertyStorageBase *prop;
        const PropertyEncoderBase *encoder;
    };
    std::vector<Property> props_;
};


} // namespace OpenVolumeMesh::IO::detail

#include <OpenVolumeMesh/IO/detail/BinaryFileWriterImplT_impl.hh>
