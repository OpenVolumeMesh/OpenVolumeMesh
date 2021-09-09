#pragma once

#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/IO/enums.hh>
#include <OpenVolumeMesh/IO/WriteOptions.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <OpenVolumeMesh/IO/detail/GeometryWriter.hh>
#include <vector>

#include <string>

namespace OpenVolumeMesh::IO::detail {

class BinaryFileWriter
{
public:
    BinaryFileWriter(
            GeometryWriterBase const& _geometry_writer,
            std::ostream &_stream,
            TopologyKernel const &_mesh,
            WriteOptions const &_options = WriteOptions(),
            PropertyCodecs const &_prop_codecs = g_default_property_codecs)
        : geometry_writer_(_geometry_writer)
        , ostream_(_stream)
        , mesh_(_mesh)
        , options_(_options)
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

protected:
    GeometryWriterBase const &geometry_writer_;
    std::ostream& ostream_;
    TopologyKernel const& mesh_;
    WriteOptions const& options_;
    PropertyCodecs const& prop_codecs_;

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

#include <OpenVolumeMesh/IO/detail/BinaryFileWriter.cc>
