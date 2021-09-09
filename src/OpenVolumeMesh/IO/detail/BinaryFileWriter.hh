#pragma once

#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/IO/PropertyCodec.hh>
#include <OpenVolumeMesh/IO/enums.hh>
#include <OpenVolumeMesh/IO/WriteOptions.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <OpenVolumeMesh/IO/detail/GeometryWriter.hh>
#include <vector>

#include <string>

namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT BinaryFileWriter
{
public:
    BinaryFileWriter(
            std::ostream &_stream,
            TopologyKernel const &_mesh,
            std::unique_ptr<GeometryWriterBase> _geometry_writer,
            WriteOptions const &_options,
            PropertyCodecs const &_prop_codecs)
        : ostream_(_stream)
        , mesh_(_mesh)
        , geometry_writer_(std::move(_geometry_writer))
        , options_(_options)
        , prop_codecs_(_prop_codecs)
    {
        // preallocate to avoid reallocations
        header_buffer_.need(64);
        chunk_buffer_.need(1024*1024*100);
    }
    WriteResult write_file();
private:
    void write_chunk(ChunkType type);

    void write_vertices(uint32_t first, uint32_t count);
    void write_edges   (uint32_t first, uint32_t count);
    void write_faces   (uint32_t first, uint32_t count);
    void write_cells   (uint32_t first, uint32_t count);

    void write_propdir();
    void write_all_props();

protected:
    std::ostream& ostream_;
    TopologyKernel const& mesh_;
    std::unique_ptr<const GeometryWriterBase> geometry_writer_;
    WriteOptions options_;
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
