#pragma once

#include <OpenVolumeMesh/IO/detail/BinaryIStream.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>
#include <OpenVolumeMesh/IO/enums.hh>
#include <OpenVolumeMesh/IO/ReadOptions.hh>
#include <OpenVolumeMesh/IO/PropertyCodec.hh>
#include <OpenVolumeMesh/IO/detail/GeometryReader.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Config/Export.hh>

#include <vector>
#include <limits>
#include <limits>
#include <any>
#include <functional>


namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT BinaryFileReader
{
public:
    BinaryFileReader(std::istream &_s,
                     std::unique_ptr<GeometryReaderBase> _geometry_reader,
                     ReadOptions const& _options,
                     PropertyCodecs const &_prop_codecs = g_default_property_codecs)
        : stream_(_s)
        , geometry_reader_(std::move(_geometry_reader))
        , options_(_options)
        , prop_codecs_(_prop_codecs)
    {}

    std::optional<TopoType> topo_type();
    std::optional<uint8_t> vertex_dim();
    std::optional<VertexEncoding> vertex_encoding();

    template<typename MeshT>
    ReadCompatibility compatibility();

    template<typename MeshT>
    ReadResult read_file(MeshT &_mesh);

private:
    ReadResult internal_read_file(TopologyKernel &_mesh);

    bool read_header();
    void read_chunk();
    void readPropDirChunk(Decoder &reader);
    void readVerticesChunk(Decoder &reader);
    void readEdgesChunk(Decoder &reader);
    void readFacesChunk(Decoder &reader);
    void readCellsChunk(Decoder &reader);
    void readPropChunk(Decoder &reader);

    template<typename HandleT, typename ReadFunc, typename AddFunc>
    void readFacesOrCells(Decoder &reader,
                          TopoChunkHeader const &header,
                          uint8_t fixed_valence,
                          IntEncoding valence_enc,
                          uint64_t n,
                          ReadFunc read_handle,
                          AddFunc add_entity);

    bool validate_span(uint64_t total, uint64_t read, ArraySpan const&span);
    std::vector<uint32_t> read_valences(Decoder &reader, IntEncoding enc, size_t count);

private:
    detail::BinaryIStream stream_;
    std::unique_ptr<GeometryReaderBase> geometry_reader_;
    TopologyKernel *mesh_;
    ReadOptions options_;
    PropertyCodecs const& prop_codecs_;

    FileHeader file_header_;

    ReadState state_ = ReadState::Init;
    bool reached_eof_chunk = false;

    uint64_t n_verts_read_ = 0;
    uint64_t n_edges_read_ = 0;
    uint64_t n_faces_read_ = 0;
    uint64_t n_cells_read_ = 0;


    struct Property {
        Property() = default;
        Property(PropertyEntity _entity,
                 std::shared_ptr<PropertyStorageBase> _prop,
                 const PropertyDecoderBase *_decoder)
                 : entity(_entity)
                 , prop(std::move(_prop))
                 , decoder(_decoder)
        {}
        PropertyEntity entity;
        std::shared_ptr<PropertyStorageBase> prop;
        const PropertyDecoderBase *decoder = nullptr;
    };

    std::vector<Property> props_;

};

} // namespace OpenVolumeMesh::IO::detail


#include <OpenVolumeMesh/IO/detail/BinaryFileReader_impl.hh>
