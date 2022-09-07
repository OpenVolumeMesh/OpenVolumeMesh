#pragma once

#include <OpenVolumeMesh/IO/detail/BinaryIStream.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>
#include <OpenVolumeMesh/IO/enums.hh>
#include <OpenVolumeMesh/IO/ReadOptions.hh>
#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
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
                     ReadOptions const& _options,
                     PropertyCodecs const &_prop_codecs = g_default_property_codecs)
        : stream_(_s)
        , options_(_options)
        , prop_codecs_(_prop_codecs)
    {}

    std::optional<TopoType> topo_type();
    std::optional<uint8_t> vertex_dim();

    template<typename MeshT>
    ReadCompatibility compatibility();

    template<typename MeshT>
    ReadResult read_file(MeshT &_mesh);

    std::string const &get_error_msg() const {return error_msg_;}

private:
    ReadResult internal_read_file(TopologyKernel &_mesh);

    bool read_header();
    void read_chunk();
    void read_propdir_chunk(Decoder &reader);
    void read_vertices_chunk(Decoder &reader);
    void read_prop_chunk(Decoder &reader);

    using ValenceVec = std::vector<uint32_t>;
    template<typename T, typename FuncMakeT>
    bool read_n_ints(Decoder &reader,
                     IntEncoding enc,
                     std::vector<T> &out_vec,
                     size_t count,
                     FuncMakeT make_t);

    void read_topo_chunk(Decoder &reader);
    void read_edges(Decoder &reader, TopoChunkHeader const &header);
    void read_faces(Decoder &reader, TopoChunkHeader const &header, const ValenceVec &_valences);
    void read_cells(Decoder &reader, TopoChunkHeader const &header, const ValenceVec &_valences);

    bool validate_span(uint64_t total, uint64_t read, ArraySpan const&span);

private:
    detail::BinaryIStream stream_;
    std::unique_ptr<GeometryReaderBase> geometry_reader_;
    TopologyKernel *mesh_ = nullptr;
    ReadOptions options_;
    PropertyCodecs const& prop_codecs_;

    FileHeader file_header_;

    ReadState state_ = ReadState::Init;
    std::string error_msg_;
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
