#pragma once

#include <OpenVolumeMesh/IO/detail/BinaryIStream.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>
#include <OpenVolumeMesh/IO/enums.hh>
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>
#include <OpenVolumeMesh/Config/Export.hh>

#include <vector>
#include <limits>
#include <limits>
#include <any>
#include <functional>


namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT BinaryFileReaderImpl
{
public:
    BinaryFileReaderImpl(std::istream &_s)
        : stream_(_s)
    {}

    template<typename MeshT>
    ReadCompatibility compatibility();

    template<typename MeshT>
    ReadResult read_file(MeshT &_mesh);

//    ReadResult state() const {return state_;}
    void enable_topology_check(bool enabled) { topology_check_ = enabled;}
    void enable_bottom_up_incidences(bool enabled) { bottom_up_incidences_ = enabled;}
private:

    void read_header();

    template<typename MeshT>
    void read_chunk(MeshT &out);

    template<typename MeshT>
    void readPropDirChunk(Decoder &reader, MeshT &out);

    template<typename PointT, typename AddVertexFunc>
    void readVerticesChunk(Decoder &reader, AddVertexFunc);

    template<typename AddEdgeFunc>
    void readEdgesChunk(Decoder &reader, AddEdgeFunc);

    template<typename AddFaceFunc>
    void readFacesChunk(Decoder &reader, AddFaceFunc);

    template<typename AddCellFunc>
    void readCellsChunk(Decoder &reader, AddCellFunc);

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
    const PropertyCodecs *prop_codecs_ = &g_default_property_codecs; // TODO: setter

    bool topology_check_ = true;
    bool bottom_up_incidences_ = true;

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


#include <OpenVolumeMesh/IO/detail/BinaryFileReaderImplT_impl.hh>
