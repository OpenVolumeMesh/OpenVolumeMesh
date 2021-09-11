#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/detail/GeometryWriter.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileWriter.hh>
#include <OpenVolumeMesh/IO/detail/TopologyType.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>


namespace OpenVolumeMesh::IO {


template<typename MeshT>
std::unique_ptr<detail::BinaryFileWriter>
make_ovmb_writer(std::ostream & _ostream,
                 MeshT const& _mesh,
                 WriteOptions _options,
                 PropertyCodecs const &_prop_codecs)
{
    if (_options.topology_type == WriteOptions::TopologyType::AutoDetect) {
        _options.topology_type = detail::detect_topology_type(_mesh);
    }
    using GeomWriter = detail::GeometryWriterT<typename MeshT::Point>;
    auto geom_writer = std::make_unique<GeomWriter>(_mesh.vertex_positions());
    return std::make_unique<detail::BinaryFileWriter>(
                _ostream,
                _mesh,
                std::move(geom_writer),
                _options,
                _prop_codecs);

}

template<typename MeshT>
WriteResult ovmb_write(std::ostream &_ostream,
                       MeshT const &_mesh,
                       WriteOptions _options,
                       PropertyCodecs const &_prop_codecs)
{
    auto writer = make_ovmb_writer(_ostream, _mesh, _options, _prop_codecs);
    return writer->write_file();
}

} // namespace OpenVolumeMesh::IO
