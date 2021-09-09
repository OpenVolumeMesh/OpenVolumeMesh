#include <OpenVolumeMesh/IO/ovmb_write.hh>
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/IO/detail/GeometryWriter.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileWriter.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>


namespace OpenVolumeMesh::IO {

template<typename MeshT>
bool mesh_is_tetrahedral(MeshT const &_mesh)
{
    // for now, just check face and cell valences, and require at least one cell.
    if (_mesh.n_cells() == 0) {
        return false;
    }
    for (const auto fh: _mesh.faces()) {
        if (_mesh.valence(fh) != 3) {
            return false;
        }
    }
    for (const auto ch: _mesh.cells()) {
        if (_mesh.valence(ch) != 4) {
            return false;
        }
    }
    return true;
    // TODO: check order of faces in cells
}
template<typename MeshT>
bool mesh_is_hexahedral(MeshT const &_mesh)
{
    if (_mesh.n_cells() == 0) {
        return false;
    }
    for (const auto fh: _mesh.faces()) {
        if (_mesh.valence(fh) != 4) {
            return false;
        }
    }
    for (const auto ch: _mesh.cells()) {
        if (_mesh.valence(ch) != 6) {
            return false;
        }
    }
    return true;
    // TODO: check order of faces in cells
}


template<typename MeshT>
WriteOptions::TopologyType detect_topology_type(MeshT const& _mesh)
{
        if (std::is_base_of<TetrahedralMeshTopologyKernel, MeshT>::value) {
            return WriteOptions::TopologyType::Tetrahedral;
        }
        if (std::is_base_of<HexahedralMeshTopologyKernel, MeshT>::value) {
            return WriteOptions::TopologyType::Hexahedral;
        }
        // This is a PolyhedralMesh or something user-defined.
        // Let's look at the contents:
        if (mesh_is_tetrahedral(_mesh)) {
            return WriteOptions::TopologyType::Tetrahedral;
        }
        if (mesh_is_hexahedral(_mesh)) {
            return WriteOptions::TopologyType::Hexahedral;
        }
        return WriteOptions::TopologyType::Polyhedral;
}

template<typename MeshT>
WriteResult ovmb_write(std::ostream &_ostream,
                       MeshT const &_mesh,
                       WriteOptions _options,
                       PropertyCodecs const &_prop_codecs)
{
    using PointT = typename MeshT::PointT;

    if (_options.topology_type == WriteOptions::TopologyType::AutoDetect) {
        _options.topology_type = detect_topology_type(_mesh);
    }
    detail::GeometryWriterT geom_writer(_mesh.vertex_positions());
    detail::BinaryFileWriter writer(
                geom_writer,
                _ostream, _mesh,
                _options, _prop_codecs);
    return writer.write_file();
}

#if 0
template OVM_EXPORT WriteResult ovmb_write(std::ostream &, GeometricPolyhedralMeshV3d const&);
template OVM_EXPORT WriteResult ovmb_write(std::ostream &, GeometricTetrahedralMeshV3d const&);
#endif

} // namespace OpenVolumeMesh::IO
