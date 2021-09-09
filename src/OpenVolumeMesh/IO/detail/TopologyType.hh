#pragma once

#include <OpenVolumeMesh/IO/WriteOptions.hh>
#include <type_traits>


namespace OpenVolumeMesh {
    class TetrahedralMeshTopologyKernel;
    class HexahedralMeshTopologyKernel;
}

namespace OpenVolumeMesh::IO::detail {

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

} // namespace OpenVolumeMesh::IO::detail
