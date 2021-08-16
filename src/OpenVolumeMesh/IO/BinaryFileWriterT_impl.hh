#include <OpenVolumeMesh/IO/BinaryFileWriter.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileWriterImpl.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

namespace OpenVolumeMesh::IO {

template<typename MeshT>
OVM_EXPORT
WriteResult ovmb_write(std::ostream &_ostream, MeshT const &_mesh) {
    detail::BinaryFileWriterImpl<MeshT> writer(_ostream, _mesh);
    return writer.write_file();
}

#if 0
template OVM_EXPORT WriteResult ovmb_write(std::ostream &, GeometricPolyhedralMeshV3d const&);
template OVM_EXPORT WriteResult ovmb_write(std::ostream &, GeometricTetrahedralMeshV3d const&);
#endif

} // namespace OpenVolumeMesh::IO
