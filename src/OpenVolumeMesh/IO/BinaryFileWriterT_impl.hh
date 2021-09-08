#include <OpenVolumeMesh/IO/BinaryFileWriter.hh>
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileWriterImpl.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

namespace OpenVolumeMesh::IO {

template<typename MeshT>
WriteResult ovmb_write(std::ostream &_ostream,
                       MeshT const &_mesh,
                       WriteOptions const& _options,
                       PropertyCodecs const &_prop_codecs)
{
    detail::BinaryFileWriterImpl<MeshT> writer(_ostream, _mesh, _options, _prop_codecs);
    return writer.write_file();
}

#if 0
template OVM_EXPORT WriteResult ovmb_write(std::ostream &, GeometricPolyhedralMeshV3d const&);
template OVM_EXPORT WriteResult ovmb_write(std::ostream &, GeometricTetrahedralMeshV3d const&);
#endif

} // namespace OpenVolumeMesh::IO
