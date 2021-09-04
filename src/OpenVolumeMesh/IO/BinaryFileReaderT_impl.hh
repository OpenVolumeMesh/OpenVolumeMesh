#include <OpenVolumeMesh/IO/BinaryFileReader.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReaderImpl.hh>
#include <fstream>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

namespace OpenVolumeMesh::IO {


template<typename MeshT>
ReadCompatibility
BinaryFileReader::compatibility() {
    return pimpl_->compatibility<MeshT>();
}

template<typename MeshT>
ReadResult BinaryFileReader::read_file(MeshT &out)
{
    return pimpl_->read_file(out);
}

#if 0
template OVM_EXPORT ReadCompatibility BinaryFileReader::compatibility<GeometricPolyhedralMeshV3d>();
template OVM_EXPORT ReadResult BinaryFileReader::read_file(GeometricPolyhedralMeshV3d&);

template OVM_EXPORT ReadCompatibility BinaryFileReader::compatibility<GeometricTetrahedralMeshV3d>();
template OVM_EXPORT ReadResult BinaryFileReader::read_file(GeometricTetrahedralMeshV3d &);
#endif

} // namespace OpenVolumeMesh::IO
