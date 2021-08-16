#include <OpenVolumeMesh/IO/BinaryFileReader.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReaderImpl.hh>
#include <fstream>
#include <OpenVolumeMesh/Mesh/PolyhedralMesh.hh>
#include <OpenVolumeMesh/Mesh/TetrahedralMesh.hh>
#include <OpenVolumeMesh/Mesh/HexahedralMesh.hh>

namespace OpenVolumeMesh::IO {

BinaryFileReader::BinaryFileReader(std::istream &_s)
    : pimpl_(std::make_unique<detail::BinaryFileReaderImpl>(_s))
{}

BinaryFileReader::BinaryFileReader(const char *filename)
    : fstream_(filename, std::ios::binary)
    , pimpl_(std::make_unique<detail::BinaryFileReaderImpl>(fstream_))
{}

BinaryFileReader::~BinaryFileReader() = default;

void BinaryFileReader::enable_topology_check(bool enabled) {
    pimpl_->enable_topology_check(enabled);
}
void BinaryFileReader::enable_bottom_up_incidences(bool enabled) {
    pimpl_->enable_bottom_up_incidences(enabled);
}

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
