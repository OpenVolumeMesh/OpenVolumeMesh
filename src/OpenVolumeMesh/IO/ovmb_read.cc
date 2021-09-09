#include <OpenVolumeMesh/IO/ovmb_read.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReader.hh>

namespace OpenVolumeMesh::IO {

BinaryFileReader::BinaryFileReader(std::istream &_s)
    : pimpl_(std::make_unique<detail::BinaryFileReader>(_s))
{}

BinaryFileReader::BinaryFileReader(const char *filename)
    : fstream_(filename, std::ios::binary)
    , pimpl_(std::make_unique<detail::BinaryFileReader>(fstream_))
{}

BinaryFileReader::~BinaryFileReader() = default;

void BinaryFileReader::enable_topology_check(bool enabled) {
    pimpl_->enable_topology_check(enabled);
}
void BinaryFileReader::enable_bottom_up_incidences(bool enabled) {
    pimpl_->enable_bottom_up_incidences(enabled);
}

} // namespace OpenVolumeMesh::IO
