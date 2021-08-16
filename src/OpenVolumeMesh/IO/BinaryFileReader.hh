#pragma once

#include <OpenVolumeMesh/IO/enums.hh>
#include <istream>
#include <fstream>
#include <memory>

namespace OpenVolumeMesh::IO {

namespace detail {
class BinaryFileReaderImpl;
}


class OVM_EXPORT BinaryFileReader
{
public:
    BinaryFileReader(std::istream &_s);
    BinaryFileReader(const char *filename);
    ~BinaryFileReader();

    template<typename MeshT>
    ReadCompatibility compatibility();

    template<typename MeshT>
    ReadResult read_file(MeshT &out);

    void enable_topology_check(bool enabled);
    void enable_bottom_up_incidences(bool enabled);
private:
    std::ifstream fstream_;
    std::unique_ptr<detail::BinaryFileReaderImpl> pimpl_;
};

} // namespace OpenVolumeMesh::IO

#ifndef OVM_DO_NOT_INCLUDE_FILE_READER_IMPL
#  include <OpenVolumeMesh/IO/BinaryFileReaderT_impl.hh>
#endif
