#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/IO/enums.hh>
#include <fstream>

namespace OpenVolumeMesh::IO {


template<typename MeshT>
OVM_EXPORT
WriteResult ovmb_write(std::ostream &, MeshT const &);

template<typename MeshT>
inline WriteResult ovmb_write(const char *_filename, MeshT const &_mesh) {
    std::ofstream f(_filename, std::ios::binary);
    if (f.good()) {
        return WriteResult::CannotOpenFile;
    }
    return ovmb_write(f, _mesh);
}

} // namespace OpenVolumeMesh::IO

#ifndef OVM_DO_NOT_INCLUDE_FILE_WRITER_IMPL
#  include <OpenVolumeMesh/IO/BinaryFileWriterT_impl.hh>
#endif
