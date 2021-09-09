#pragma once

#include <OpenVolumeMesh/IO/ovmb_read.hh>
#include <OpenVolumeMesh/FileManager/FileManager.hh>

namespace OpenVolumeMesh::IO {

/// Read a mesh in ovm or ovmb file format depending on the file ending.
template<typename MeshT>
bool read_file(std::string const&_filename, MeshT &_mesh,
               bool _topo_check = true, bool _bottom_up_incidences = true)
{
    if (_filename.size() < 5) {
        return false;
    }
    std::string_view last4{_filename.end() - 4, _filename.end()};
    if (last4 == ".ovm") {
        FileManager file_manager;
        return file_manager.readFile(_filename, _mesh, _topo_check, _bottom_up_incidences);
    } else if (last4 == "ovmb") {
        ReadOptions options;
        options.topology_check = _topo_check;
        options.bottom_up_incidences = _bottom_up_incidences;
        auto result = ovmb_read(_filename.c_str(), _mesh, options);
        return result == ReadResult::Ok;
    } else {
        // unknown file extension
        return false;
    }
}

} // namespace OpenVolumeMesh::IO
