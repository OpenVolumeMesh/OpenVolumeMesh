#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/IO/enums.hh>
#include <OpenVolumeMesh/IO/WriteOptions.hh>
#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileWriter.hh>
#include <fstream>

namespace OpenVolumeMesh::IO {

template<typename MeshT>
std::unique_ptr<detail::BinaryFileWriter>
make_ovmb_writer(std::ostream & _ostream,
            MeshT const& _mesh,
            WriteOptions options,
            PropertyCodecs const &_prop_codecs);

/// The _istream MUST be opened in using std::ios::binary!
template<typename MeshT>
WriteResult ovmb_write(std::ostream & _ostream,
                       MeshT const& _mesh,
                       WriteOptions options = WriteOptions(),
                       PropertyCodecs const &_prop_codecs = g_default_property_codecs);

template<typename MeshT>
WriteResult ovmb_write(const char *_filename,
                              MeshT const& _mesh,
                              WriteOptions _options = WriteOptions(),
                              PropertyCodecs const &_prop_codecs = g_default_property_codecs)

{
    std::ofstream f(_filename, std::ios::binary);
    if (!f.good()) {
        return WriteResult::CannotOpenFile;
    }
    return ovmb_write(f, _mesh, _options, _prop_codecs);
}

} // namespace OpenVolumeMesh::IO

#ifndef OVM_DO_NOT_INCLUDE_FILE_WRITER_IMPL
#  include <OpenVolumeMesh/IO/ovmb_write_impl.hh>
#endif
