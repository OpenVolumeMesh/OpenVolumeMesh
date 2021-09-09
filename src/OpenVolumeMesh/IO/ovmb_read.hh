#pragma once

#include <OpenVolumeMesh/IO/enums.hh>
#include <OpenVolumeMesh/IO/ReadOptions.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReader.hh>
#include <istream>
#include <fstream>
#include <memory>

namespace OpenVolumeMesh::IO {

template<typename MeshT>
std::unique_ptr<detail::BinaryFileReader>
make_ovmb_reader(std::istream & _ostream,
                 ReadOptions const &_options,
                 PropertyCodecs const &_prop_codecs);

template<typename MeshT>
ReadResult ovmb_read(std::istream & _istream,
                 MeshT & _mesh,
                 ReadOptions _options = ReadOptions(),
                 PropertyCodecs const &_prop_codecs = g_default_property_codecs);

template<typename MeshT>
ReadResult ovmb_read(const char *_filename,
                 MeshT & _mesh,
                 ReadOptions _options = ReadOptions(),
                 PropertyCodecs const &_prop_codecs = g_default_property_codecs)
{

    std::ifstream f(_filename, std::ios::binary);
    if (!f.good()) {
        return ReadResult::CannotOpenFile;
    }
    return ovmb_read(f, _mesh, _options, _prop_codecs);
}

} // namespace OpenVolumeMesh::IO

#ifndef OVM_DO_NOT_INCLUDE_FILE_READER_IMPL
#  include <OpenVolumeMesh/IO/ovmb_read_impl.hh>
#endif
