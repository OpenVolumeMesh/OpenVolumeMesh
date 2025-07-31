#pragma once

#include <OpenVolumeMesh/IO/enums.hh>
#include <OpenVolumeMesh/IO/ReadOptions.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReader.hh>
#include <istream>
#include <fstream>
#include <memory>
#include <filesystem>

namespace OpenVolumeMesh::IO {

std::unique_ptr<detail::BinaryFileReader>
inline make_ovmb_reader(std::istream & _istream,
            const ReadOptions &_options,
            PropertyCodecs const &_prop_codecs)
{
    return std::make_unique<detail::BinaryFileReader>(
                _istream,
                _options,
                _prop_codecs);

}

/// The _istream MUST be opened in using std::ios::binary!
template<typename MeshT>
ReadResult ovmb_read(std::istream &_istream,
                     MeshT &_mesh,
                     ReadOptions _options = ReadOptions(),
                     PropertyCodecs const &_prop_codecs = g_default_property_codecs)
{
    static_assert(std::is_base_of_v<TopologyKernel, MeshT>);
    auto reader = make_ovmb_reader(_istream, _options, _prop_codecs);
    auto result =  reader->read_file(_mesh);
#ifndef NDEBUG
    if (result != ReadResult::Ok) {
        std::cerr << "Error reading ovmb file: "
                  << to_string(result)
                  << " ("
                  << reader->get_error_msg()
                  << ")"
                  << std::endl;
    }
#endif
    return result;
}

template<typename MeshT>
ReadResult ovmb_read(std::filesystem::path const&_filename,
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
