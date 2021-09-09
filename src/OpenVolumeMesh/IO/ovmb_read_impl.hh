#include <OpenVolumeMesh/IO/ovmb_read.hh>
#include <OpenVolumeMesh/IO/PropertyCodec.hh>
#include <OpenVolumeMesh/IO/detail/GeometryReader.hh>
#include <OpenVolumeMesh/IO/detail/BinaryFileReader.hh>


namespace OpenVolumeMesh::IO {

std::unique_ptr<detail::BinaryFileReader>
make_ovmb_reader(std::istream & _istream,
            const ReadOptions &_options,
            PropertyCodecs const &_prop_codecs)
{
    return std::make_unique<detail::BinaryFileReader>(
                _istream,
                _options,
                _prop_codecs);

}

template<typename MeshT>
ReadResult ovmb_read(std::istream &_istream,
                     MeshT &_mesh,
                     ReadOptions _options,
                     PropertyCodecs const &_prop_codecs)
{
    auto reader = make_ovmb_reader(_istream, _options, _prop_codecs);
    return reader->read_file(_mesh);
}

} // namespace OpenVolumeMesh::IO
