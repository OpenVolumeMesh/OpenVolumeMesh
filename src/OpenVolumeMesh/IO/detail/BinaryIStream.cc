#include <OpenVolumeMesh/IO/detail/BinaryIStream.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>

namespace OpenVolumeMesh::IO::detail {

static uint64_t stream_length(std::istream &_s) {
    auto orig = _s.tellg();
    _s.seekg(0, std::ios_base::end);
    auto res =_s.tellg();
    _s.seekg(orig, std::ios_base::beg);
    return res;
}

BinaryIStream::BinaryIStream(std::istream &_s)
    : BinaryIStream(_s, stream_length(_s))
{
}

BinaryIStream::BinaryIStream(std::istream &_s, uint64_t _size)
    : s_(_s)
    , size_(_size)
{
}

Decoder BinaryIStream::make_decoder(size_t n)
{
    if (remaining_bytes() < n) {
        throw parse_error("make_reader: not enough bytes left.");
    }
    std::vector<uint8_t> vec(n);
    s_.read(reinterpret_cast<char*>(vec.data()), n);
    pos_ += n;
    return Decoder(std::move(vec));
}

}
