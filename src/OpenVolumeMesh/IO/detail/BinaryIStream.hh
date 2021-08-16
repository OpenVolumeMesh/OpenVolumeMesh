#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <istream>
#include <ostream>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>
#include <limits>


namespace OpenVolumeMesh::IO::detail {
class Decoder;

class OVM_EXPORT BinaryIStream {
public:
    explicit BinaryIStream(std::istream &_s);
    BinaryIStream(std::istream &_s,
                 uint64_t _size);
    uint64_t remaining_bytes() const {
        return size_ - pos_;
    }
    /// sub-readers share their istream; do not interleave use!
    Decoder make_decoder(size_t n);
private:
    std::istream &s_;
    uint64_t size_;
    uint64_t pos_ = 0;
};

} // namespace OpenVolumeMesh::IO::detail
