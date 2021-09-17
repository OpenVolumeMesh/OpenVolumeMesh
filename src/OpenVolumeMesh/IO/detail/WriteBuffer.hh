#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <ostream>
#include <vector>

namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT WriteBuffer {
public:
    void reset();
    std::vector<uint8_t> vec() const;;

    void write_to_stream(std::ostream &s);
    size_t size() const {return pos_;}
    size_t allocated_size() const {return data_.capacity();}

    /// reserve additional `n` bytes
    void need(size_t n);

    void write(const uint8_t *s, size_t n);
    uint8_t* bytes_to_write(size_t n);
private:
    std::vector<uint8_t> data_;
    size_t pos_ = 0;
};

} // namespace OpenVolumeMesh::IO::detail
