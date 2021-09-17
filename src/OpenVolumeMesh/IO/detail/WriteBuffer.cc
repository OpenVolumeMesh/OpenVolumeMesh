#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <cassert>
#include <iostream>

namespace OpenVolumeMesh::IO::detail {

void WriteBuffer::need(size_t n)
{
    size_t remaining = data_.size() - pos_;
    if (n > remaining) {
        data_.resize(pos_ + n);
        assert(n == data_.size() - pos_);
    }
}
void WriteBuffer::write_to_stream(std::ostream &s)
{
    s.write(reinterpret_cast<const char*>(data_.data()), pos_);
    reset();
}

void WriteBuffer::reset()
{
    pos_ = 0;
}

std::vector<uint8_t> WriteBuffer::vec() const {
    return {data_.begin(), data_.begin() + pos_};
}

void WriteBuffer::write(const uint8_t *s, size_t n)
{
    need(n);
    std::copy(s, s + n, &data_[pos_]);
    pos_ += n;
}

uint8_t* WriteBuffer::bytes_to_write(size_t n) {
    need(n);
    uint8_t* ret = &data_[pos_];
    pos_ += n;
    return ret;
}

} // namespace OpenVolumeMesh::IO::detail
