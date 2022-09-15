#include <OpenVolumeMesh/IO/detail/Decoder.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>
#include <OpenVolumeMesh/Core/Entities.hh>

#include <cassert>
#include <cstring>
#include <limits>
#include <iostream>

namespace OpenVolumeMesh::IO::detail {


uint8_t Decoder::u8()
{
    return *(cur_++);

}
uint16_t Decoder::u16()
{
    uint16_t res = (cur_[0] + (uint16_t(cur_[1]) << 8));
    cur_ += 2;
    return res;
}
uint32_t Decoder::u32()
{
    uint32_t res =  (cur_[0]
            + (uint32_t(cur_[1]) << 8)
            + (uint32_t(cur_[2]) << 16)
            + (uint32_t(cur_[3]) << 24));
    cur_ += 4;
    return res;
}

uint64_t Decoder::u64()
{
    uint64_t res = (cur_[0]
            + (uint64_t(cur_[1]) << 8)
            + (uint64_t(cur_[2]) << 16)
            + (uint64_t(cur_[3]) << 24)
            + (uint64_t(cur_[4]) << 32)
            + (uint64_t(cur_[5]) << 40)
            + (uint64_t(cur_[6]) << 48)
            + (uint64_t(cur_[7]) << 56)
            );
    cur_ += 8;
    return res;
}

double Decoder::dbl()
{
    uint64_t tmp {u64()};
    double val = 0.;
    static_assert(sizeof(tmp) == sizeof(val));
    // in c++20, we will be able to use bit_cast
    ::memcpy(&val, &tmp, sizeof(val));
    return val;
}

float Decoder::flt()
{
    uint32_t tmp {u32()};
    float val = 0.f;
    static_assert(sizeof(tmp) == sizeof(val));
    // in c++20, we will be able to use bit_cast
    ::memcpy(&val, &tmp, sizeof(val));
    return val;
}

void Decoder::read(std::string &v)
{
   auto len = u32();
   need(len);
   v.resize(len);
   read(v.data(), len);
}


void Decoder::read(uint8_t *s, size_t n)
{
    std::memcpy(s, cur_, n);
    cur_ += n;
}

void Decoder::read(char *s, size_t n)
{
    std::memcpy(s, cur_, n);
    cur_ += n;
}


size_t Decoder::remaining_bytes() const {
    assert (cur_ <= end_);
    return end_ - cur_;
}

void Decoder::seek(size_t off) {
    assert(off <= size());
    cur_ = data_.data() + off;
}

void Decoder::need(size_t n)
{
   if (remaining_bytes() < n) {
       throw parse_error("read beyond buffer");
   }
}

void Decoder::padding(uint8_t n)
{
    for (int i = 0; i < n; ++i) {
        if (cur_[i] != 0) {
            throw parse_error("padding not 0");
        }
    }
    cur_ += n;
}


} // namespace OpenVolumeMesh::IO::detail
