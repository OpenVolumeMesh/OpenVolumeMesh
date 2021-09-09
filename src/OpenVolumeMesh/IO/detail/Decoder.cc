#include <OpenVolumeMesh/IO/detail/Decoder.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
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
    double ret = 0.;
    read(reinterpret_cast<uint8_t*>(&ret), sizeof(double));
    return ret;
}

float Decoder::flt()
{
    float ret = 0.;
    read(reinterpret_cast<uint8_t*>(&ret), sizeof(float));
    return ret;
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
       throw std::runtime_error("read beyond buffer");
   }
}

void Decoder::padding(uint8_t n)
{
    for (int i = 0; i < n; ++i) {
        if (cur_[i] != 0) {
            throw std::runtime_error("padding not 0");
        }
    }
    cur_ += n;
}


} // namespace OpenVolumeMesh::IO::detail