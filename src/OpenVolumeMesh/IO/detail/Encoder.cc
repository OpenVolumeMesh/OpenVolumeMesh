
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/Core/Entities.hh>

#include <cassert>
#include <cstring>
#include <limits>
#include <iostream>

namespace OpenVolumeMesh::IO::detail {


void Encoder::u8(uint8_t val)
{
    uint8_t *out = buffer_.bytes_to_write(1);
    *out = val;
}

void Encoder::u16(uint16_t val)
{
    uint8_t *out = buffer_.bytes_to_write(2);
    out[0] =  val       & 0xff;
    out[1] = (val >> 8) & 0xff;
}

void Encoder::u32(uint32_t val)
{
    uint8_t *out = buffer_.bytes_to_write(4);
    for (int i = 0; i < 4; ++i) {
        out[i] = (val >> (8*i)) & 0xff;
    }
}

void Encoder::u64(uint64_t val)
{
    uint8_t *out = buffer_.bytes_to_write(8);
    for (int i = 0; i < 8; ++i) {
        out[i] = (val >> (8*i)) & 0xff;
    }
}

void Encoder::dbl(double val)
{
    uint8_t *out = buffer_.bytes_to_write(sizeof(double));
    std::memcpy(out, &val, sizeof(double));
}

void Encoder::flt(float val)
{
    uint8_t *out = buffer_.bytes_to_write(sizeof(float));
    std::memcpy(out, &val, sizeof(float));
}
void Encoder::padding(size_t n)
{
    uint8_t *out = buffer_.bytes_to_write(n);
    std::memset(out, 0, n);
}


void Encoder::write(const uint8_t *s, size_t n)
{
    buffer_.write(s, n);
}

void Encoder::write(const char *s, size_t n)
{
    buffer_.write(reinterpret_cast<const uint8_t*>(s), n);
}

} // namespace OpenVolumeMesh::IO::detail