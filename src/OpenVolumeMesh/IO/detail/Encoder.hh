#pragma once

#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>
#include <istream>
#include <ostream>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>
#include <limits>
#include <cstring>
#include <string>


namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT Encoder {
public:
    Encoder(WriteBuffer &_buffer)
        : buffer_(_buffer)
    {}
public:
    //WriteBuffer &buffer() {return buffer_;}
    // basic types
    void u8(uint8_t);
    void u16(uint16_t);
    void u32(uint32_t);

    void u64(uint64_t);
    void dbl(double);
    void flt(float);
    inline void write(uint8_t  v) {u8(v);}
    inline void write(uint16_t v) {u16(v);}
    inline void write(uint32_t v) {u32(v);}
    inline void write(uint64_t v) {u64(v);}
    inline void write(int8_t   v) {u8(v);}
    inline void write(int16_t  v) {u16(v);}
    inline void write(int32_t  v) {u32(v);}
    inline void write(int64_t  v) {u64(v);}
    inline void write(float v)    {flt(v);}
    inline void write(double v)   {dbl(v);}
    void write(std::string const&v);

    // helpers
    void padding(size_t n);
    template<uint8_t N> void reserved();
    template<size_t N> void write(std::array<uint8_t, N> const &arr);
    void write(const uint8_t *s, size_t n);
    void write(const char *s, size_t n);

    template<typename LengthT, typename Vec>
    void writeVec(Vec const &vec);
private:
    WriteBuffer &buffer_;
};


template<size_t N>
void Encoder::write(const std::array<uint8_t, N> &arr)
{
    write(arr.data(), arr.size());
}

OVM_EXPORT extern const std::array<uint8_t, 256> zero_buf;
template<uint8_t N>
void Encoder::reserved()
{
    uint8_t *out = buffer_.bytes_to_write(N);
    std::memset(out, 0, N);
}

template<typename LengthT, typename Vec>
void Encoder::writeVec(Vec const &vec)
{
    if (vec.size() > std::numeric_limits<LengthT>::max()) {
        throw write_error("vector too long for length data type: " +
                          std::to_string(vec.size()));
        // TODO: maybe just truncate?
    }
    auto len = static_cast<LengthT>(vec.size());
    write(len);
    write(vec.data(), len);
}
} // namespace OpenVolumeMesh::IO::detail
