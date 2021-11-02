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
#include <string>


namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT Decoder {
public:
    Decoder(std::vector<uint8_t> _data)
        : data_(std::move(_data))
        , cur_(data_.data())
        , end_(cur_ + data_.size())
    {}
public:
// file position handling:
    bool finished() const {return cur_ == end_;}
    inline size_t size() const {return data_.size();};
    size_t remaining_bytes() const;
    inline size_t pos() const {return cur_ - data_.data();}
    void seek(size_t off);
    inline void skip() {seek(size());};

    void need(size_t n);

// basic types, WARNING: length unchecked!
    uint8_t  u8();
    uint16_t u16();
    uint32_t u32();
    uint64_t u64();
    double   dbl();
    float    flt();
    inline void read(uint8_t  &v) {v = u8();}
    inline void read(uint16_t &v) {v = u16();}
    inline void read(uint32_t &v) {v = u32();}
    inline void read(uint64_t &v) {v = u64();}
    inline void read(int8_t   &v) {v = u8();}
    inline void read(int16_t  &v) {v = u16();}
    inline void read(int32_t  &v) {v = u32();}
    inline void read(int64_t  &v) {v = u64();}
    inline void read(double &v)   {v = dbl();}
    inline void read(float &v)    {v = flt();}
    void read(std::string &v);

// basic types, WARNING: length unchecked!
    void     padding(uint8_t n);
    template<uint8_t N> void reserved();
    template<size_t N> void read(std::array<uint8_t, N> &arr);
    void read(uint8_t *s, size_t n);
    void read(char *s, size_t n);

    template<typename LengthT, typename Vec>
    void readVec(Vec &vec);

    std::vector<uint8_t> data_;
    uint8_t *cur_;
    uint8_t *end_;
};
template<typename LengthT, typename Vec>
void Decoder::readVec(Vec &vec)
{
    need(sizeof(LengthT));
    LengthT len = 0;
    read(len);
    need(len);
    vec.resize(len);
    read(vec.data(), len);
}


template<size_t N>
void Decoder::read(std::array<uint8_t, N> &arr)
{
    read(arr.data(), arr.size());
}

template<uint8_t N>
void Decoder::reserved()
{
    std::array<uint8_t, N> buf;
    read(buf);
    for (const auto ch: buf) {
        if (ch != 0) {
            throw parse_error("reserved entry != 0 at position " +
                              std::to_string(pos()));
        }
    }
}
} // namespace OpenVolumeMesh::IO::detail
