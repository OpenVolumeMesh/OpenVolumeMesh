#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/detail/internal_type_name.hh>
#include <memory>
#include <bitset>
#include <array>
#include <vector>

#include <OpenVolumeMesh/IO/PropertyCodecsT_impl.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

namespace OpenVolumeMesh::IO {
using namespace detail;

namespace Codecs {


struct BoolPropCodec {
    using T = bool;
    static void encode_one(Encoder &enc, const T &val) {
        enc.u8(val);
    }
    static void decode_one(Decoder &reader, T &val) {
        uint8_t v = reader.u8();
        if (v != 0 && v != 1) {
            throw parse_error("invalid bool encoding");
        }
        val = !!v;
    }
    static void encode_n(Encoder &enc, std::vector<T> const&vec, size_t idx_begin, size_t idx_end) {
        for (size_t i = idx_begin; i < idx_end; i += 8) {
            std::bitset<8> bitset;
            size_t n_bits = std::min(size_t(8), idx_end-i);
            for (size_t bit = 0; bit < n_bits; ++bit) {
                bitset.set(bit, vec[i+bit]);
            }
            enc.u8(static_cast<uint8_t>(bitset.to_ulong()));
        }
    }
    static void decode_n(Decoder &decoder, std::vector<T> &vec, size_t idx_begin, size_t idx_end) {
        decoder.need((idx_end-idx_begin+7)/8);
        for (size_t i = idx_begin; i < idx_end; i += 8) {
            uint8_t tmp = decoder.u8();
            std::bitset<8> bitset(tmp);
            size_t n_bits = std::min(size_t(8), idx_end-i);
            for (size_t bit = 0; bit < n_bits; ++bit) {
                vec[i+bit] = bitset[bit];
            }
        }
    }
};
}


PropertyCodecs PropertyCodecs::with_all_default_types() {
    PropertyCodecs codecs;
    codecs.add_default_types();
    codecs.add_ovm_vector_types();
    return codecs;
}

void PropertyCodecs::add_default_types()
{
    using namespace Codecs;
    register_codec<BoolPropCodec>("b");
    register_codec<SimplePropCodec<Primitive<uint8_t>>> ("u8");
    register_codec<SimplePropCodec<Primitive<uint16_t>>>("u16");
    register_codec<SimplePropCodec<Primitive<uint32_t>>>("u32");
    register_codec<SimplePropCodec<Primitive<uint64_t>>>("u64");
    register_codec<SimplePropCodec<Primitive<int8_t>>>  ("i8");
    register_codec<SimplePropCodec<Primitive<int16_t>>> ("i16");
    register_codec<SimplePropCodec<Primitive<int32_t>>> ("i32");
    register_codec<SimplePropCodec<Primitive<int64_t >>>("i64");
    register_codec<SimplePropCodec<Primitive<float>>>   ("f");
    register_codec<SimplePropCodec<Primitive<double>>>  ("d");

    register_codec<SimplePropCodec<Primitive<std::string>>>("s32");

    register_codec<SimplePropCodec<OVMHandle<VH>>> ("vh");
    register_codec<SimplePropCodec<OVMHandle<EH>>> ("eh");
    register_codec<SimplePropCodec<OVMHandle<HEH>>>("heh");
    register_codec<SimplePropCodec<OVMHandle<FH>>> ("fh");
    register_codec<SimplePropCodec<OVMHandle<HFH>>>("hfh");
    register_codec<SimplePropCodec<OVMHandle<CH>>> ("ch");
}

void PropertyCodecs::add_ovm_vector_types()
{
    register_arraylike<VectorT<double, 2>>("2d");
    register_arraylike<VectorT<double, 3>>("3d");
    register_arraylike<VectorT<double, 4>>("4d");

    register_arraylike<VectorT<float, 2>>("2f");
    register_arraylike<VectorT<float, 3>>("3f");
    register_arraylike<VectorT<float, 4>>("4f");

    register_arraylike<VectorT<uint32_t, 2>>("2u32");
    register_arraylike<VectorT<uint32_t, 3>>("3u32");
    register_arraylike<VectorT<uint32_t, 4>>("4u32");

    register_arraylike<VectorT<int32_t, 2>>("2i32");
    register_arraylike<VectorT<int32_t, 3>>("3i32");
    register_arraylike<VectorT<int32_t, 4>>("4i32");
}


const PropertyDecoderBase* PropertyCodecs::get_decoder(const std::string &ovmb_type_name) const
{
    auto it = decoders_.find(ovmb_type_name);
    if (it == decoders_.end()) return nullptr;
    return it->second.get();
}

const PropertyEncoderBase* PropertyCodecs::get_encoder(const std::string &internal_type_name) const
{
    auto it = encoders_.find(internal_type_name);
    if (it == encoders_.end()) return nullptr;
    return it->second.get();

}

const PropertyCodecs g_default_property_codecs = PropertyCodecs::with_all_default_types();


} // namespace OpenVolumeMesh::IO

