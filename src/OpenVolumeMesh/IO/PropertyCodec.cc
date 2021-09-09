#include <OpenVolumeMesh/IO/PropertyCodec.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/detail/internal_type_name.hh>
#include <memory>
#include <bitset>
#include <array>
#include <vector>

#include <OpenVolumeMesh/IO/PropertyCodecT_impl.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>

namespace OpenVolumeMesh::IO {
using namespace detail;

namespace Codecs {


struct BoolPropCodec {
    using T = bool;
    static void encode_one(Encoder &enc, const T &val) {
        uint8_t v = val ? 1 : 0;
        enc.write(v);
    }
    static void decode_one(Decoder &reader, T &val) {
        uint8_t v = val ? 1 : 0;
        reader.read(v);
        if (v != 0 && v != 1) {
            throw std::runtime_error("invalid bool encoding");
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
            uint8_t tmp = static_cast<uint8_t>(bitset.to_ulong());
            enc.write(tmp);
        }
    }
    static void decode_n(Decoder &decoder, std::vector<T> &vec, size_t idx_begin, size_t idx_end) {
        decoder.need((idx_end-idx_begin+7)/8);
        for (size_t i = idx_begin; i < idx_end; i += 8) {
            uint8_t tmp = 0;
            decoder.read(tmp);
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
    register_type<bool, BoolPropCodec>("b");
    register_type<uint8_t,  SimplePropCodec<Primitive<uint8_t>>> ("u8");
    register_type<uint16_t, SimplePropCodec<Primitive<uint16_t>>>("u16");
    register_type<uint32_t, SimplePropCodec<Primitive<uint32_t>>>("u32");
    register_type<uint64_t, SimplePropCodec<Primitive<uint64_t>>>("u64");
    register_type<int8_t,   SimplePropCodec<Primitive<int8_t>>>  ("i8");
    register_type<int16_t,  SimplePropCodec<Primitive<int16_t>>> ("i16");
    register_type<int32_t,  SimplePropCodec<Primitive<int32_t>>> ("i32");
    register_type<int64_t,  SimplePropCodec<Primitive<int64_t >>>("i64");
    register_type<float,    SimplePropCodec<Primitive<float>>>   ("f");
    register_type<double,   SimplePropCodec<Primitive<double>>>  ("d");

    register_type<VH,       SimplePropCodec<OVMHandle<VH>>>      ("vh");
    register_type<EH,       SimplePropCodec<OVMHandle<EH>>>      ("eh");
    register_type<HEH,      SimplePropCodec<OVMHandle<HEH>>>     ("heh");
    register_type<FH,       SimplePropCodec<OVMHandle<FH>>>      ("fh");
    register_type<HFH,      SimplePropCodec<OVMHandle<HFH>>>     ("hfh");
    register_type<CH,       SimplePropCodec<OVMHandle<CH>>>      ("ch");
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

template<typename T, size_t N>
void PropertyCodecs::register_arraylike(const std::string &ovmb_type_name)
{
    using namespace Codecs;
    register_type<T, SimplePropCodec<ArrayLike<T, N>>>(ovmb_type_name);
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

template<typename T, typename Codec>
void PropertyCodecs::register_type(const std::string &ovmb_type_name)
{
    // TODO: check if codec is already registered
    encoders_[OpenVolumeMesh::detail::internal_type_name<T>()] = std::make_unique<PropertyEncoderT<T, Codec>>(ovmb_type_name);
    decoders_[ovmb_type_name] = std::move(std::make_unique<PropertyDecoderT<T, Codec>>());
}

const PropertyCodecs g_default_property_codecs = PropertyCodecs::with_all_default_types();


} // namespace OpenVolumeMesh::IO

