#pragma once

#include <OpenVolumeMesh/IO/PropertyCodecs.hh>
#include <OpenVolumeMesh/Core/detail/internal_type_name.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <OpenVolumeMesh/IO/detail/exceptions.hh>


namespace OpenVolumeMesh::IO {
using namespace detail;

template<typename T, typename Codec>
class PropertyEncoderT : public PropertyEncoderBase {
public:
    using PropertyEncoderBase::PropertyEncoderBase;
    void serialize_default(const PropertyStorageBase *prop, WriteBuffer&) const final;
    void serialize(const PropertyStorageBase *prop, WriteBuffer&, size_t idx_begin, size_t idx_end) const final;
};


template<typename T, typename Codec>
void PropertyEncoderT<T, Codec>::serialize_default(
        const PropertyStorageBase *prop,
        WriteBuffer &_write_buffer) const
{
    Encoder encoder(_write_buffer);
    Codec::encode_one(encoder, prop->cast_to_StorageT<T>()->def());
}

template<typename T, typename Codec>
void PropertyEncoderT<T, Codec>::serialize(
        const PropertyStorageBase *prop_base,
        WriteBuffer& _buffer,
        size_t idx_begin,
        size_t idx_end) const
{
    const PropertyStorageT<T> *prop = prop_base->cast_to_StorageT<T>();
    assert (idx_begin < prop->size());
    assert (idx_begin <= idx_end && idx_end <= prop->size());
    Encoder encoder(_buffer);
    Codec::encode_n(encoder, prop->data_vector(), idx_begin, idx_end);
}


template<typename T, typename Codec>
class PropertyDecoderT : public PropertyDecoderBase {
public:
    std::shared_ptr<PropertyStorageBase> request_property(
            ResourceManager &resman,
            EntityType type,
            std::string const &name,
            const std::vector<uint8_t> &encoded_def) const final;
    void deserialize(PropertyStorageBase*, Decoder&, size_t idx_begin, size_t idx_end) const final;
};


template<typename T, typename Codec>
std::shared_ptr<PropertyStorageBase>
PropertyDecoderT<T, Codec>::
request_property(
        ResourceManager &resman,
        EntityType type,
        std::string const &name,
        const std::vector<uint8_t> &encoded_def) const
{
    T def;
    Decoder decoder(encoded_def);
    Codec::decode_one(decoder, def);

    return entitytag_dispatch(type, [&](auto entity_tag)
    {
        auto prop = resman.request_property<T, decltype(entity_tag)>(name, def);
        resman.set_persistent(prop);
        // TODO: cast shoudl not be needed:
        return static_cast<PropertyStoragePtr<T>*>(&prop)->storage();
    });
}

template<typename T, typename Codec>
void PropertyDecoderT<T, Codec>::
deserialize(
        PropertyStorageBase* _prop_base,
        Decoder &_decoder,
        size_t idx_begin,
        size_t idx_end) const
{
    PropertyStorageT<T> *prop = _prop_base->cast_to_StorageT<T>();
    if (idx_begin > idx_end
            || idx_end > prop->size())
    {
        throw parse_error("invalid prop range");
    }
    Codec::decode_n(_decoder, prop->data_vector(), idx_begin, idx_end);
}

namespace Codecs {
template<typename Codec>
struct SimplePropCodec {
    using T = typename Codec::T;

    static void encode_one(Encoder &enc, const T &val) {
        Codec::encode(enc, val);
    }
    static void decode_one(Decoder &dec, T &val) {
        Codec::decode(dec, val);
    }
    static void decode_n(Decoder &decoder, std::vector<T> &vec, size_t idx_begin, size_t idx_end) {
        assert(idx_begin < idx_end);
        assert(idx_end <= vec.size());
        for (size_t i = idx_begin; i < idx_end; ++i) {
            Codec::decode(decoder, vec[i]);
        }
    }
    static void encode_n(Encoder &enc, std::vector<T> const&vec, size_t idx_begin, size_t idx_end) {
        assert(idx_begin < idx_end);
        assert(idx_end <= vec.size());
        for (size_t i = idx_begin; i < idx_end; ++i) {
            Codec::encode(enc, vec[i]);
        }
    }

};

template<typename _T>
struct Primitive {
    using T = _T;
    static void encode(Encoder &enc, const T &val) {
        enc.write(val);
    }
    static void decode(Decoder &reader, T &val) {
        reader.read(val);
    }
};

template<typename _T>
struct OVMHandle {
    static_assert(is_handle_v<_T>);
    using T = _T;
    static void encode(Encoder &enc, const T &val) {
        enc.write(val.idx());
    }
    static void decode(Decoder &reader, T &val) {
        reader.read(val.idx_mutable());
    }
};

template<typename _T, size_t N, typename SubCodec=Primitive<typename _T::value_type>>
struct ArrayLike {
    using T = _T;
    static void encode(Encoder &enc, const T &val)
    {
        for (size_t i = 0; i < N; ++i) {
            SubCodec::encode(enc, val[i]);
        }
    }

    static void decode(Decoder &reader, T &val)
    {
        for (size_t i = 0; i < N; ++i) {
            SubCodec::decode(reader, val[i]);
        }
    }
};

template<typename _T, size_t _Rows, size_t _Cols, typename SubCodec=Primitive<typename _T::value_type>>
struct MatrixLike {
    using T = _T;
    static void encode(Encoder &enc, const T &val)
    {
        for (size_t r = 0; r < _Rows; ++r) {
            for (size_t c = 0; c < _Cols; ++c) {
                SubCodec::encode(enc, val(r, c));
            }
        }
    }

    static void decode(Decoder &reader, T &val)
    {
        for (size_t r = 0; r < _Rows; ++r) {
            for (size_t c = 0; c < _Cols; ++c) {
                SubCodec::decode(reader, val(r, c));
            }
        }
    }
};

} // namespace Codecs

template<typename T, size_t N>
void PropertyCodecs::register_arraylike(const std::string &ovmb_type_name)
{
    using namespace Codecs;
    register_codec<SimplePropCodec<ArrayLike<T, N>>>(ovmb_type_name);
}

template<typename T, size_t Rows, size_t Cols>
void PropertyCodecs::register_matrixlike(const std::string &ovmb_type_name)
{
    using namespace Codecs;
    register_codec<SimplePropCodec<MatrixLike<T, Rows, Cols>>>(ovmb_type_name);
}

template<typename Codec>
void PropertyCodecs::register_codec(const std::string &ovmb_type_name)
{
    using T = typename Codec::T;
    // TODO: check if codec is already registered
    encoders_[OpenVolumeMesh::detail::internal_type_name<T>()] = std::make_shared<PropertyEncoderT<T, Codec>>(ovmb_type_name);
    decoders_[ovmb_type_name] = std::make_shared<PropertyDecoderT<T, Codec>>();
}



} // namespace OpenVolumeMesh::IO
