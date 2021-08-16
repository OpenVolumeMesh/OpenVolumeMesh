
#pragma once
#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/Core/detail/internal_type_name.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>


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
        throw std::runtime_error("invalid prop range");
    }
    Codec::decode_n(_decoder, prop->data_vector(), idx_begin, idx_end);
}


} // namespace OpenVolumeMesh::IO
