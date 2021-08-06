#pragma once
#include <OpenVolumeMesh/IO/PropertySerialization.hh>


namespace OpenVolumeMesh::IO {

/// public (filename) name -> decoder
extern std::map<std::string, std::unique_ptr<PropertyDecoderBase>> property_dec;

/// internal type name -> encoder
extern std::map<std::string, std::unique_ptr<PropertyEncoderBase>> property_enc;

template<typename T, typename Codec>
class PropertyEncoderT : public PropertyEncoderBase {
public:
    //void get_prop(ResourceManager &resman, EntityType type, std::string const &name) override;
    void serialize(PropertyStorageBase *prop, StreamWriter&, size_t idx_begin, size_t idx_end) override;
private:
    std::unique_ptr<PropertyStoragePtr<T>> prop_;
};


template<typename T, typename Codec>
class PropertyDecoderT : public PropertyDecoderBase {
public:
    void request_prop(ResourceManager &resman,
                      EntityType type,
                      std::string const &name,
                      const std::vector<uint8_t> &encoded_def) override;
    void deserialize(BufferReader&, size_t idx_begin, size_t idx_end) override;
private:
    std::unique_ptr<PropertyStoragePtr<T>> prop_;
};



template<typename T, typename Codec>
void PropertyDecoderT<T, Codec>::request_prop(
        ResourceManager &resman,
        EntityType type,
        const std::string &name,
        std::vector<uint8_t> const &encoded_def)
{
    T def = Codec::decode(encoded_def);

    prop_ = entitytag_dispatch(type, [&](auto entity_tag)
    {
        return resman.request_property<T, decltype(entity_tag)>(name, def);
    });
}

template<typename T, typename Codec>
void PropertyDecoderT<T, Codec>::deserialize(BufferReader &reader, size_t idx_begin, size_t idx_end)
{
    assert (idx_begin >=0 && idx_begin < prop_.size());
    assert (idx_begin <= idx_end && idx_end <= prop_.size());
    for (size_t idx = idx_begin; idx < idx_end; ++idx) {
        prop_[idx] = Codec::decode_one(reader);
    }
}

template<typename T, typename Codec>
void PropertyEncoderT<T, Codec>::serialize(
        PropertyStorageBase *prop_base,
        StreamWriter& writer,
        size_t idx_begin,
        size_t idx_end)
{
    PropertyStorageT<T> *prop = prop_base->cast_to_StorageT<T>();
    assert (idx_begin >=0 && idx_begin < prop->size());
    assert (idx_begin <= idx_end && idx_end <= prop->size());
    for (size_t idx = idx_begin; idx < idx_end; ++idx) {
        Codec::encode_one(writer,  (*prop)[idx]);
    }
}


template<typename T, typename Codec>
void register_prop_codec(std::string &public_name)
{
    property_enc[get_type_name<T>()] = std::make_unique<PropertyEncoderT<T, Codec>>();
    property_dec[public_name] = std::make_unique<PropertyDecoderT<T, Codec>>();
}

} // namespace OpenVolumeMesh::IO
