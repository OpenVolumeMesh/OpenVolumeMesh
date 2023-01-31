#pragma once

#include <map>
#include <string>
#include <optional>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <OpenVolumeMesh/Core/Properties/PropertyPtr.hh>
#include <OpenVolumeMesh/Core/EntityUtils.hh>

namespace OpenVolumeMesh::IO {

class OVM_EXPORT PropertyDecoderBase {
public:
    virtual ~PropertyDecoderBase() = default;
    virtual std::shared_ptr<PropertyStorageBase> request_property(
            ResourceManager &resman,
            EntityType type,
            std::string const &name,
            const std::vector<uint8_t> &encoded_def) const = 0;
    virtual void deserialize(PropertyStorageBase*, detail::Decoder&, size_t idx_begin, size_t idx_end) const = 0;
};

class OVM_EXPORT PropertyEncoderBase {
public:
    PropertyEncoderBase(std::string _ovmb_type_name)
        : ovmb_type_name_(std::move(_ovmb_type_name))
    {}
    virtual ~PropertyEncoderBase() = default;
    virtual void serialize_default(const PropertyStorageBase *prop, detail::WriteBuffer&) const = 0;
    virtual void serialize(const PropertyStorageBase *prop, detail::WriteBuffer&, size_t idx_begin, size_t idx_end) const = 0;
    std::string const &ovmb_type_name() const {return ovmb_type_name_;};
private:
    std::string ovmb_type_name_;
};


class OVM_EXPORT PropertyCodecs {
public:
    PropertyCodecs() = default;
    static PropertyCodecs with_all_default_types();
    void add_default_types();
    void add_ovm_vector_types();

    /// Register a data type for property (de)serialisation.
    /// \tparam Codec          a type with static {en,de}code_{one,n} methods
    /// \param ovmb_type_name  the name as used in the ovmb file format
    template<typename Codec>
    void register_codec(std::string const &ovmb_type_name);

    template<typename T, size_t N = T::size()>
    void register_arraylike(std::string const &ovmb_type_name);
    template<typename T, size_t Rows, size_t Cols>
    void register_matrixlike(std::string const &ovmb_type_name);

    const PropertyDecoderBase* get_decoder(std::string const &ovmb_type_name) const;
    const PropertyEncoderBase* get_encoder(std::string const &internal_type_name) const;
private:
    // TODO: add namespaces for types, to never mix up user types and OVM types
    // NB: a unique_ptr would suffice, however MSVC gives compile errors for map with unique_ptr as value type in a DLL build
    std::map<std::string, std::shared_ptr<PropertyDecoderBase>> decoders_;
    std::map<std::string, std::shared_ptr<PropertyEncoderBase>> encoders_;
};

OVM_EXPORT extern const PropertyCodecs g_default_property_codecs;


} // namespace OpenVolumeMesh::IO

