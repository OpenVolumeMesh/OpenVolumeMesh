#pragma once

#include <map>
#include <string>
#include <OpenVolumeMesh/IO/BinaryIO.hh>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>
#include <OpenVolumeMesh/Core/PropertyPtr.hh>
#include <OpenVolumeMesh/Core/PropertyStorageT.hh>
#include <OpenVolumeMesh/Core/EntityUtils.hh>
#include <any>

namespace OpenVolumeMesh::IO {


class PropertyEncoderBase {
public:
    PropertyEncoderBase(std::string _ovmb_type_name)
        : ovmb_type_name_(std::move(_ovmb_type_name))
    {}
    virtual ~PropertyEncoderBase() = default;
    //virtual void get_prop(ResourceManager &resman, EntityType type, std::string const &name) = 0;
    virtual void serialize_default(PropertyStorageBase *prop, StreamWriter&) = 0;
    virtual void serialize(PropertyStorageBase *prop, StreamWriter&, size_t idx_begin, size_t idx_end) = 0;
    std::string const &ovmb_type_name() {return ovmb_type_name_;};
private:
    std::string ovmb_type_name_;
};

class PropertyDecoderBase {
public:
    virtual ~PropertyDecoderBase() = default;
    virtual void request_prop(
            ResourceManager &resman,
            EntityType type,
            std::string const &name,
            const std::vector<uint8_t> &encoded_def) = 0;
    virtual void deserialize(BufferReader&, size_t idx_begin, size_t idx_end) = 0;
};

template<typename T, typename Codec>
void register_prop_codec(std::string const &ovmb_type_name);


} // namespace OpenVolumeMesh::IO

#include <OpenVolumeMesh/IO/PropertySerializationT.hh>
