#pragma once

#include <map>
#include <string>
#include <OpenVolumeMesh/IO/BinaryIO.hh>
#include <OpenVolumeMesh/Core/Entities.hh>
#include <OpenVolumeMesh/Core/ResourceManager.hh>

namespace OpenVolumeMesh::IO {


#if 0
class PropertyEncoder {
    virtual bool read(StreamWriter&, size_t base, size_t count);
private:
    const std::string data_type_name_;
};

class PropertyDecoder {
    virtual bool create_prop(ResourceManager &resman, EntityType type, std::string const &name) = 0;
    virtual bool read(BufferReader&, size_t base, size_t count);
private:
};

template<typename T>
class PropertyDecoderT : public PropertyDecoder {
    bool create_prop(ResourceManager &resman, EntityType type, std::string const &name) override;
    bool read(BufferReader&, size_t base, size_t count) override;
private:
    // ref to prop

};

class PropertyEncoderFactory {
    template<typename MeshT>
    PropertyEncoder make(const MeshT &mesh, property...);
};

class PropertyDecoderFactory {
};

template<typename T>
void register_prop_codec(std::string name);
#endif


} // namespace OpenVolumeMesh::IO

#include <OpenVolumeMesh/IO/PropertySerializationT.hh>
