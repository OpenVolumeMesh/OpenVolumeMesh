#include <OpenVolumeMesh/IO/PropertySerialization.hh>
#include <OpenVolumeMesh/Config/Export.hh>

namespace OpenVolumeMesh::IO {

/// public (filename) type name -> serializer
std::map<std::string, std::unique_ptr<PropertyDecoderBase>> property_dec;

/// internal type name -> serializer
std::map<std::string, std::unique_ptr<PropertyEncoderBase>> property_enc;

} // namespace OpenVolumeMesh::IO

