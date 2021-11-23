#pragma once

#include <OpenVolumeMesh/Config/Export.hh>
#include <string>
#include <typeinfo>

namespace OpenVolumeMesh::detail {

/// Get an internal name for a type. Important: this differs between
/// compilers and versions, do NOT use in file formats!
/// We need this in order to provide property type safety when
/// only limited RTTI support is available.
OVM_EXPORT std::string internal_type_name(std::type_info const &ti);

template<typename T>
std::string internal_type_name() {return internal_type_name(typeid(T));}

} //namespace OpenVolumeMesh::detail
