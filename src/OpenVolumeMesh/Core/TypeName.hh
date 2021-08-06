#pragma once

#include <string>
#include <typeinfo>

/// Get an internal name for a type. Important: this differs between
/// compilers and versions, do NOT use in file formats!
/// We need this in order to provide property type safety when
/// only limited RTTI support is available.
std::string get_type_name(std::type_info const &ti);

template<typename T>
std::string get_type_name() {return get_type_name(typeid(T));}
