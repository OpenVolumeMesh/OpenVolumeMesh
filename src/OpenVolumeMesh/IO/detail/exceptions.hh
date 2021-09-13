#pragma once
#include <OpenVolumeMesh/Config/Export.hh>
#include <stdexcept>

namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT io_error : public std::runtime_error
{
    using std::runtime_error::runtime_error;
};

class OVM_EXPORT parse_error : public io_error
{
    using io_error::io_error;
};

class OVM_EXPORT write_error : public io_error
{
    using io_error::io_error;
};


} // namespace OpenVolumeMesh::IO::detail
