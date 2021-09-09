#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/IO/detail/GeometryWriter.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>

namespace OpenVolumeMesh::IO::detail {

template<typename VecT>
constexpr size_t GeometryWriterT<VecT>::dim() const
{
   return VecT::dim();
}

template<typename VecT>
constexpr VertexEncoding GeometryWriterT<VecT>::vertex_encoding() const
{
    return VertexEncoding::Double;
}


template<typename VecT>
constexpr size_t GeometryWriterT<VecT>::elem_size() const
{
    return dim() * ::OpenVolumeMesh::IO::detail::elem_size(vertex_encoding());
}

template<typename VecT>
void GeometryWriterT<VecT>::write(WriteBuffer &_writebuf, uint32_t first, uint32_t count) const
{

    Encoder encoder(_writebuf);
    auto end = first + count;
    assert(end <= geometry_kernel_.size());

    auto write_all = [&](auto write_one)
    {
        for (uint32_t i = first; i < end; ++i) {
            auto vh = VertexHandle::from_unsigned(i);
            const auto &pos = geometry_kernel_[vh];
            for (size_t dim = 0; dim < this->dim(); ++dim) {
                write_one(encoder, pos[dim]);
            }
        }
    };
    call_with_encoder(vertex_encoding(), write_all);
}

} // namespace OpenVolumeMesh::IO::detail
