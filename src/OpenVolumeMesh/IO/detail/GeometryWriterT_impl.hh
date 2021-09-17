#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/IO/detail/GeometryWriter.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>

namespace OpenVolumeMesh::IO::detail {

template<typename VecT>
size_t GeometryWriterT<VecT>::dim() const
{
   return VecT::dim();
}

template<typename VecT>
VertexEncoding GeometryWriterT<VecT>::vertex_encoding() const
{
    return VertexEncoding::Double;
}


template<typename VecT>
size_t GeometryWriterT<VecT>::elem_size() const
{
    return dim() * ::OpenVolumeMesh::IO::detail::elem_size(vertex_encoding());
}

template<typename VecT>
void GeometryWriterT<VecT>::write(WriteBuffer & _writebuf,
                                  ArraySpan const& _span) const
{

    Encoder encoder(_writebuf);
    auto end = _span.first + _span.count;
    assert(end <= geometry_kernel_.size());

    auto write_all = [&](auto write_one)
    {
        for (uint32_t i = _span.first; i < end; ++i) {
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
