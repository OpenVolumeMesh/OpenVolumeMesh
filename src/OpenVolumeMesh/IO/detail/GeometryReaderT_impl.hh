#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/IO/detail/GeometryReader.hh>
#include <OpenVolumeMesh/IO/detail/Encoder.hh>

namespace OpenVolumeMesh::IO::detail {

template<typename VecT>
size_t GeometryReaderT<VecT>::dim() const
{
   return VecT::dim();
}

template<typename VecT>
void GeometryReaderT<VecT>::read(
        Decoder & _decoder,
        VertexEncoding _vertex_encoding,
        uint32_t first, uint32_t count)
const
{
    auto read_all = [&](auto read_one)
    {
        for (size_t i = first; i < first+count; ++i) {
            auto vh = VH::from_unsigned(i);
            VecT &point = geometry_kernel_[vh];
            for (size_t d = 0; d < dim(); ++d) {
                point[d] = read_one(_decoder);
            }
        }
    };

    call_with_decoder(_vertex_encoding, read_all);
}

} // namespace OpenVolumeMesh::IO::detail
