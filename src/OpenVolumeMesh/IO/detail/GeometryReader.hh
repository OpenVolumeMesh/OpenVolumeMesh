#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/IO/detail/Decoder.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>


namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT GeometryReaderBase {
public:
    virtual ~GeometryReaderBase() = default;
    virtual size_t dim() const = 0;
    virtual void read(Decoder & _decoder, VertexEncoding encoding, uint32_t first, uint32_t count) const = 0;
};

template<typename VecT>
class OVM_EXPORT GeometryReaderT : public GeometryReaderBase
{
public:
    GeometryReaderT(GeometryKernelT<VecT> & _geometry_kernel)
        : geometry_kernel_(_geometry_kernel)
    {}
    size_t dim() const override;
    void read(Decoder & _decoder, VertexEncoding encoding, uint32_t first, uint32_t count) const override;

private:
    GeometryKernelT<VecT> & geometry_kernel_;
};

extern template class OVM_EXPORT GeometryReaderT<Vec3f>;
extern template class OVM_EXPORT GeometryReaderT<Vec3d>;

} // namespace OpenVolumeMesh::IO::detail


#include <OpenVolumeMesh/IO/detail/GeometryReaderT_impl.hh>
