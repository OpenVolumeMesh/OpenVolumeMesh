#pragma once

#include <OpenVolumeMesh/Core/GeometryKernel.hh>
#include <OpenVolumeMesh/IO/detail/WriteBuffer.hh>
#include <OpenVolumeMesh/IO/detail/ovmb_format.hh>
#include <OpenVolumeMesh/Geometry/VectorT.hh>


namespace OpenVolumeMesh::IO::detail {

class OVM_EXPORT GeometryWriterBase {
public:
    virtual ~GeometryWriterBase() = default;
    virtual size_t dim() const = 0;
    virtual VertexEncoding vertex_encoding() const = 0;
    virtual size_t elem_size() const = 0;
    virtual void write(WriteBuffer & _writebuf, uint32_t first, uint32_t count) const = 0;
};

template<typename VecT>
class GeometryWriterT : public GeometryWriterBase
{
public:
    GeometryWriterT(GeometryKernelT<VecT> const& _geometry_kernel)
        : geometry_kernel_(_geometry_kernel)
    {}
    constexpr VertexEncoding vertex_encoding() const override;
    constexpr size_t dim() const override;
    constexpr size_t elem_size() const override;
    void write(WriteBuffer & _writebuf, uint32_t first, uint32_t count) const override;

private:
    GeometryKernelT<VecT> const& geometry_kernel_;
};

extern template class OVM_EXPORT GeometryWriterT<Vec3f>;
extern template class OVM_EXPORT GeometryWriterT<Vec3d>;

} // namespace OpenVolumeMesh::IO::detail
