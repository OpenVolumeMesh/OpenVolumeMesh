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
    virtual void write(WriteBuffer & _writebuf, ArraySpan const& _span) const = 0;
};

template<typename VecT>
class GeometryWriterT : public GeometryWriterBase
{
public:
    using VPositionProp = typename GeometryKernel<VecT>::VPositionProp;
    GeometryWriterT(VPositionProp const& _vertex_pos)
        : vertex_pos_(_vertex_pos)
    {}
    VertexEncoding vertex_encoding() const override;
    size_t dim() const override;
    size_t elem_size() const override;
    void write(WriteBuffer & _writebuf, ArraySpan const& _span) const override;

private:
    VPositionProp const& vertex_pos_;
};

extern template class OVM_EXPORT GeometryWriterT<Vec3f>;
extern template class OVM_EXPORT GeometryWriterT<Vec3d>;


} // namespace OpenVolumeMesh::IO::detail
#include <OpenVolumeMesh/IO/detail/GeometryWriterT_impl.hh>
