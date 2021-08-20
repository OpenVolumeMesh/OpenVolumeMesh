#include <OpenVolumeMesh/Core/Iterators/detail/VertexIHalfEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {


VertexIHalfEdgeIterImpl::VertexIHalfEdgeIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    voh_iter_(_mesh->voh_iter(_ref_h, _max_laps))
{
    BaseIter::valid(voh_iter_.valid());
    if (BaseIter::valid()) {
        BaseIter::cur_handle(mesh()->opposite_halfedge_handle(*voh_iter_));
    }
}

VertexIHalfEdgeIterImpl& VertexIHalfEdgeIterImpl::operator--() {
    --voh_iter_;
    BaseIter::lap(voh_iter_.lap());
    BaseIter::valid(voh_iter_.valid());
    BaseIter::cur_handle(mesh()->opposite_halfedge_handle(*voh_iter_));
    return *this;
}

VertexIHalfEdgeIterImpl& VertexIHalfEdgeIterImpl::operator++() {
    ++voh_iter_;
    BaseIter::lap(voh_iter_.lap());
    BaseIter::valid(voh_iter_.valid());
    BaseIter::cur_handle(mesh()->opposite_halfedge_handle(*voh_iter_));
    return *this;
}


} // namespace OpenVolumeMesh::detail
