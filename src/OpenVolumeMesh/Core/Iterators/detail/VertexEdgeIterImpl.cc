#include <OpenVolumeMesh/Core/Iterators/detail/VertexEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

VertexEdgeIterImpl::VertexEdgeIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    voh_iter_(_mesh->voh_iter(_ref_h, _max_laps))
{
    BaseIter::valid(voh_iter_.valid());
    if (BaseIter::valid()) {
        BaseIter::cur_handle(mesh()->edge_handle(*voh_iter_));
    }
}

VertexEdgeIterImpl& VertexEdgeIterImpl::operator--() {
    --voh_iter_;
    BaseIter::lap(voh_iter_.lap());
    BaseIter::valid(voh_iter_.valid());
    BaseIter::cur_handle(mesh()->edge_handle(*voh_iter_));
    return *this;
}

VertexEdgeIterImpl& VertexEdgeIterImpl::operator++() {
    ++voh_iter_;
    BaseIter::lap(voh_iter_.lap());
    BaseIter::valid(voh_iter_.valid());
    BaseIter::cur_handle(mesh()->edge_handle(*voh_iter_));
    return *this;
}


} // namespace OpenVolumeMesh::detail
