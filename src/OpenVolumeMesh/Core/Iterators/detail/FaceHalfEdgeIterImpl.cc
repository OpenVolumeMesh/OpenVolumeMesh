#include <OpenVolumeMesh/Core/Iterators/detail/FaceHalfEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

FaceHalfEdgeIterImpl::FaceHalfEdgeIterImpl(
    const FaceHandle& _ref_h,
    const TopologyKernel* _mesh,
    int _max_laps)
        : BaseIter(_mesh, _ref_h, _max_laps)
        , halfedges_(_mesh->face(_ref_h).halfedges())
{
    BaseIter::cur_handle(halfedges_[cur_index_]);
}

FaceHalfEdgeIterImpl& FaceHalfEdgeIterImpl::operator++()
{
    cur_index_++;
    if (cur_index_ >= halfedges_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_) {
            BaseIter::valid(false);
        }
    }
    BaseIter::cur_handle(halfedges_[cur_index_]);
    return *this;
}

FaceHalfEdgeIterImpl& FaceHalfEdgeIterImpl::operator--()
{
    if (cur_index_ == 0) {
        cur_index_ = halfedges_.size();
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    --cur_index_;
    BaseIter::cur_handle(halfedges_[cur_index_]);
    return *this;
}


} // namespace OpenVolumeMesh::detail
