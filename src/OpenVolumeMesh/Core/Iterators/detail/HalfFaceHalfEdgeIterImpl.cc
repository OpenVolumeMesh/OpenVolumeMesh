#include <OpenVolumeMesh/Core/Iterators/detail/HalfFaceHalfEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

HalfFaceHalfEdgeIterImpl::HalfFaceHalfEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    halfedges_(_mesh->halfface(_ref_h).halfedges()),
    cur_index_(0)
{
    assert(halfedges_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(halfedges_[cur_index_]);
    }
}

HalfFaceHalfEdgeIterImpl& HalfFaceHalfEdgeIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = halfedges_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(halfedges_[cur_index_]);
    return *this;
}

HalfFaceHalfEdgeIterImpl& HalfFaceHalfEdgeIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= halfedges_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(halfedges_[cur_index_]);
    return *this;
}



} // namespace OpenVolumeMesh::detail
