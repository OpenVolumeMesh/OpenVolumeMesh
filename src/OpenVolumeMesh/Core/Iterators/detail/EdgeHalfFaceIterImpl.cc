#include <OpenVolumeMesh/Core/Iterators/detail/EdgeHalfFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

EdgeHalfFaceIterImpl::EdgeHalfFaceIterImpl(const EdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    const HalfEdgeHandle he = _mesh->halfedge_handle(_ref_h, 0);
    for (HalfEdgeHalfFaceIter hehf_iter = _mesh->hehf_iter(he); hehf_iter.valid(); ++hehf_iter) {
        halffaces_.push_back(*hehf_iter);
        halffaces_.push_back(_mesh->opposite_halfface_handle(*hehf_iter));
    }

    BaseIter::valid(halffaces_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(halffaces_[cur_index_]);
    }
}

EdgeHalfFaceIterImpl& EdgeHalfFaceIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = halffaces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(halffaces_[cur_index_]);
    return *this;
}

EdgeHalfFaceIterImpl& EdgeHalfFaceIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= halffaces_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(halffaces_[cur_index_]);
    return *this;
}

} // namespace OpenVolumeMesh::detail
