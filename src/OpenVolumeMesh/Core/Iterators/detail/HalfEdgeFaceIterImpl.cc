#include <OpenVolumeMesh/Core/Iterators/detail/HalfEdgeFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

HalfEdgeFaceIterImpl::HalfEdgeFaceIterImpl(const HalfEdgeHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    for (HalfEdgeHalfFaceIter hehf_iter = _mesh->hehf_iter(_ref_h); hehf_iter.valid(); ++hehf_iter) {
        faces_.push_back(_mesh->face_handle(*hehf_iter));
    }

    // Remove all duplicate entries
    std::sort(faces_.begin(), faces_.end());
    faces_.resize(std::unique(faces_.begin(), faces_.end()) - faces_.begin());

    BaseIter::valid(faces_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(faces_[cur_index_]);
    }
}

HalfEdgeFaceIterImpl& HalfEdgeFaceIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = faces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(faces_[cur_index_]);
    return *this;
}

HalfEdgeFaceIterImpl& HalfEdgeFaceIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= faces_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(faces_[cur_index_]);
    return *this;
}



} // namespace OpenVolumeMesh::detail
