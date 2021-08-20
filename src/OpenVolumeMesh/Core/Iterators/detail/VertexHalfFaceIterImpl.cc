#include <OpenVolumeMesh/Core/Iterators/detail/VertexHalfFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

VertexHalfFaceIterImpl::VertexHalfFaceIterImpl(const VertexHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    for (VertexEdgeIter ve_iter = _mesh->ve_iter(_ref_h); ve_iter.valid(); ++ve_iter) {
        for (EdgeHalfFaceIter ehf_iter = _mesh->ehf_iter(*ve_iter); ehf_iter.valid(); ++ehf_iter) {
            halffaces_.push_back(*ehf_iter);
        }
    }

    // Remove all duplicate entries
    std::sort(halffaces_.begin(), halffaces_.end());
    halffaces_.resize(std::unique(halffaces_.begin(), halffaces_.end()) - halffaces_.begin());

    BaseIter::valid(halffaces_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(halffaces_[cur_index_]);
    }
}

VertexHalfFaceIterImpl& VertexHalfFaceIterImpl::operator--() {
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

VertexHalfFaceIterImpl& VertexHalfFaceIterImpl::operator++() {
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
