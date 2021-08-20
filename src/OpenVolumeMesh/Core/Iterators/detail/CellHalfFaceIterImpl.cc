#include <OpenVolumeMesh/Core/Iterators/detail/CellHalfFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

CellHalfFaceIterImpl::CellHalfFaceIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    hf_begin_(BaseIter::mesh()->cell(_ref_h).halffaces().begin()),
    hf_iter_(hf_begin_),
    hf_end_(BaseIter::mesh()->cell(_ref_h).halffaces().end())
{
    BaseIter::valid(hf_iter_ != hf_end_);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(*hf_iter_);
    }
}

CellHalfFaceIterImpl& CellHalfFaceIterImpl::operator--() {
    if (hf_iter_ == hf_begin_) {
        hf_iter_ = hf_end_;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    --hf_iter_;
    BaseIter::cur_handle(*hf_iter_);
    return *this;
}

CellHalfFaceIterImpl& CellHalfFaceIterImpl::operator++() {
    ++hf_iter_;
    if (hf_iter_ == hf_end_) {
        hf_iter_ = hf_begin_;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(*hf_iter_);
    return *this;
}

} // namespace OpenVolumeMesh::detail
