#include <OpenVolumeMesh/Core/Iterators/detail/CellFaceIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

CellFaceIterImpl::CellFaceIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    hf_iter_(BaseIter::mesh()->cell(_ref_h).halffaces().begin())
{
    BaseIter::valid(hf_iter_ != BaseIter::mesh()->cell(_ref_h).halffaces().end());
    if (BaseIter::valid()) {
        BaseIter::cur_handle(BaseIter::mesh()->face_handle(*hf_iter_));
    }
}

CellFaceIterImpl& CellFaceIterImpl::operator--() {
    const std::vector<HalfFaceHandle>& halffaces =
        BaseIter::mesh()->cell(ref_handle_).halffaces();
    if (hf_iter_ == halffaces.begin()) {
        hf_iter_ = halffaces.end();
        --lap_;
        if (lap_ < 0) {
            BaseIter::valid(false);
            return *this;
        }
    }
    --hf_iter_;
    BaseIter::cur_handle(BaseIter::mesh()->face_handle(*hf_iter_));
    return *this;
}

CellFaceIterImpl& CellFaceIterImpl::operator++() {
    ++hf_iter_;
    const std::vector<HalfFaceHandle>& halffaces =
        BaseIter::mesh()->cell(ref_handle_).halffaces();
    if (hf_iter_ == halffaces.end()) {
        hf_iter_ = halffaces.begin();
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(BaseIter::mesh()->face_handle(*hf_iter_));
    return *this;
}


} // namespace OpenVolumeMesh::detail
