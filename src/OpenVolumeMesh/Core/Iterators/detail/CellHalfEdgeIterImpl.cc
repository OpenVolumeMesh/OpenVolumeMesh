#include <OpenVolumeMesh/Core/Iterators/detail/CellHalfEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

CellHalfEdgeIterImpl::CellHalfEdgeIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    for (CellHalfFaceIter chf_iter =_mesh->chf_iter(_ref_h); chf_iter.valid(); ++chf_iter) {
        for (HalfFaceHalfEdgeIter hfhe_iter =_mesh->hfhe_iter(*chf_iter); hfhe_iter.valid(); ++hfhe_iter) {
            halfedges_.push_back(*hfhe_iter);
        }
    }
    BaseIter::valid(halfedges_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(halfedges_[cur_index_]);
    }
}

CellHalfEdgeIterImpl& CellHalfEdgeIterImpl::operator--() {
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

CellHalfEdgeIterImpl& CellHalfEdgeIterImpl::operator++() {
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
