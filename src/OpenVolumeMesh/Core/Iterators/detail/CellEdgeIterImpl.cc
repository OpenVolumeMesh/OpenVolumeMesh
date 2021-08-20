#include <OpenVolumeMesh/Core/Iterators/detail/CellEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

CellEdgeIterImpl::CellEdgeIterImpl(const CellHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    for (CellHalfEdgeIter che_iter = _mesh->che_iter(_ref_h); che_iter.valid(); ++che_iter) {
        edges_.push_back(_mesh->edge_handle(*che_iter));
    }

    // Remove all duplicate entries
    std::sort(edges_.begin(), edges_.end());
    edges_.resize(std::unique(edges_.begin(), edges_.end()) - edges_.begin());

    BaseIter::valid(edges_.size() > 0);
    if (BaseIter::valid()) {
        BaseIter::cur_handle(edges_[cur_index_]);
    }
}


CellEdgeIterImpl& CellEdgeIterImpl::operator--() {
    if (cur_index_ == 0) {
        cur_index_ = edges_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(edges_[cur_index_]);
    return *this;
}

CellEdgeIterImpl& CellEdgeIterImpl::operator++() {
    ++cur_index_;
    if (cur_index_ >= edges_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(edges_[cur_index_]);
    return *this;
}



} // namespace OpenVolumeMesh::detail
