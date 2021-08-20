#include <OpenVolumeMesh/Core/Iterators/detail/HalfFaceEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

HalfFaceEdgeIterImpl::HalfFaceEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    BaseIter::valid(_ref_h.is_valid() && _mesh->halfface(_ref_h).halfedges().size() > 0);
    if (BaseIter::valid()) {
        HalfEdgeHandle he = _mesh->halfface(_ref_h).halfedges()[cur_index_];
        BaseIter::cur_handle(_mesh->edge_handle(he));
    }
}

HalfFaceEdgeIterImpl& HalfFaceEdgeIterImpl::operator--() {
    const std::vector<HalfEdgeHandle> halfedges =
        mesh()->halfface(ref_handle()).halfedges();
    if (cur_index_ == 0) {
        cur_index_ = halfedges.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(mesh()->edge_handle(halfedges[cur_index_]));
    return *this;
}

HalfFaceEdgeIterImpl& HalfFaceEdgeIterImpl::operator++() {
    const std::vector<HalfEdgeHandle> halfedges =
        mesh()->halfface(ref_handle()).halfedges();
    ++cur_index_;
    if (cur_index_ >= halfedges.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(mesh()->edge_handle(halfedges[cur_index_]));
    return *this;
}


} // namespace OpenVolumeMesh::detail
