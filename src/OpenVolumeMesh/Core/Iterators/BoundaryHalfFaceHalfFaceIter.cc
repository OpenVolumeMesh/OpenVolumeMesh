#include <OpenVolumeMesh/Core/Iterators/BoundaryHalfFaceHalfFaceIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

BoundaryHalfFaceHalfFaceIter::BoundaryHalfFaceHalfFaceIter(const HalfFaceHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    // Go over all incident halfedges
//    const std::vector<HalfEdgeHandle> halfedges = _mesh->halfface(_ref_h).halfedges();
    const OpenVolumeMeshFace& halfface = _mesh->halfface(_ref_h);
    const std::vector<HalfEdgeHandle>& halfedges = halfface.halfedges();
    for(std::vector<HalfEdgeHandle>::const_iterator he_it = halfedges.begin();
            he_it != halfedges.end(); ++he_it) {

        // Get outside halffaces
        OpenVolumeMesh::HalfEdgeHalfFaceIter hehf_it = _mesh->hehf_iter(_mesh->opposite_halfedge_handle(*he_it));
        for(; hehf_it.valid(); ++hehf_it) {

            if(_mesh->is_boundary(*hehf_it)) {
                neighbor_halffaces_.push_back(*hehf_it);
                common_edges_.push_back(_mesh->edge_handle(*he_it));
            }
        }
    }

    cur_index_ = 0;
    BaseIter::valid(neighbor_halffaces_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    }
}


BoundaryHalfFaceHalfFaceIter& BoundaryHalfFaceHalfFaceIter::operator--() {

    if (cur_index_ == 0)
    {
        cur_index_ = neighbor_halffaces_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else{
        --cur_index_;
    }

    BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    return *this;
}


BoundaryHalfFaceHalfFaceIter& BoundaryHalfFaceHalfFaceIter::operator++() {

    ++cur_index_;
    if (cur_index_ == neighbor_halffaces_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(neighbor_halffaces_[cur_index_]);
    return *this;
}

} // namespace OpenVolumeMesh
