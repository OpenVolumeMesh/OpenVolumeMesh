#include <OpenVolumeMesh/Core/Iterators/HalfEdgeHalfFaceIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

HalfEdgeHalfFaceIter::HalfEdgeHalfFaceIter(const HalfEdgeHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

    if(!_mesh->has_edge_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) {
        BaseIter::valid(false);
    }

    if(BaseIter::valid()) {
        if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[_ref_h].size()) {
            BaseIter::valid(false);
        }
    }

    if(BaseIter::valid()) {
        BaseIter::cur_handle((
                BaseIter::mesh()->incident_hfs_per_he_[_ref_h])[cur_index_]);
    }
}


HalfEdgeHalfFaceIter& HalfEdgeHalfFaceIter::operator--() {

    size_t n_outgoing_halffaces = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()].size();

    if (cur_index_ == 0) {
        cur_index_ = n_outgoing_halffaces-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle((BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()])[cur_index_]);

    return *this;
}


HalfEdgeHalfFaceIter& HalfEdgeHalfFaceIter::operator++() {


    size_t n_outgoing_halffaces = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()].size();

    ++cur_index_;

    if (cur_index_ == n_outgoing_halffaces) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle((BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()])[cur_index_]);

    return *this;
}


} // namespace OpenVolumeMesh
