#include <OpenVolumeMesh/Core/Iterators/VertexCellIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

VertexCellIter::VertexCellIter(const VertexHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_mesh->has_full_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
        BaseIter::valid(false);
        return;
    }

    // Build up cell list
    const std::vector<HalfEdgeHandle>& incidentHalfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h];
    for(std::vector<HalfEdgeHandle>::const_iterator it = incidentHalfedges.begin(); it != incidentHalfedges.end(); ++it) {

        if(!it->is_valid() || (unsigned int)it->idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) continue;
        const std::vector<HalfFaceHandle>& incidentHalfFaces = BaseIter::mesh()->incident_hfs_per_he_[*it];

        for(std::vector<HalfFaceHandle>::const_iterator hf_it = incidentHalfFaces.begin();
                hf_it != incidentHalfFaces.end(); ++hf_it) {
            if((unsigned int)hf_it->idx() < BaseIter::mesh()->incident_cell_per_hf_.size()) {
                CellHandle c_idx = BaseIter::mesh()->incident_cell_per_hf_[*hf_it];
                if(c_idx != TopologyKernel::InvalidCellHandle)
                    cells_.push_back(c_idx);
            }
        }
    }
    // Remove all duplicate entries
    std::sort(cells_.begin(), cells_.end());
    cells_.resize(std::unique(cells_.begin(), cells_.end()) - cells_.begin());

    // Remove invalid handles
    if ((cells_.size() > 0) && !cells_.front().is_valid())
        cells_.erase(cells_.begin());

    cur_index_ = 0;
    BaseIter::valid(cells_.size()>0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(cells_[cur_index_]);
    }
}

VertexCellIter& VertexCellIter::operator--() {

    if(cur_index_ == 0) {
        cur_index_ = cells_.size()-1;
        --lap_;
        if (lap_ < 0) {
            BaseIter::valid(false);
        }
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}


VertexCellIter& VertexCellIter::operator++() {

    ++cur_index_;
    if(cur_index_ == cells_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_) {
            BaseIter::valid(false);
        }
    }
    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}

} // namespace OpenVolumeMesh
