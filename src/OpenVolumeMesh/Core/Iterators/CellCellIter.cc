#include <OpenVolumeMesh/Core/Iterators/CellCellIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

CellCellIter::CellCellIter(const CellHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    std::vector<HalfFaceHandle>::const_iterator hf_iter = BaseIter::mesh()->cell(_ref_h).halffaces().begin();
    std::vector<HalfFaceHandle>::const_iterator hf_end  = BaseIter::mesh()->cell(_ref_h).halffaces().end();
    for(; hf_iter != hf_end; ++hf_iter) {

        HalfFaceHandle opp_hf = BaseIter::mesh()->opposite_halfface_handle(*hf_iter);
        CellHandle ch = BaseIter::mesh()->incident_cell_per_hf_[opp_hf];
        if(ch != TopologyKernel::InvalidCellHandle) {
            adjacent_cells_.push_back(ch);
        }
    }

    // Remove all duplicate entries
    std::sort(adjacent_cells_.begin(), adjacent_cells_.end());
    adjacent_cells_.resize(std::unique(adjacent_cells_.begin(), adjacent_cells_.end()) - adjacent_cells_.begin());

    cur_index_ = 0;
    BaseIter::valid(adjacent_cells_.size()>0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(adjacent_cells_[cur_index_]);
    }
}


CellCellIter& CellCellIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_  = adjacent_cells_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(adjacent_cells_[cur_index_]);
    return *this;
}


CellCellIter& CellCellIter::operator++() {

    ++cur_index_;
    if (cur_index_ == adjacent_cells_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(adjacent_cells_[cur_index_]);
    return *this;
}

} // namespace OpenVolumeMesh
