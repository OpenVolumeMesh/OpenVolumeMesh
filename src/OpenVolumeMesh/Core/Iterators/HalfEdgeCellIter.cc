#include <OpenVolumeMesh/Core/Iterators/HalfEdgeCellIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

#include <set>

namespace OpenVolumeMesh {

HalfEdgeCellIter::HalfEdgeCellIter(const HalfEdgeHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

    if(!_mesh->has_edge_bottom_up_incidences() || !_mesh->has_face_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

    if((unsigned int)_ref_h.idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) {

        BaseIter::valid(false);
        return;
    }
    if((unsigned int)cur_index_ >= BaseIter::mesh()->incident_hfs_per_he_[_ref_h].size()) {

        BaseIter::valid(false);
        return;
    }
    if((unsigned int)((BaseIter::mesh()->incident_hfs_per_he_[_ref_h])[cur_index_]).idx() >=
            BaseIter::mesh()->incident_cell_per_hf_.size()) {

        BaseIter::valid(false);
        return;
    }

    // collect cell handles
    const std::vector<HalfFaceHandle>& incidentHalffaces = BaseIter::mesh()->incident_hfs_per_he_[_ref_h];
    std::set<CellHandle> cells;
    for (unsigned int i = 0; i < incidentHalffaces.size(); ++i)
    {
        CellHandle ch = getCellHandle(i);
        if (ch.is_valid()) {
            if(cells.count(ch) == 0) {
                cells_.push_back(ch);
            }
            cells.insert(ch);
        }
    }

    BaseIter::valid(cells_.size() > 0);

    if(BaseIter::valid()) {
        BaseIter::cur_handle(cells_[cur_index_]);
    }
}


HalfEdgeCellIter& HalfEdgeCellIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = cells_.size()-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}


HalfEdgeCellIter& HalfEdgeCellIter::operator++() {

    ++cur_index_;

    if (cur_index_ == cells_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(cells_[cur_index_]);
    return *this;
}

CellHandle HalfEdgeCellIter::getCellHandle(int _cur_index) const
{
    const std::vector<HalfFaceHandle>& halffacehandles = BaseIter::mesh()->incident_hfs_per_he_[BaseIter::ref_handle()];
    HalfFaceHandle currentHalfface = halffacehandles[_cur_index];
    if(!currentHalfface.is_valid()) return CellHandle(-1);
    CellHandle cellhandle = BaseIter::mesh()->incident_cell_per_hf_[currentHalfface];
    return cellhandle;
}



} // namespace OpenVolumeMesh
