#include <OpenVolumeMesh/Core/Iterators/CellVertexIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

CellVertexIter::CellVertexIter(const CellHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    for(const auto hfh: BaseIter::mesh()->cell_halffaces(_ref_h)) {
        for (const auto vh: BaseIter::mesh()->face_vertices(BaseIter::mesh()->face_handle(hfh))) {
            incident_vertices_.push_back(vh);
        }
    }

    // Remove all duplicate entries
    std::sort(incident_vertices_.begin(), incident_vertices_.end());
    incident_vertices_.erase(std::unique(incident_vertices_.begin(), incident_vertices_.end()), incident_vertices_.end());

    cur_index_ = 0;
    BaseIter::valid(incident_vertices_.size() > 0);

    if(BaseIter::valid()) {
        BaseIter::cur_handle(incident_vertices_[cur_index_]);
    }
}


CellVertexIter& CellVertexIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_ = incident_vertices_.size()-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(incident_vertices_[cur_index_]);
    return *this;
}


CellVertexIter& CellVertexIter::operator++() {

    ++cur_index_;
    if (cur_index_ == incident_vertices_.size()){
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }

    BaseIter::cur_handle(incident_vertices_[cur_index_]);
    return *this;
}

} // namespace OpenVolumeMesh
