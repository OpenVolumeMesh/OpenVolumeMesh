#include <OpenVolumeMesh/Core/Iterators/VertexFaceIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

VertexFaceIter::VertexFaceIter(const VertexHandle& _ref_h,
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

    // Build up face list
    const std::vector<HalfEdgeHandle>& incidentHalfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h];
    for(std::vector<HalfEdgeHandle>::const_iterator it = incidentHalfedges.begin(); it != incidentHalfedges.end(); ++it) {

        if(!it->is_valid() || (unsigned int)it->idx() >= BaseIter::mesh()->incident_hfs_per_he_.size()) continue;
            const std::vector<HalfFaceHandle>& incidentHalfFaces = BaseIter::mesh()->incident_hfs_per_he_[*it];

        for (std::vector<HalfFaceHandle>::const_iterator hf_it = incidentHalfFaces.begin();
                hf_it != incidentHalfFaces.end(); ++hf_it) {
            faces_.push_back(BaseIter::mesh()->face_handle(*hf_it));
        }
    }
    // Remove all duplicate entries
    std::sort(faces_.begin(), faces_.end());
    faces_.resize(std::unique(faces_.begin(), faces_.end()) - faces_.begin());

    // Remove invalid handles
    if ((faces_.size() > 0) && !faces_.front().is_valid())
        faces_.erase(faces_.begin());

    cur_index_ = 0;
    BaseIter::valid(faces_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(faces_[cur_index_]);
    }
}

VertexFaceIter& VertexFaceIter::operator--() {

    if(cur_index_ == 0) {
        cur_index_ = faces_.size()-1;
        --lap_;
        if (lap_ < 0) {
            BaseIter::valid(false);
        }
    } else {
        --cur_index_;
    }

    BaseIter::cur_handle(faces_[cur_index_]);
    return *this;
}


VertexFaceIter& VertexFaceIter::operator++() {

    ++cur_index_;
    if(cur_index_ == faces_.size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_) {
            BaseIter::valid(false);
        }
  }
    BaseIter::cur_handle(faces_[cur_index_]);
  return *this;
}



} // namespace OpenVolumeMesh
