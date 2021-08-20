#include <OpenVolumeMesh/Core/Iterators/VertexVertexIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

VertexVertexIter::VertexVertexIter(const VertexHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps),
cur_index_(0) {

  if(!_mesh->has_vertex_bottom_up_incidences()) {
#ifndef NDEBUG
        std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
#endif
        BaseIter::valid(false);
        return;
    }

  if((size_t)_ref_h.idx() >= BaseIter::mesh()->outgoing_hes_per_vertex_.size()) {
    BaseIter::valid(false);
  }

  if(BaseIter::valid()) {
    if((size_t)cur_index_ >= BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h].size()) {
      BaseIter::valid(false);
    }
  }

  if(BaseIter::valid()) {
    HalfEdgeHandle heh = BaseIter::mesh()->outgoing_hes_per_vertex_[_ref_h][cur_index_];
    BaseIter::cur_handle(BaseIter::mesh()->halfedge(heh).to_vertex());
  }
}


VertexVertexIter& VertexVertexIter::operator--() {

    size_t n_outgoing_halfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle()].size();

    if (cur_index_ == 0) {
        cur_index_ = n_outgoing_halfedges-1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    HalfEdgeHandle heh = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle()][cur_index_];
    BaseIter::cur_handle(BaseIter::mesh()->halfedge(heh).to_vertex());

  return *this;
}


VertexVertexIter& VertexVertexIter::operator++() {

    size_t n_outgoing_halfedges = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle()].size();

    ++cur_index_;

    if (cur_index_ == n_outgoing_halfedges) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }


    HalfEdgeHandle heh = BaseIter::mesh()->outgoing_hes_per_vertex_[BaseIter::ref_handle()][cur_index_];
    BaseIter::cur_handle(BaseIter::mesh()->halfedge(heh).to_vertex());

  return *this;
}



} // namespace OpenVolumeMesh
