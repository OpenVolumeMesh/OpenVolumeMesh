#include <OpenVolumeMesh/Core/Iterators/HalfFaceVertexIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {


HalfFaceVertexIter::HalfFaceVertexIter(const HalfFaceHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_ref_h.is_valid()) return;

    const OpenVolumeMeshFace& halfface = _mesh->halfface(_ref_h);
    const std::vector<HalfEdgeHandle>& hes = halfface.halfedges();
    for(std::vector<HalfEdgeHandle>::const_iterator he_it = hes.begin();
            he_it != hes.end(); ++he_it) {
        vertices_.push_back(_mesh->halfedge(*he_it).from_vertex());
    }

    cur_index_ = 0;

    BaseIter::valid(vertices_.size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(vertices_[cur_index_]);
    }
}


HalfFaceVertexIter& HalfFaceVertexIter::operator--() {

    if (cur_index_ == 0) {
        cur_index_  = vertices_.size() - 1;
        --lap_;
        if (lap_ < 0)
            BaseIter::valid(false);
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(vertices_[cur_index_]);
    return *this;
}


HalfFaceVertexIter& HalfFaceVertexIter::operator++() {

    ++cur_index_;
    if (cur_index_ == vertices_.size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_)
            BaseIter::valid(false);
    }
    BaseIter::cur_handle(vertices_[cur_index_]);
    return *this;
}

} // namespace OpenVolumeMesh
