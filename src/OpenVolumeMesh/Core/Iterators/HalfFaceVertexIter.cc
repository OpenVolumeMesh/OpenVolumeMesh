#include <OpenVolumeMesh/Core/Iterators/HalfFaceVertexIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {


HalfFaceVertexIter::HalfFaceVertexIter(const HalfFaceHandle& _ref_h,
        const TopologyKernel* _mesh, int _max_laps) :
BaseIter(_mesh, _ref_h, _max_laps) {

    if(!_ref_h.is_valid()) return;

    cur_index_ = 0;

    BaseIter::valid(f_hehs().size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(cur_vh());
    }
}


HalfFaceVertexIter& HalfFaceVertexIter::operator--()
{
    if (cur_index_ == 0) {
        cur_index_  = f_hehs().size() - 1;
        --lap_;
        if (lap_ < 0) {
            BaseIter::valid(false);
        }
    }
    else {
        --cur_index_;
    }

    BaseIter::cur_handle(cur_vh());
    return *this;
}


HalfFaceVertexIter& HalfFaceVertexIter::operator++() {

    ++cur_index_;
    if (cur_index_ == f_hehs().size())
    {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_) {
            BaseIter::valid(false);
        }
    }
    BaseIter::cur_handle(cur_vh());
    return *this;
}

const std::vector<HalfEdgeHandle>& HalfFaceVertexIter::f_hehs() const
{
    return mesh()->face(ref_handle_.face_handle()).halfedges();
}

VertexHandle HalfFaceVertexIter::cur_vh() const
{
    const auto& hehs = f_hehs();
    if (ref_handle_.subidx()==1)
    {
        // Reversed Face
        return mesh()->to_vertex_handle(hehs[hehs.size()-1-cur_index_]);
    }
    else
    {
        // Original Face
        return mesh()->from_vertex_handle(hehs[cur_index_]);
    }
}

} // namespace OpenVolumeMesh
