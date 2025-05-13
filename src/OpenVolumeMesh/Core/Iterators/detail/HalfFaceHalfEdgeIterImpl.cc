#include <OpenVolumeMesh/Core/Iterators/detail/HalfFaceHalfEdgeIterImpl.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh::detail {

HalfFaceHalfEdgeIterImpl::HalfFaceHalfEdgeIterImpl(const HalfFaceHandle& _ref_h, const TopologyKernel* _mesh, int _max_laps) :
    BaseIter(_mesh, _ref_h, _max_laps),
    cur_index_(0)
{
    BaseIter::valid(f_hehs().size() > 0);
    if(BaseIter::valid()) {
        BaseIter::cur_handle(cur_heh());
    }
}

HalfFaceHalfEdgeIterImpl& HalfFaceHalfEdgeIterImpl::operator--()
{
    if (cur_index_ == 0) {
        cur_index_ = f_hehs().size() - 1;
        --lap_;
        if (lap_ < 0) {
            BaseIter::valid(false);
        }
    }
    else {
        --cur_index_;
    }
    BaseIter::cur_handle(cur_heh());
    return *this;
}

HalfFaceHalfEdgeIterImpl& HalfFaceHalfEdgeIterImpl::operator++()
{
    ++cur_index_;
    if (cur_index_ >= f_hehs().size()) {
        cur_index_ = 0;
        ++lap_;
        if (lap_ >= max_laps_) {
            BaseIter::valid(false);
        }
    }
    BaseIter::cur_handle(cur_heh());
    return *this;
}

const std::vector<HalfEdgeHandle>& HalfFaceHalfEdgeIterImpl::f_hehs() const
{
    return mesh()->face(ref_handle_.face_handle()).halfedges();
}

HalfEdgeHandle HalfFaceHalfEdgeIterImpl::cur_heh() const
{
    if (ref_handle_.subidx()==1)
    {
        // Reversed Face
        return f_hehs()[f_hehs().size()-1-cur_index_].opposite_handle();
    }
    else
    {
        // Original Face
        return f_hehs()[cur_index_];
    }
}


} // namespace OpenVolumeMesh::detail
