#include <OpenVolumeMesh/Core/Iterators/HalfEdgeIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {


HalfEdgeIter::HalfEdgeIter(const TopologyKernel* _mesh, const HalfEdgeHandle& _heh) :
BaseIter(_mesh, _heh),
cur_index_(_heh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->edges_.size() * 2 && BaseIter::mesh()->is_deleted(HalfEdgeHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size() * 2) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
}


HalfEdgeIter& HalfEdgeIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(HalfEdgeHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
    return *this;
}


HalfEdgeIter& HalfEdgeIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->edges_.size() * 2 && BaseIter::mesh()->is_deleted(HalfEdgeHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size() * 2) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfEdgeHandle(cur_index_));
    return *this;
}

} // namespace OpenVolumeMesh
