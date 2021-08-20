#include <OpenVolumeMesh/Core/Iterators/HalfFaceIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {


HalfFaceIter::HalfFaceIter(const TopologyKernel* _mesh, const HalfFaceHandle& _hfh) :
BaseIter(_mesh, _hfh),
cur_index_(_hfh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->faces_.size() * 2 && BaseIter::mesh()->is_deleted(HalfFaceHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size() * 2) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfFaceHandle(cur_index_));
}


HalfFaceIter& HalfFaceIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(HalfFaceHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfFaceHandle(cur_index_));
    return *this;
}


HalfFaceIter& HalfFaceIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->faces_.size() * 2 && BaseIter::mesh()->is_deleted(HalfFaceHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size() * 2) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(HalfFaceHandle(cur_index_));
    return *this;
}


} // namespace OpenVolumeMesh
