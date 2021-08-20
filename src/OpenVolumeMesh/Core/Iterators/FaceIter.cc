#include <OpenVolumeMesh/Core/Iterators/FaceIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {


FaceIter::FaceIter(const TopologyKernel* _mesh, const FaceHandle& _fh) :
BaseIter(_mesh, _fh),
cur_index_(_fh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->faces_.size() && BaseIter::mesh()->is_deleted(FaceHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(FaceHandle(cur_index_));
}


FaceIter& FaceIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(FaceHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(FaceHandle(cur_index_));
    return *this;
}


FaceIter& FaceIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->faces_.size() && BaseIter::mesh()->is_deleted(FaceHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->faces_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(FaceHandle(cur_index_));
    return *this;
}


} // namespace OpenVolumeMesh
