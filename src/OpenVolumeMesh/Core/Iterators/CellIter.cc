#include <OpenVolumeMesh/Core/Iterators/CellIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {


CellIter::CellIter(const TopologyKernel* _mesh, const CellHandle& _ch) :
BaseIter(_mesh, _ch),
cur_index_(_ch.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->cells_.size() && BaseIter::mesh()->is_deleted(CellHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->cells_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(CellHandle(cur_index_));
}


CellIter& CellIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(CellHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(CellHandle(cur_index_));
    return *this;
}


CellIter& CellIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->cells_.size() && BaseIter::mesh()->is_deleted(CellHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->cells_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(CellHandle(cur_index_));
    return *this;
}


} // namespace OpenVolumeMesh
