#include <OpenVolumeMesh/Core/Iterators/EdgeIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {


EdgeIter::EdgeIter(const TopologyKernel* _mesh, const EdgeHandle& _eh) :
BaseIter(_mesh, _eh),
cur_index_(_eh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->edges_.size() && BaseIter::mesh()->is_deleted(EdgeHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(EdgeHandle(cur_index_));
}


EdgeIter& EdgeIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(EdgeHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(EdgeHandle(cur_index_));
    return *this;
}


EdgeIter& EdgeIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->edges_.size() && BaseIter::mesh()->is_deleted(EdgeHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->edges_.size()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(EdgeHandle(cur_index_));
    return *this;
}


} // namespace OpenVolumeMesh
