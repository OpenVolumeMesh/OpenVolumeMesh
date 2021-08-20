#include <OpenVolumeMesh/Core/Iterators/VertexIter.hh>
#include <OpenVolumeMesh/Core/TopologyKernel.hh>

namespace OpenVolumeMesh {

VertexIter::VertexIter(const TopologyKernel* _mesh, const VertexHandle& _vh) :
BaseIter(_mesh, _vh),
cur_index_(_vh.idx()) {

    while ((unsigned int)cur_index_ < BaseIter::mesh()->n_vertices() && BaseIter::mesh()->is_deleted(VertexHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->n_vertices()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(VertexHandle(cur_index_));
}


VertexIter& VertexIter::operator--() {

    --cur_index_;
    while (cur_index_ >= 0 && BaseIter::mesh()->is_deleted(VertexHandle(cur_index_)))
        --cur_index_;
    if(cur_index_ < 0) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(VertexHandle(cur_index_));
    return *this;
}


VertexIter& VertexIter::operator++() {

    ++cur_index_;
    while ((unsigned int)cur_index_ < BaseIter::mesh()->n_vertices() && BaseIter::mesh()->is_deleted(VertexHandle(cur_index_)))
        ++cur_index_;
    if((unsigned int)cur_index_ >= BaseIter::mesh()->n_vertices()) {
        BaseIter::valid(false);
    }
    BaseIter::cur_handle(VertexHandle(cur_index_));
    return *this;
}

} // namespace OpenVolumeMesh
