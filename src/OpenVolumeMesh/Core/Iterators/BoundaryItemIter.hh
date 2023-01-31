#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>
#include <OpenVolumeMesh/Core/Iterators/VertexIter.hh>
#include <OpenVolumeMesh/Core/Iterators/EdgeIter.hh>
#include <OpenVolumeMesh/Core/Iterators/HalfEdgeIter.hh>
#include <OpenVolumeMesh/Core/Iterators/FaceIter.hh>
#include <OpenVolumeMesh/Core/Iterators/HalfFaceIter.hh>
#include <OpenVolumeMesh/Core/Iterators/CellIter.hh>


#include <iterator>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

template <class Iter, class Handle>
class OVM_EXPORT BoundaryItemIter : public BaseIterator<Handle> {
public:
    typedef BaseIterator<Handle> BaseIter;


    explicit BoundaryItemIter(const TopologyKernel* _mesh) :
    BaseIter(_mesh),
    it_(_mesh, Handle(0)),
    it_begin_(_mesh, Handle(0)),
    it_end_(_mesh, Handle((int)n_items())) {

        if(!has_incidences()) {
    #ifndef NDEBUG
            std::cerr << "This iterator needs bottom-up incidences!" << std::endl;
    #endif
            BaseIter::valid(false);
            return;
        }

        while(it_ != it_end_ && !BaseIter::mesh()->is_boundary(*it_)){
            ++it_;
        }
        BaseIter::valid(it_ != it_end_);
        if(BaseIter::valid()) {
            BaseIter::cur_handle(*it_);
        }
    }

    // Post increment/decrement operator
    BoundaryItemIter operator++(int) {
        BoundaryItemIter cpy = *this;
        ++(*this);
        return cpy;
    }
    BoundaryItemIter operator--(int) {
        BoundaryItemIter cpy = *this;
        --(*this);
        return cpy;
    }
    BoundaryItemIter operator+(int _n) {
        BoundaryItemIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    BoundaryItemIter operator-(int _n) {
        BoundaryItemIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    BoundaryItemIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    BoundaryItemIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    BoundaryItemIter& operator--() {
        --it_;
        while(it_ >= it_begin_ && !BaseIter::mesh()->is_boundary(*it_)){
            --it_;
        }
        if(it_ >= it_begin_) {
            BaseIter::cur_handle(*it_);
        } else {
            BaseIter::valid(false);
        }
        return *this;
    }

    BoundaryItemIter& operator++() {
        ++it_;
        while(it_ != it_end_ && !BaseIter::mesh()->is_boundary(*it_)){
            ++it_;
        }
        if(it_ != it_end_) {
            BaseIter::cur_handle(*it_);
        } else {
            BaseIter::valid(false);
        }
        return *this;
    }

private:
    size_t n_items() const;
    bool has_incidences() const;

private:
    Iter it_;
    const Iter it_begin_;
    const Iter it_end_;
};

//===========================================================================
typedef BoundaryItemIter<VertexIter, VertexHandle> BoundaryVertexIter;
typedef BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle> BoundaryHalfEdgeIter;
typedef BoundaryItemIter<EdgeIter, EdgeHandle> BoundaryEdgeIter;
typedef BoundaryItemIter<HalfFaceIter, HalfFaceHandle> BoundaryHalfFaceIter;
typedef BoundaryItemIter<FaceIter, FaceHandle> BoundaryFaceIter;
typedef BoundaryItemIter<CellIter, CellHandle> BoundaryCellIter;

// declare specializations so we can declare explicit instantiations:
//
template<> bool BoundaryItemIter<VertexIter, VertexHandle>::has_incidences() const;
template<> bool BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>::has_incidences() const;
template<> bool BoundaryItemIter<EdgeIter, EdgeHandle>::has_incidences() const;
template<> bool BoundaryItemIter<HalfFaceIter, HalfFaceHandle>::has_incidences() const;
template<> bool BoundaryItemIter<FaceIter, FaceHandle>::has_incidences() const;
template<> bool BoundaryItemIter<CellIter, CellHandle>::has_incidences() const;

// explicit instantiations:
extern template class OVM_EXPORT BoundaryItemIter<VertexIter, VertexHandle>;
extern template class OVM_EXPORT BoundaryItemIter<HalfEdgeIter, HalfEdgeHandle>;
extern template class OVM_EXPORT BoundaryItemIter<EdgeIter, EdgeHandle>;
extern template class OVM_EXPORT BoundaryItemIter<HalfFaceIter, HalfFaceHandle>;
extern template class OVM_EXPORT BoundaryItemIter<FaceIter, FaceHandle>;
extern template class OVM_EXPORT BoundaryItemIter<CellIter, CellHandle>;


} // namespace OpenVolumeMesh
