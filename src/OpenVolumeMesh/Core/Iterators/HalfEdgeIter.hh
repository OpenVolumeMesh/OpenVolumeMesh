#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT HalfEdgeIter : public BaseIterator<HalfEdgeHandle> {
public:
    typedef BaseIterator<HalfEdgeHandle> BaseIter;


    HalfEdgeIter(const TopologyKernel* _mesh, const HalfEdgeHandle& _heh = HalfEdgeHandle(0));

    // Post increment/decrement operator
    HalfEdgeIter operator++(int) {
        HalfEdgeIter cpy = *this;
        ++(*this);
        return cpy;
    }
    HalfEdgeIter operator--(int) {
        HalfEdgeIter cpy = *this;
        --(*this);
        return cpy;
    }
    HalfEdgeIter operator+(int _n) {
        HalfEdgeIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    HalfEdgeIter operator-(int _n) {
        HalfEdgeIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    HalfEdgeIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    HalfEdgeIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    HalfEdgeIter& operator++();
    HalfEdgeIter& operator--();

private:
    int cur_index_;
};


} // namespace OpenVolumeMesh
