#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT HalfFaceIter : public BaseIterator<HalfFaceHandle> {
public:
    typedef BaseIterator<HalfFaceHandle> BaseIter;


    HalfFaceIter(const TopologyKernel* _mesh, const HalfFaceHandle& _hfh = HalfFaceHandle(0));

    // Post increment/decrement operator
    HalfFaceIter operator++(int) {
        HalfFaceIter cpy = *this;
        ++(*this);
        return cpy;
    }
    HalfFaceIter operator--(int) {
        HalfFaceIter cpy = *this;
        --(*this);
        return cpy;
    }
    HalfFaceIter operator+(int _n) {
        HalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    HalfFaceIter operator-(int _n) {
        HalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    HalfFaceIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    HalfFaceIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    HalfFaceIter& operator++();
    HalfFaceIter& operator--();

private:
    int cur_index_;
};


} // namespace OpenVolumeMesh
