#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseIterator.hh>

#ifndef NDEBUG
#include <iostream>
#endif


namespace OpenVolumeMesh {
class OVM_EXPORT FaceIter : public BaseIterator<FaceHandle> {
public:
    typedef BaseIterator<FaceHandle> BaseIter;

    FaceIter(const TopologyKernel* _mesh, const FaceHandle& _fh = FaceHandle(0));

    // Post increment/decrement operator
    FaceIter operator++(int) {
        FaceIter cpy = *this;
        ++(*this);
        return cpy;
    }
    FaceIter operator--(int) {
        FaceIter cpy = *this;
        --(*this);
        return cpy;
    }
    FaceIter operator+(int _n) {
        FaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    FaceIter operator-(int _n) {
        FaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    FaceIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    FaceIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    FaceIter& operator++();
    FaceIter& operator--();

private:
    int cur_index_;
};


} // namespace OpenVolumeMesh
