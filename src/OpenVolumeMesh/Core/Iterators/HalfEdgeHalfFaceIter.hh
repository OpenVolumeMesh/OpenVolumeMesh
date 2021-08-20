#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT HalfEdgeHalfFaceIter : public BaseCirculator<
    HalfEdgeHandle,
    HalfFaceHandle> {
public:
    typedef BaseCirculator<
            HalfEdgeHandle,
            HalfFaceHandle> BaseIter;


    HalfEdgeHalfFaceIter(const HalfEdgeHandle& _heIdx, const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    HalfEdgeHalfFaceIter operator++(int) {
        HalfEdgeHalfFaceIter cpy = *this;
        ++(*this);
        return cpy;
    }
    HalfEdgeHalfFaceIter operator--(int) {
        HalfEdgeHalfFaceIter cpy = *this;
        --(*this);
        return cpy;
    }
    HalfEdgeHalfFaceIter operator+(int _n) {
        HalfEdgeHalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    HalfEdgeHalfFaceIter operator-(int _n) {
        HalfEdgeHalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    HalfEdgeHalfFaceIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    HalfEdgeHalfFaceIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    HalfEdgeHalfFaceIter& operator++();
    HalfEdgeHalfFaceIter& operator--();

private:
    size_t cur_index_;
};

} // namespace OpenVolumeMesh
