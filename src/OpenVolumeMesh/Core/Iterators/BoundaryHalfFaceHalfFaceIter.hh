#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT BoundaryHalfFaceHalfFaceIter : public BaseCirculator<HalfFaceHandle,
    HalfFaceHandle> {
private:
    typedef BaseCirculator<HalfFaceHandle,
            HalfFaceHandle> BaseIter;
public:
    BoundaryHalfFaceHalfFaceIter(const HalfFaceHandle& _ref_h,
            const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    BoundaryHalfFaceHalfFaceIter operator++(int) {
        BoundaryHalfFaceHalfFaceIter cpy = *this;
        ++(*this);
        return cpy;
    }
    BoundaryHalfFaceHalfFaceIter operator--(int) {
        BoundaryHalfFaceHalfFaceIter cpy = *this;
        --(*this);
        return cpy;
    }
    BoundaryHalfFaceHalfFaceIter operator+(int _n) {
        BoundaryHalfFaceHalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    BoundaryHalfFaceHalfFaceIter operator-(int _n) {
        BoundaryHalfFaceHalfFaceIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    BoundaryHalfFaceHalfFaceIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    BoundaryHalfFaceHalfFaceIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    const EdgeHandle& common_edge() const { return common_edges_[cur_index_]; }

    BoundaryHalfFaceHalfFaceIter& operator++();
    BoundaryHalfFaceHalfFaceIter& operator--();

private:
    std::vector<HalfFaceHandle> neighbor_halffaces_;
    std::vector<EdgeHandle> common_edges_;
    size_t cur_index_;
};

} // namespace OpenVolumeMesh
