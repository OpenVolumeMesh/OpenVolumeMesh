#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {
class OVM_EXPORT HalfFaceVertexIter : public BaseCirculator<
    HalfFaceHandle,
    VertexHandle> {
public:
    typedef BaseCirculator<
            HalfFaceHandle,
            VertexHandle> BaseIter;

    HalfFaceVertexIter(const HalfFaceHandle& _hIdx, const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    HalfFaceVertexIter operator++(int) {
        HalfFaceVertexIter cpy = *this;
        ++(*this);
        return cpy;
    }
    HalfFaceVertexIter operator--(int) {
        HalfFaceVertexIter cpy = *this;
        --(*this);
        return cpy;
    }
    HalfFaceVertexIter operator+(int _n) {
        HalfFaceVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    HalfFaceVertexIter operator-(int _n) {
        HalfFaceVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    HalfFaceVertexIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    HalfFaceVertexIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    HalfFaceVertexIter& operator++();
    HalfFaceVertexIter& operator--();

private:
    size_t cur_index_;

    const std::vector<HalfEdgeHandle>& f_hehs() const;

    VertexHandle cur_vh() const;
};

} // namespace OpenVolumeMesh
