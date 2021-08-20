#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT VertexCellIter : public BaseCirculator<
    VertexHandle,
    CellHandle> {
public:
    typedef BaseCirculator<
            VertexHandle,
            CellHandle> BaseIter;

    VertexCellIter(const VertexHandle& _vIdx, const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    VertexCellIter operator++(int) {
        VertexCellIter cpy = *this;
        ++(*this);
        return cpy;
    }
    VertexCellIter operator--(int) {
        VertexCellIter cpy = *this;
        --(*this);
        return cpy;
    }
    VertexCellIter operator+(int _n) {
        VertexCellIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    VertexCellIter operator-(int _n) {
        VertexCellIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    VertexCellIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    VertexCellIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    VertexCellIter& operator++();
    VertexCellIter& operator--();

private:
    std::vector<CellHandle> cells_;
    size_t cur_index_;
};


} // namespace OpenVolumeMesh
