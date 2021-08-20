#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT HalfEdgeCellIter : public BaseCirculator<
    HalfEdgeHandle,
    CellHandle> {
public:
    typedef BaseCirculator<
            HalfEdgeHandle,
            CellHandle> BaseIter;


    HalfEdgeCellIter(const HalfEdgeHandle& _heIdx, const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    HalfEdgeCellIter operator++(int) {
        HalfEdgeCellIter cpy = *this;
        ++(*this);
        return cpy;
    }
    HalfEdgeCellIter operator--(int) {
        HalfEdgeCellIter cpy = *this;
        --(*this);
        return cpy;
    }
    HalfEdgeCellIter operator+(int _n) {
        HalfEdgeCellIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    HalfEdgeCellIter operator-(int _n) {
        HalfEdgeCellIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    HalfEdgeCellIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    HalfEdgeCellIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    HalfEdgeCellIter& operator++();
    HalfEdgeCellIter& operator--();

private:
    CellHandle getCellHandle(int _cur_index) const;

private:
    std::vector<CellHandle> cells_;
    size_t cur_index_;
};

} // namespace OpenVolumeMesh
