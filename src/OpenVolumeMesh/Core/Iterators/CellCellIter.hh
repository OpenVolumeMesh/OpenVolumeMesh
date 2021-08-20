#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT CellCellIter : public BaseCirculator<
    CellHandle,
    CellHandle> {
public:
    typedef BaseCirculator<
            CellHandle,
            CellHandle> BaseIter;

    CellCellIter(const CellHandle& _cIdx, const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    CellCellIter operator++(int) {
        CellCellIter cpy = *this;
        ++(*this);
        return cpy;
    }
    CellCellIter operator--(int) {
        CellCellIter cpy = *this;
        --(*this);
        return cpy;
    }
    CellCellIter operator+(int _n) {
        CellCellIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    CellCellIter operator-(int _n) {
        CellCellIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    CellCellIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    CellCellIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    CellCellIter& operator++();
    CellCellIter& operator--();

private:
    std::vector<CellHandle> adjacent_cells_;
    size_t cur_index_;
};

} // namespace OpenVolumeMesh
