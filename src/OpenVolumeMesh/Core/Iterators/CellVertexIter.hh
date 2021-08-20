#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {
class TopologyKernel;

class OVM_EXPORT CellVertexIter : public BaseCirculator<
    CellHandle,
    VertexHandle> {
public:
    typedef BaseCirculator<
            CellHandle,
            VertexHandle> BaseIter;

    CellVertexIter(const CellHandle& _cIdx, const TopologyKernel* _mesh, int _max_laps = 1);

    // Post increment/decrement operator
    CellVertexIter operator++(int) {
        CellVertexIter cpy = *this;
        ++(*this);
        return cpy;
    }
    CellVertexIter operator--(int) {
        CellVertexIter cpy = *this;
        --(*this);
        return cpy;
    }
    CellVertexIter operator+(int _n) {
        CellVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            ++cpy;
        }
        return cpy;
    }
    CellVertexIter operator-(int _n) {
        CellVertexIter cpy = *this;
        for(int i = 0; i < _n; ++i) {
            --cpy;
        }
        return cpy;
    }
    CellVertexIter& operator+=(int _n) {
        for(int i = 0; i < _n; ++i) {
            ++(*this);
        }
        return *this;
    }
    CellVertexIter& operator-=(int _n) {
        for(int i = 0; i < _n; ++i) {
            --(*this);
        }
        return *this;
    }

    CellVertexIter& operator++();
    CellVertexIter& operator--();

private:
    std::vector<VertexHandle> incident_vertices_;
    size_t cur_index_;
};

} // namespace OpenVolumeMesh
