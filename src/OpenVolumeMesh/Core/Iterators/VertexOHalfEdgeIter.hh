#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {
class TopologyKernel;

class OVM_EXPORT VertexOHalfEdgeIter :
    public BaseCirculator<
  VertexHandle,
  HalfEdgeHandle> {
public:
    typedef BaseCirculator<
      VertexHandle,
      HalfEdgeHandle> BaseIter;


  VertexOHalfEdgeIter(const VertexHandle& _vIdx,
            const TopologyKernel* _mesh, int _max_laps = 1);

  // Post increment/decrement operator
  VertexOHalfEdgeIter operator++(int) {
    VertexOHalfEdgeIter cpy = *this;
    ++(*this);
    return cpy;
  }
  VertexOHalfEdgeIter operator--(int) {
    VertexOHalfEdgeIter cpy = *this;
    --(*this);
    return cpy;
  }
  VertexOHalfEdgeIter operator+(int _n) {
    VertexOHalfEdgeIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      ++cpy;
    }
    return cpy;
  }
  VertexOHalfEdgeIter operator-(int _n) {
    VertexOHalfEdgeIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      --cpy;
    }
    return cpy;
  }
  VertexOHalfEdgeIter& operator+=(int _n) {
    for(int i = 0; i < _n; ++i) {
      ++(*this);
    }
    return *this;
  }
  VertexOHalfEdgeIter& operator-=(int _n) {
    for(int i = 0; i < _n; ++i) {
      --(*this);
    }
    return *this;
  }

  VertexOHalfEdgeIter& operator++();
  VertexOHalfEdgeIter& operator--();

private:

    size_t cur_index_;
};

} // namespace OpenVolumeMesh
