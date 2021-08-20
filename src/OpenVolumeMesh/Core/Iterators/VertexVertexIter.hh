#pragma once

#include <OpenVolumeMesh/Core/Handles.hh>
#include <OpenVolumeMesh/Config/Export.hh>
#include <OpenVolumeMesh/Core/Iterators/BaseCirculator.hh>

#ifndef NDEBUG
#include <iostream>
#endif

namespace OpenVolumeMesh {

class OVM_EXPORT VertexVertexIter :
    public BaseCirculator<
  VertexHandle,
  VertexHandle> {
public:
    typedef BaseCirculator<
      VertexHandle,
      VertexHandle> BaseIter;


  VertexVertexIter(const VertexHandle& _vIdx,
            const TopologyKernel* _mesh, int _max_laps = 1);

  // Post increment/decrement operator
  VertexVertexIter operator++(int) {
    VertexVertexIter cpy = *this;
    ++(*this);
    return cpy;
  }
  VertexVertexIter operator--(int) {
    VertexVertexIter cpy = *this;
    --(*this);
    return cpy;
  }
  VertexVertexIter operator+(int _n) {
    VertexVertexIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      ++cpy;
    }
    return cpy;
  }
  VertexVertexIter operator-(int _n) {
    VertexVertexIter cpy = *this;
    for(int i = 0; i < _n; ++i) {
      --cpy;
    }
    return cpy;
  }
  VertexVertexIter& operator+=(int _n) {
    for(int i = 0; i < _n; ++i) {
      ++(*this);
    }
    return *this;
  }
  VertexVertexIter& operator-=(int _n) {
    for(int i = 0; i < _n; ++i) {
      --(*this);
    }
    return *this;
  }

  VertexVertexIter& operator++();
  VertexVertexIter& operator--();

private:

    size_t cur_index_;
};



} // namespace OpenVolumeMesh
