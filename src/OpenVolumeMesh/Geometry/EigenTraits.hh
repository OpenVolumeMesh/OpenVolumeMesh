#pragma once
#include <OpenVolumeMesh/Geometry/Vector11T.hh>
#include <Eigen/Core>

namespace OpenVolumeMesh::Geometry
{

template<typename _Scalar, int _Cols>
struct vector_dim<Eigen::Matrix<_Scalar, 1, _Cols>> {
    constexpr static int dim = _Cols;
};
template<typename _Scalar, int _Rows>
struct vector_dim<Eigen::Matrix<_Scalar, _Rows, 1>> {
    constexpr static int dim = _Rows;
};
}

