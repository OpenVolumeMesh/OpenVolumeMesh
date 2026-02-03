#pragma once
#include <OpenVolumeMesh/Geometry/Vector11T.hh>
#include <OpenVolumeMesh/Core/Properties/Defaults.hh>

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

namespace OpenVolumeMesh {
template<typename Scalar, int _Rows, int _Cols>
struct default_prop<Eigen::Matrix<Scalar, _Rows, _Cols>> {
    using MatrixT = Eigen::Matrix<Scalar, _Rows, _Cols>;
    static inline const MatrixT value = MatrixT::Zero();
};
}
