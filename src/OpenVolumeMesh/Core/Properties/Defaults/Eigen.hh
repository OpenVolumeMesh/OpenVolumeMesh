#include <OpenVolumeMesh/Core/Properties/Defaults.hh>
#include <OpenVolumeMesh/Geometry/Vector11T.hh>
#include <Eigen/Core>

namespace OpenVolumeMesh {

template<typename Scalar, int _Rows, int _Cols>
struct default_prop<Eigen::Matrix<Scalar, _Rows, _Cols>
                   ,std::enable_if_t<_Rows >= 0 && _Cols >= 0 >>
                   {
    using MatrixT = Eigen::Matrix<Scalar, _Rows, _Cols>;
    static inline const MatrixT value = MatrixT::Zero();
};

template<typename Scalar, int _Rows>
struct default_prop<Eigen::Matrix<Scalar, _Rows, Eigen::Dynamic>
                   ,std::enable_if_t<_Rows >= 0>>
{
    using MatrixT = Eigen::Matrix<Scalar, _Rows, Eigen::Dynamic>;
    static inline const MatrixT value = MatrixT::Zero(_Rows, 1);
};

template<typename Scalar, int _Cols>
struct default_prop<Eigen::Matrix<Scalar, Eigen::Dynamic, _Cols>
                   ,std::enable_if_t<_Cols >= 0>>
{
    using MatrixT = Eigen::Matrix<Scalar, Eigen::Dynamic, _Cols>;
    static inline const MatrixT value = MatrixT::Zero(1, _Cols);
};
template<typename Scalar>
struct default_prop<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>> {
    using MatrixT = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
    static inline const MatrixT value = MatrixT::Zero(1, 1);
};

template<typename Scalar, int Options>
struct default_prop<Eigen::Quaternion<Scalar, Options>>
{
    using QuatT = Eigen::Quaternion<Scalar, Options>;
    static inline const QuatT value = QuatT(
            default_prop<Scalar>::value,
            default_prop<Scalar>::value,
            default_prop<Scalar>::value,
            default_prop<Scalar>::value);
};

} // namespace OpenVolumeMesh
