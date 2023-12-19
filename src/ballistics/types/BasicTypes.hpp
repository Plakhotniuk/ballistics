//
// Created by Арсений Плахотнюк on 15.09.2022.
//

#ifndef BALLISTICS_PLAKHOTNIUK_ARSENIY_BASICTYPES_H
#define BALLISTICS_PLAKHOTNIUK_ARSENIY_BASICTYPES_H
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace Ballistics {

    using scalar = double;
    using integer = std::int32_t;
    using index = uint;
    using Vector3d = Eigen::Vector<scalar, 3>;
    template<typename T, integer size> using Vector = Eigen::Vector<T, size>;

    using Matrix3x3 = Eigen::Matrix<scalar, 3, 3>;
    using Matrix6x6 = Eigen::Matrix<scalar, 6, 6>;

    template<typename T, index n, index m> using Matrix = Eigen::Matrix<T, n, m>;
    template<typename T, index n> using MatrixXX = Eigen::Matrix<T, n, n>;

    using Quaternion = Eigen::Quaternion<scalar>;
    using Quaterniond = Eigen::Quaterniond;


    template<typename T, index size>
    static bool isEqualByRelativeError(const Vector<T, size> &l, const Vector<T, size> &r, const scalar tolerance)
    {
        return (l - r).norm() <= tolerance * r.norm();
    }

    template<typename T, index n, index m>
    using diagMatrix = Eigen::DiagonalMatrix<T, n, m>;

    template<typename T, integer size>
    using vector = Eigen::Vector<T, size>;

    template<typename T, index n, index m>
    using matrix = Eigen::Matrix<T, n, m>;

    using quaterniond = Eigen::Quaterniond;

}

#endif //BALLISTICS_PLAKHOTNIUK_ARSENIY_BASICTYPES_H
