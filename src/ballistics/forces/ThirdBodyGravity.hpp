//
// Created by Арсений Плахотнюк on 17.11.2022.
//

#ifndef BALLISTICS_THIRDBODYGRAVITY_HPP
#define BALLISTICS_THIRDBODYGRAVITY_HPP

#include "ballistics/types/BasicTypes.hpp"
#include <utility>
#include <vector>
namespace Ballistics::Forces {

    struct  ThirdBodyGravityForce{

        [[nodiscard]] static Vector3d calcForce(const scalar mass, const Vector3d& positionGCRS,
                                                const std::vector<Vector3d>& otherBodyPosGCRS,
                                                const std::vector<scalar>& gravParam) {
            /**
             * Функция рассчитывает гравитационную силу для отличных от Земли небесных тел
             * @param positionGCRS - вектор положения космического аппарата в GCRS
             * @param mass - масса космического аппарата
             * @param otherBodyPosGCRS - вектор векторов положения других небесных тел в GCRS
             * @param gravParam - вектор гравитационных параметров для других тел
             * @return вектор силы
             */

            Vector3d resultForce = Vector3d::Zero();
            for (int i = 0; i < otherBodyPosGCRS.size(); ++i) {
                const Vector3d ri_r = otherBodyPosGCRS[i] - positionGCRS;
                const scalar sqr1 = (ri_r).squaredNorm();
                const scalar norm1 = std::sqrt(sqr1);
                const scalar sqr2 = otherBodyPosGCRS[i].squaredNorm();
                const scalar norm2 = std::sqrt(sqr2);
                resultForce += gravParam[i] * ((ri_r) / (sqr1 * norm1) - otherBodyPosGCRS[i] / (sqr2 * norm2));
            }

            return resultForce * mass;
        };

    };

}

#endif //BALLISTICS_THIRDBODYGRAVITY_HPP
