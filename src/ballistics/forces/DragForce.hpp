#ifndef PROJECT_BALLISTICS_DRAGFORCE_HPP
#define PROJECT_BALLISTICS_DRAGFORCE_HPP

#include "../types/BasicTypes.hpp"

namespace Ballistics::Forces
{
    struct DragForce{

        [[nodiscard]] static Vector3d calcForce(const Vector3d &velocity, const scalar density,
                                                const scalar S, const scalar Cd = 2.2) noexcept{
         /** Метод для расчета силы трения КА об атмосферу
         * @param density - плотность атмосферы
         * @param velocity - скорость КА
         * @param S - площадь
         * @param Cd - коэффициент трения об атмосферу
         * @return вектор силы трения
         */
            return - (1. / 2.) * density * velocity.norm() * Cd * S *velocity;
        }
    };
}

#endif //PROJECT_BALLISTICS_DRAGFORCE_HPP