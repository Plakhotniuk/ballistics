//
// Created by Арсений Плахотнюк on 09.03.2023.
//

#ifndef BALLISTICS_BUTCHERTABLEERK4_HPP
#define BALLISTICS_BUTCHERTABLEERK4_HPP
#include "ballistics/types/BasicTypes.hpp"

namespace Ballistics::Integrators::ExplicitRK {

    // The classical Runge-Kutta method
    class ButcherTableERK4 {
    public:
        constexpr static index size = 4;
        constexpr static index numOfSubDiagEl = (size * size - size) / 2;
        constexpr static std::array<scalar, size> column = {
                static_cast<scalar>(0),
                static_cast<scalar>(1) / static_cast<scalar>(2),
                static_cast<scalar>(1) / static_cast<scalar>(2),
                static_cast<scalar>(1)
        };
        constexpr static std::array<scalar, size> row = {
                static_cast<scalar>(1) / static_cast<scalar>(6),
                static_cast<scalar>(2) / static_cast<scalar>(6),
                static_cast<scalar>(2) / static_cast<scalar>(6),
                static_cast<scalar>(1) / static_cast<scalar>(6)
        };
        constexpr static std::array<scalar, numOfSubDiagEl> matrix = {
                static_cast<scalar>(1) / static_cast<scalar>(2),
                static_cast<scalar>(0),
                static_cast<scalar>(1) / static_cast<scalar>(2),
                static_cast<scalar>(0),
                static_cast<scalar>(0),
                static_cast<scalar>(1)
        };
    };
}


#endif //BALLISTICS_BUTCHERTABLEERK4_HPP
