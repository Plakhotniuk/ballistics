//
// Created by Арсений Плахотнюк on 09.03.2023.
//

#ifndef BALLISTICS_BUTCHERTABLEIRK6_HPP
#define BALLISTICS_BUTCHERTABLEIRK6_HPP

#include "ballistics/types/BasicTypes.hpp"

namespace Ballistics::Integrators::ImplicitRK {
    // The Gauss method
    class ButcherTableIRK6 {
    private:
        constexpr static scalar sqrt15 = static_cast<scalar>(3.872983346207416885179);
    public:
        constexpr static index size = 3;
        constexpr static std::array<scalar, size> column = {
                static_cast<scalar>(1) / static_cast<scalar>(2) - sqrt15 / static_cast<scalar>(10),
                static_cast<scalar>(1) / static_cast<scalar>(2),
                static_cast<scalar>(1) / static_cast<scalar>(2) + sqrt15 / static_cast<scalar>(10)
        };
        constexpr static std::array<scalar, size> row = {
                static_cast<scalar>(5) / static_cast<scalar>(18),
                static_cast<scalar>(4) / static_cast<scalar>(9),
                static_cast<scalar>(5) / static_cast<scalar>(18)
        };
        constexpr static inline std::array<std::array<scalar, size>, size> matrix { { {
                                                                                              static_cast<scalar>(5) / static_cast<scalar>(36),
                                                                                              static_cast<scalar>(2) / static_cast<scalar>(9) - sqrt15  / static_cast<scalar>(15),
                                                                                              static_cast<scalar>(5) / static_cast<scalar>(36) - sqrt15 / static_cast<scalar>(30)
                                                                                      }, {
                                                                                              static_cast<scalar>(5) / static_cast<scalar>(36) + sqrt15 / static_cast<scalar>(24),
                                                                                              static_cast<scalar>(2) / static_cast<scalar>(9),
                                                                                              static_cast<scalar>(5) / static_cast<scalar>(36) - sqrt15 / static_cast<scalar>(24)
                                                                                      }, {
                                                                                              static_cast<scalar>(5) / static_cast<scalar>(36) + sqrt15 / static_cast<scalar>(30),
                                                                                              static_cast<scalar>(2) / static_cast<scalar>(9) + sqrt15  / static_cast<scalar>(15),
                                                                                              static_cast<scalar>(5) / static_cast<scalar>(36)
                                                                                      } } };
    };
}



#endif //BALLISTICS_BUTCHERTABLEIRK6_HPP
