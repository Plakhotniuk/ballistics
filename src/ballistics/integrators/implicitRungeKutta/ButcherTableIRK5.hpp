//
// Created by Арсений Плахотнюк on 09.03.2023.
//

#ifndef BALLISTICS_BUTCHERTABLEIRK5_HPP
#define BALLISTICS_BUTCHERTABLEIRK5_HPP

#include "ballistics/types/BasicTypes.hpp"

namespace Ballistics::Integrators::ImplicitRK {
    // The Radau method
    class ButcherTableIRK5 {
    private:
        constexpr static scalar sqrt6 = static_cast<scalar>(2.4494897427831780981972);
    public:
        constexpr static index size = 3;
        constexpr static std::array<scalar, size> column = {
                (static_cast<scalar>(4) - sqrt6) / static_cast<scalar>(10),
                (static_cast<scalar>(4) + sqrt6) / static_cast<scalar>(10),
                static_cast<scalar>(1)
        };
        constexpr static std::array<scalar, size> row = {
                (static_cast<scalar>(16) - sqrt6) / static_cast<scalar>(36),
                (static_cast<scalar>(16) + sqrt6) / static_cast<scalar>(36),
                static_cast<scalar>(1)            / static_cast<scalar>(9)
        };
        constexpr static inline std::array<std::array<scalar, size>, size> matrix { { {
                                                                                              (static_cast<scalar>(88)  - static_cast<scalar>(7) * sqrt6)   / static_cast<scalar>(360),
                                                                                              (static_cast<scalar>(296) - static_cast<scalar>(169) * sqrt6) / static_cast<scalar>(1800),
                                                                                              (static_cast<scalar>(-2)  + static_cast<scalar>(3) * sqrt6)   / static_cast<scalar>(255)
                                                                                      }, {
                                                                                              (static_cast<scalar>(296) + static_cast<scalar>(169) * sqrt6) / static_cast<scalar>(1800),
                                                                                              (static_cast<scalar>(88)  + static_cast<scalar>(7) * sqrt6)   / static_cast<scalar>(360),
                                                                                              (static_cast<scalar>(-2)  - static_cast<scalar>(3) * sqrt6)   / static_cast<scalar>(255)
                                                                                      }, {
                                                                                              (static_cast<scalar>(16) - sqrt6) / static_cast<scalar>(36),
                                                                                              (static_cast<scalar>(16) + sqrt6) / static_cast<scalar>(36),
                                                                                              static_cast<scalar>(1)            / static_cast<scalar>(9)
                                                                                      } } };
    };
}

#endif //BALLISTICS_BUTCHERTABLEIRK5_HPP
