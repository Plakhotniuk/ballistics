//
// Created by Арсений Плахотнюк on 09.03.2023.
//

#ifndef BALLISTICS_EXPLICITRUNGEKUTTA_HPP
#define BALLISTICS_EXPLICITRUNGEKUTTA_HPP

#include "ballistics/types/BasicTypes.hpp"

namespace Ballistics::Integrators::ExplicitRK {

    template<typename RightPart, typename T>
    struct IntegrationParameters {
        const scalar step;
        const T tFinal;
    };

    template<typename ButcherTable>
    class ExplicitRungeKutta {
    public:
        template<typename RightPart, typename T>
        [[nodiscard]] static inline std::vector<typename RightPart::State> calc(const typename RightPart::State& initState,
                                                                                const typename RightPart::Params& params,
                                                                                const IntegrationParameters<RightPart, T>& intParams) {
            const index numOfIter = static_cast<index>(std::abs(initState.t - intParams.tFinal) / intParams.step);
            Vector<scalar, RightPart::size> uTmp = initState.u;
            typename RightPart::TimeState tTmp = initState.t;
            std::vector<typename RightPart::State> solution = {initState};

            for (integer j = 0; j < numOfIter; ++j) {
                // Calculating k[i]
                std::vector<Vector<scalar, RightPart::size>> kMas = {RightPart::calc({uTmp, tTmp}, params)};
                index n = 0;
                for (integer i = 1; i < ButcherTable::size; ++i) {
                    Vector<scalar, RightPart::size> sum = Vector<scalar, RightPart::size>::Zero();
                    for (integer k = 0; k < i; ++k)
                    {
                        sum += ButcherTable::matrix[n + k] * kMas[k];
                    }
                    n += i;

                    kMas.push_back(RightPart::calc({uTmp + sum * intParams.step, tTmp + ButcherTable::column[i] * intParams.step},
                                                   params));
                }
                // Calculating sum k[i] * b[i];
                Vector<scalar, RightPart::size> sum = Vector<scalar, RightPart::size>::Zero();
                for (integer i = 0; i < ButcherTable::size; ++i) {
                    sum += kMas[i] * ButcherTable::row[i];
                }

                uTmp = uTmp + intParams.step * sum;
                tTmp = tTmp + intParams.step;

                solution.push_back({uTmp, tTmp});
            }

            return solution;
        }
    };
}



#endif //BALLISTICS_EXPLICITRUNGEKUTTA_HPP
