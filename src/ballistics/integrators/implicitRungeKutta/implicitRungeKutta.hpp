//
// Created by Арсений Плахотнюк on 09.03.2023.
//

#ifndef BALLISTICS_IMPLICITRUNGEKUTTA_HPP
#define BALLISTICS_IMPLICITRUNGEKUTTA_HPP
#include "ballistics/types/BasicTypes.hpp"

namespace Ballistics::Integrators::ImplicitRK {

    template<typename RightPart>
    struct IntegrationParameters {
        const scalar step;
        const typename RightPart::Time tFinal;
        const scalar toleranceOfSNAE;
    };

    template<typename ButcherTable>
    class ImplicitRungeKutta {
    public:
        template<typename RightPart>
        [[nodiscard]] static inline std::vector<typename RightPart::State> calc(const typename RightPart::State &initState,
                                                                                const typename RightPart::Params &params,
                                                                                const IntegrationParameters<RightPart> &intParams) {

            const auto numOfIter = static_cast<index>(std::abs(initState.t - intParams.tFinal) / intParams.step);
            Vector<scalar, RightPart::size> uTmp = initState.u;
            typename RightPart::Time tTmp = initState.t;
            std::vector<typename RightPart::State> solution = {initState};

            Matrix<scalar, ButcherTable::size, RightPart::size> k =
                    Matrix<scalar, ButcherTable::size, RightPart::size>::Zero();
            for (integer l = 0; l < numOfIter; ++l) {

                // Simple iteration method for systems of nonlinear algebraic equations
                bool key = false;
                while (!key) {
                    Matrix<scalar, ButcherTable::size, RightPart::size> kTmp;
                    // Here we go through each row of the matrix k
                    for (integer i = 0; i < ButcherTable::size; ++i) {
                        // calculating sum of a[i][j] * k[j], where j in range from 0 to ButcherTable::size
                        Vector<scalar, RightPart::size> sum = Vector<scalar, RightPart::size>::Zero();
                        for (integer j = 0; j < ButcherTable::size; ++j) {
                            sum += ButcherTable::matrix[i][j] * k.row(j);
                        }
                        kTmp.row(i) = RightPart::calc({uTmp + intParams.step * sum,
                                                       tTmp + intParams.step * ButcherTable::column[i]}, params);
                    }
                    // Here we will check whether we have achieved the required accuracy in solving the system
                    key = true;
                    for (integer i = 0; i < ButcherTable::size; ++i) {
                        key *= isEqualByRelativeError<scalar,RightPart::size>(kTmp.row(i),  k.row(i), intParams.toleranceOfSNAE);
                    }
                    k = kTmp;
                }

                // Calculating sum of b[i] * k[i], where i in range from 0 to ButcherTable::size
                Vector<scalar, RightPart::size> sum = Vector<scalar, RightPart::size>::Zero();
                for (integer i = 0; i < ButcherTable::size; ++i) {
                    sum += ButcherTable::row[i] * k.row(i);
                }

                uTmp = uTmp + intParams.step * sum;
                tTmp = tTmp + intParams.step;
                solution.push_back({uTmp, tTmp});
            }

            return solution;
        }
    };
}



#endif //BALLISTICS_IMPLICITRUNGEKUTTA_HPP
