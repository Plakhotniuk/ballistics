//
// Created by Арсений Плахотнюк on 07.06.2023.
//

#ifndef BALLISTICS_DIFFERENTIALCORRECTION_HPP
#define BALLISTICS_DIFFERENTIALCORRECTION_HPP

#include "ballistics/types/BasicTypes.hpp"

namespace Ballistics::Optimization::DiffCorrection
{
    template<typename RightPart>
    struct DiffCorrectionParameters
    {
        const typename RightPart::State initState;
        const std::vector<diagMatrix<scalar, RightPart::size, RightPart::size>> W;
        const vector<scalar, RightPart::size> derivativeSteps;
        const scalar integrationStep;
        const scalar tolerance;
    };

    template<typename Integrator, typename RightPart>
    [[nodiscard]] std::pair<vector<scalar, RightPart::size>, matrix<scalar, RightPart::size, RightPart::size>>
    Optimize(const std::vector<typename RightPart::State> &measurements,
             const typename RightPart::Params &modelParams,
             const DiffCorrectionParameters<RightPart> &diffCorrParams)
    {
        vector<scalar, RightPart::size> updatedPosition = diffCorrParams.initState.u;
        matrix<scalar, RightPart::size, RightPart::size> blockA;
        matrix<scalar, RightPart::size, RightPart::size> covarianceMatrix;
        scalar RMS;
        scalar tmpRMS = 0.;

        do
        {
            RMS = tmpRMS;
            tmpRMS = 0.;

            matrix<scalar, RightPart::size, RightPart::size> firstPart =
                    matrix<scalar, RightPart::size, RightPart::size>::Zero();
            vector<scalar, RightPart::size> secondPart = vector<scalar, RightPart::size>::Zero();
            typename RightPart::State tmpState = {updatedPosition, diffCorrParams.initState.t};
            vector<typename RightPart::State, RightPart::size> tmpStatesVec;

            for (integer j = 0; j < RightPart::size; ++j)
            {
                vector<scalar, RightPart::size> increment = vector<scalar, RightPart::size>::Zero();
                increment[j] += diffCorrParams.derivativeSteps[j];
                tmpStatesVec[j] = {updatedPosition + increment, diffCorrParams.initState.t};
            }

            // Here we go through our measurements
            for (integer i = 0; i < measurements.size(); ++i)
            {
                const std::vector<typename RightPart::State> masOfUnperturbedPositions =
                        Integrator::template calc<RightPart>(tmpState,
                                                             modelParams,
                                                             {diffCorrParams.integrationStep, measurements[i].t});
                const vector<scalar, RightPart::size> unperturbedPosition = masOfUnperturbedPositions.back().u;
                tmpState = masOfUnperturbedPositions.back();

                // Calculating of the matrix A
                // Calculating of the i-th value of unperturbed solution
                for (integer j = 0; j < RightPart::size; ++j)
                {
                    const std::vector<typename RightPart::State> masOfPerturbedPositions =
                            Integrator::template calc<RightPart>(tmpStatesVec[j],
                                                                 modelParams,
                                                                 {diffCorrParams.integrationStep, measurements[i].t});
                    const vector<scalar, RightPart::size> perturbedPosition = masOfPerturbedPositions.back().u;
                    tmpStatesVec[j] = masOfPerturbedPositions.back();
                    blockA.col(j) = (perturbedPosition - unperturbedPosition) / diffCorrParams.derivativeSteps[j];
                }
                // End of calculating matrix A

                const vector<scalar, RightPart::size> deltaY = measurements[i].u - unperturbedPosition;
                firstPart += blockA.transpose() * diffCorrParams.W[i] * blockA;
                secondPart += blockA.transpose() * diffCorrParams.W[i] * deltaY;
                tmpRMS += sqrt((deltaY.transpose() * diffCorrParams.W[i]).dot(deltaY) / static_cast<double>(deltaY.size()));
            }

            covarianceMatrix = firstPart.inverse();
            updatedPosition += covarianceMatrix * secondPart;

        } while(std::abs(RMS - tmpRMS) > diffCorrParams.tolerance * tmpRMS);

        return {updatedPosition, covarianceMatrix};
    }
}



#endif //BALLISTICS_DIFFERENTIALCORRECTION_HPP
