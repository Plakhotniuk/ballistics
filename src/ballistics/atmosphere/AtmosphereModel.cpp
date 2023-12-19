#include "AtmosphereModel.hpp"

namespace Ballistics::Atmosphere
{
    scalar AtmosphereModel::calcDensity(const Eigen::Vector3d &r,
                                        const Time::Time<Time::Scale::UTC> &t,
                                        const scalar h,
                                        const scalar F10_7,
                                        const scalar F81,
                                        const integer F0,
                                        const scalar kpp,
                                        const scalar alpha,
                                        const scalar delta,
                                        const scalar S,
                                        const scalar d)
    {
        if (h < 0)
        {
            std::stringstream buff;
            buff << "The height above the Earth must be greater than zero. The received value: " << h << ". File: " << __FILE__ << ". Line: " << __LINE__;
            throw ProjectException(buff.str());
        }

        scalar density = -1.;

        if (0 <= h && h < 120)
        {
            integer i;
            for (i = 0; i < ConstsLowAltitudeThresholds.size() - 1; ++i)
            {
                if (h >= ConstsLowAltitudeThresholds[i] && h < ConstsLowAltitudeThresholds[i + 1])
                {
                    break;
                }
            }
            density = a0i[i] * exp(k1i[i] * (h - ConstsLowAltitudeThresholds[i]) + k2i[i] *
                        (h - ConstsLowAltitudeThresholds[i]) * (h - ConstsLowAltitudeThresholds[i]));
        }

        if (120 <= h && h <= 1500)
        {
            std::array<integer, 6> altitudeIndexes{};
            const std::array<integer, 6> altitudeSizes = {ConstsSizes::AISize, ConstsSizes::BISize,
                                                          ConstsSizes::CISize, ConstsSizes::DISize,
                                                          ConstsSizes::EISize, ConstsSizes::LISize};
            const integer F0Index = getF0Index(F0);

            for (int i = 0; i < 6; ++i)
            {
                if(h < ConstsAltitudeThresholds[F0Index][i] / 1000)
                {
                    altitudeIndexes[i] = 0;
                } else
                {
                    altitudeIndexes[i] = altitudeSizes[i];
                }
            }

            const std::vector<scalar> coefficientsAI = {ai[altitudeIndexes[ConstsIndex::AI] + 0][F0Index],
                                                        ai[altitudeIndexes[ConstsIndex::AI] + 1][F0Index],
                                                        ai[altitudeIndexes[ConstsIndex::AI] + 2][F0Index],
                                                        ai[altitudeIndexes[ConstsIndex::AI] + 3][F0Index],
                                                        ai[altitudeIndexes[ConstsIndex::AI] + 4][F0Index],
                                                        ai[altitudeIndexes[ConstsIndex::AI] + 5][F0Index],
                                                        ai[altitudeIndexes[ConstsIndex::AI] + 6][F0Index]};

            const scalar nightDensity = InitialNightDensity * exp(calcPartialSumPowerSeries(coefficientsAI, h));

            const std::vector<scalar> coefficientsLI = {li[altitudeIndexes[ConstsIndex::LI] + 0][F0Index],
                                                        li[altitudeIndexes[ConstsIndex::LI] + 1][F0Index],
                                                        li[altitudeIndexes[ConstsIndex::LI] + 2][F0Index],
                                                        li[altitudeIndexes[ConstsIndex::LI] + 3][F0Index],
                                                        li[altitudeIndexes[ConstsIndex::LI] + 4][F0Index]};

            const scalar K0 = 1 + (calcPartialSumPowerSeries(coefficientsLI, h)) * (F81 - F0) / F0;

            const scalar timeInSeconds = (t - SECONDS_IN_DAY / 2.).get_jdPart() * SECONDS_IN_DAY;
            const scalar beta = alpha - S - EarthOmega * timeInSeconds + phi1[F0Index];

            const scalar cos_phi = (1 / r.norm()) * (r.z() * sin(delta) + cos(delta) *
                                    (r.x() * cos(beta) + r.y() * sin(beta)));

            const scalar cos_phi_2 = sqrt((1 + cos_phi) / 2);

            const std::vector<scalar> coefficientsCI = {ci[altitudeIndexes[ConstsIndex::CI] + 0][F0Index],
                                                        ci[altitudeIndexes[ConstsIndex::CI] + 1][F0Index],
                                                        ci[altitudeIndexes[ConstsIndex::CI] + 2][F0Index],
                                                        ci[altitudeIndexes[ConstsIndex::CI] + 3][F0Index],
                                                        ci[altitudeIndexes[ConstsIndex::CI] + 4][F0Index]};

            const scalar K1 = (calcPartialSumPowerSeries(coefficientsCI, h)) *
                              pow(cos_phi_2, ni[0] + ni[1] * h + ni[2] * h * h);

            const std::vector<scalar> coefficientsADI = {Ai[0], Ai[1], Ai[2], Ai[3], Ai[4],
                                                         Ai[5], Ai[6], Ai[7], Ai[8]};

            const std::vector<scalar> coefficientsDI = {di[altitudeIndexes[ConstsIndex::DI] + 0][F0Index],
                                                        di[altitudeIndexes[ConstsIndex::DI] + 1][F0Index],
                                                        di[altitudeIndexes[ConstsIndex::DI] + 2][F0Index],
                                                        di[altitudeIndexes[ConstsIndex::DI] + 3][F0Index],
                                                        di[altitudeIndexes[ConstsIndex::DI] + 4][F0Index]};

            const scalar K2 = calcPartialSumPowerSeries(coefficientsDI, h) *
                              calcPartialSumPowerSeries(coefficientsADI, static_cast<scalar>(d));

            const std::vector<scalar> coefficientsBI = {bi[altitudeIndexes[ConstsIndex::BI] + 0][F0Index],
                                                        bi[altitudeIndexes[ConstsIndex::BI] + 1][F0Index],
                                                        bi[altitudeIndexes[ConstsIndex::BI] + 2][F0Index],
                                                        bi[altitudeIndexes[ConstsIndex::BI] + 3][F0Index],
                                                        bi[altitudeIndexes[ConstsIndex::BI] + 4][F0Index]};

            const scalar K3 = calcPartialSumPowerSeries(coefficientsBI, h) * (F10_7 - F81) /
                              (F81 + abs(F10_7 - F81));

            const std::vector<scalar> coefficientsEI = {ei[altitudeIndexes[ConstsIndex::EI] + 0][F0Index],
                                                        ei[altitudeIndexes[ConstsIndex::EI] + 1][F0Index],
                                                        ei[altitudeIndexes[ConstsIndex::EI] + 2][F0Index],
                                                        ei[altitudeIndexes[ConstsIndex::EI] + 3][F0Index],
                                                        ei[altitudeIndexes[ConstsIndex::EI] + 4][F0Index]};

            const std::vector<scalar> coefficientsETI = {eti[0][F0Index],
                                                         eti[1][F0Index],
                                                         eti[2][F0Index],
                                                         eti[3][F0Index]};

            const scalar K4 = calcPartialSumPowerSeries(coefficientsEI, h) *
                              calcPartialSumPowerSeries(coefficientsETI, kpp);


            density = nightDensity * K0 * (1 + K1 + K2 + K3 + K4);
        }

        if (h > 1500)
        {
            density = 0.;
        }

        return density;
    }
}
