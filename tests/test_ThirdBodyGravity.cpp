#include "ballistics/forces/ThirdBodyGravity.hpp"
#include "../../third_party/include/rapidcsv.h"
#include "ballistics/types/BasicTypes.hpp"
#include "ballistics/exception/ProjectException.hpp"
#include "gtest/gtest.h"

bool isEqualByAbsoluteError(const Eigen::Vector3d &l,const Eigen::Vector3d &r, const Ballistics::scalar tolerance)
{
    return (l - r).norm() <= tolerance*r.norm();
}


TEST(THIRDBODYGRAVITY, TEST)
{
    using namespace Ballistics;
    const std::string FILE_PATH = __FILE__;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 25);
    std::string filename = dir_path + "utility/third_body_gravity.csv";
    rapidcsv::Document doc(filename);

    std::vector<scalar> xCol = doc.GetColumn<scalar>("x");
    std::vector<scalar> yCol = doc.GetColumn<scalar>("y");
    std::vector<scalar> zCol = doc.GetColumn<scalar>("z");
    std::vector<scalar> fxCol = doc.GetColumn<scalar>("fx");
    std::vector<scalar> fyCol = doc.GetColumn<scalar>("fy");
    std::vector<scalar> fzCol = doc.GetColumn<scalar>("fz");

    integer n = xCol.size();
    scalar tolerance = 1.e-5;

    const std::vector<Vector3d> SunAndMoonPositions = {{-9.66955e+10, -1.02917e+11, -4.46137e+10},
                                                       {-3.63564e+08, 6.23951e+06,3.86832e+07}};
    const std::vector<scalar> SunAndMoonParameters = {1.32712e+20, 4.9028e+12};
    const scalar mass = 1;

    for (int i = 0; i < n; ++i)
    {
        const Vector3d force = Ballistics::Forces::ThirdBodyGravity(SunAndMoonPositions,
                                                                  SunAndMoonParameters,
                                                                  {xCol[i], yCol[i], zCol[i]},
                                                                  mass);
//        std::cout<<force<<std::endl;
//        Vector3d ans = {fxCol[i], fyCol[i], fzCol[i]};
//        std::cout<<ans<<std::endl;
        ASSERT_TRUE(isEqualByAbsoluteError(force, {fxCol[i], fyCol[i], fzCol[i]}, tolerance));
    }
}