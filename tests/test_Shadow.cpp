#include "ballistics/shadow/Shadow.hpp"
#include "../../third_party/include/rapidcsv.h"
#include "ballistics/types/BasicTypes.hpp"
#include "ballistics/exception/ProjectException.hpp"
#include "gtest/gtest.h"

TEST(SHADOW, TEST)
{
    using namespace Ballistics;
    const std::string FILE_PATH = __FILE__;
    const std::string dir_path = FILE_PATH.substr(0, FILE_PATH.size() - 15);
    std::string filename = dir_path + "utility/penumbra_shadow.csv";
    rapidcsv::Document doc(filename);
    std::vector<scalar> xCol = doc.GetColumn<scalar>("x");
    std::vector<scalar> yCol = doc.GetColumn<scalar>("y");
    std::vector<scalar> zCol = doc.GetColumn<scalar>("z");
    std::vector<scalar> shadowCol = doc.GetColumn<scalar>("shadow");

    integer n = xCol.size();
    scalar tolerance = 5e-14;

    for (int i = 0; i < n; ++i)
    {
        ASSERT_TRUE(abs(Shadow::calcShadowFunction({xCol[i], yCol[i], zCol[i]}, {-1.47098291e10, 0., 0.}, 695.700e6,  6.371e6)
            - shadowCol[i]) <= tolerance);
    }
}