#ifndef PROJECT_BALLISTICS_ATMOSPHERECONSTANTS_HPP
#define PROJECT_BALLISTICS_ATMOSPHERECONSTANTS_HPP

#include "../types/BasicTypes.hpp"
#include <array>

namespace Ballistics::Atmosphere
{
    // Константы и табличные параметры к модели атмосферы Земли по ГОСТ Р25645.166-2004
    // Ссылка на документацию:
    // https://docs.cntd.ru/document/1200036026

    constexpr scalar MIN_DISTANCE = 120000;  // нижняя граница высотного диапазона модели, м
    constexpr scalar MAX_DISTANCE = 1500000;  // верхняя граница высотного диапазона модели, м

    // Фиксированные уровни солнечной активности
    const std::array<int, 7> F0Levels = {75, 100, 125, 150, 175, 200, 250};

    // Индексы фиксированных уровней солнечной активности
    constexpr integer inline getF0Index(int F0) {
        if (F0 == 250) return 6;
        return (F0 - 75) / 25;
    }

    // Угловая скорость вращения Земли, рад/с
    constexpr scalar EarthOmega = 7.292115e-5;
    // Плотность ночной атмосферы на высоте 120 км, кг/м3
    constexpr scalar InitialNightDensity = 1.58868E-8;

    // Индексы массивов констант модели
    enum ConstsIndex { AI = 0, BI = 1, CI = 2, DI = 3, EI = 4, LI = 5 };
    // Размеры массивов констант модели (не физическая размерность)
    enum ConstsSizes {
        AISize = 7,
        ADISize = 9,
        BISize = 5,
        CISize = 5,
        DISize = 5,
        EISize = 9,
        EIKpSize = 4,
        LISize = 5,
        NISize = 3
    };

    constexpr scalar F81Weight = 60.75;         // посчитано ручками заранее
    constexpr scalar timeObsF10_7 = 20. / 24.;  // часть суток, на которое наблюдается F10.7.
    constexpr scalar timeObsKpDaily = 12. / 24.;  // часть суток, на которое наблюдается среднесуточный Kp.
    constexpr scalar timeObsKpThreeHours = 1.5 / 24.;  // часть суток, на которое наблюдается первый трехчасовой kp.

    // Коэффициент перевода спектральной плотности потока солнечного электромагнитного излучения из sfu (solar flux unit) в
    // СИ (10^-22 Вт/м^2Гц)
    constexpr scalar fromSFUToSI = 1e-22;

    // Время (в днях) запаздывания влияния индексов солнечной и геомагнитной активности на распределение плотности
    // в соответствующий им момент времени
    constexpr scalar F107TimeLag = 1.7;
    constexpr scalar ApTimeLag = 0.6;     // при использовании среднесуточных индексов
    constexpr scalar Ap3hTimeLag = 0.25;  // при использовании 3-часовых индексов

    constexpr scalar rPosDelta = 0.3;
    constexpr scalar rNegDelta = 0.7;

    // Средний Ap за последние 4 года (отсчитывая от 18.07.2021), используемый для прогноза геомагнитной активности
        constexpr scalar PropagatedAp = 7.747945205479452;
    // Todo: после того как напишешь модуль, посчитай более современное значение

    // Коэффициенты аппроксимации зависимости полугодового эффекта от времени года:
    constexpr std::array<scalar, 9> Ai = std::array<scalar, 9>{
    -2.53418e-2, -2.44075e-3, 3.08389e-6, 2.90115e-6, -4.99606e-8, 3.36327e-10, -1.0966e-12, 1.73227e-15, -1.06271e-18};

    // Данные для следующих ниже коэффициентов расположены так:
    // сначала половина коэффициентов для первого высотного диапазона, оставшаяся половина - для второго высотного диапазона
    // i строка отвечает i индексу коэффициента
    // j число в i строке отвечает значению коэффициента с индексом i при F0Levels[j]

    // Граница применимости между нижним и верхним диапазоном для различных коэффициентов.
    // Порядок определяется enum ConstsIndex
    constexpr std::array<std::array<scalar, 6>, 7> ConstsAltitudeThresholds = {
            std::array<scalar, 6>{500000, 600000, 640000, 1500000, 600000, 640000},
            std::array<scalar, 6>{500000, 660000, 700000, 1500000, 700000, 660000},
            std::array<scalar, 6>{500000, 760000, 760000, 1500000, 780000, 740000},
            std::array<scalar, 6>{500000, 800000, 820000, 1500000, 800000, 800000},
            std::array<scalar, 6>{500000, 860000, 860000, 1500000, 800000, 860000},
            std::array<scalar, 6>{500000, 900000, 920000, 1500000, 900000, 900000},
            std::array<scalar, 6>{500000, 1000000, 980000, 1500000, 760000, 900000}};

    constexpr std::array<std::array<scalar, 7>, 14> ai = {
            std::array<scalar, 7>{26.8629, 27.4598, 28.6395, 29.6418, 30.1671, 29.7578, 30.7854},
            std::array<scalar, 7>{-0.451674, -0.463668, -0.490987, -0.514957, -0.527837, -0.517915, -0.545695},
            std::array<scalar, 7>{0.00290397, 0.002974, 0.00320649, 0.00341926, 0.00353211, 0.00342699, 0.00370328},
            std::array<scalar, 7>{-1.06953e-5, -1.0753e-5, -1.1681e-5, -1.25785e-5, -1.30227e-5, -1.24137e-5, -1.37072e-5},
            std::array<scalar, 7>{2.21598e-8, 2.17059e-8, 2.36847e-8, 2.5727e-8, 2.66455e-8, 2.48209e-8, 2.80614e-8},
            std::array<scalar, 7>{-2.42941e-11, -2.30249e-11, -2.51809e-11, -2.75874e-11, -2.85432e-11, -2.58413e-11,
                                  -3.00184e-11},
            std::array<scalar, 7>{1.09926e-14, 1.00123e-14, 1.09536e-14, 1.21091e-14, 1.25009e-14, 1.09383e-14, 1.31142e-14},

            std::array<scalar, 7>{17.8781, -2.54909, -13.9599, -23.3079, -14.7264, -4.912, -5.40952},
            std::array<scalar, 7>{-0.132025, 0.0140064, 0.0844951, 0.135141, 0.0713256, 0.0108326, 0.00550749},
            std::array<scalar, 7>{0.000227717, -0.00016946, -0.000328875, -0.000420802, -0.000228015, -8.10546e-5, -3.78851e-5},
            std::array<scalar, 7>{-2.2543e-7, 3.27196e-7, 5.05918e-7, 5.73717e-7, 2.8487e-7, 1.15712e-7, 2.4808e-8},
            std::array<scalar, 7>{1.33574e-10, -2.8763e-10, -3.92299e-10, -4.03238e-10, -1.74383e-10, -8.13296e-11,
                                  4.92183e-12},
            std::array<scalar, 7>{-4.50458e-14, 1.22625e-13, 1.52279e-13, 1.42846e-13, 5.08071e-14, 3.04913e-14, -8.65011e-15},
            std::array<scalar, 7>{6.72086e-18, -2.05736e-17, -2.35576e-17, -2.01726e-17, -5.34955e-18, -4.94989e-18,
                                  1.9849e-18}};

    constexpr std::array<std::array<scalar, 7>, 10> bi = {
            std::array<scalar, 7>{0.0687894, 0.15073, 0.0479451, 0.0223448, -0.00326391, -0.0514749, -0.107255},
            std::array<scalar, 7>{-0.00284077, -0.00400889, -0.00239453, -0.0019798, -0.00159869, -0.000921059, -0.000174343},
            std::array<scalar, 7>{1.83922e-5, 2.43937e-5, 1.70335e-5, 1.54101e-5, 1.40443e-5, 1.15147e-5, 9.02759e-6},
            std::array<scalar, 7>{9.19605e-9, -9.92772e-9, -1.31626e-9, -2.3543e-9, -3.02287e-9, -1.22901e-9, -3.16512e-10},
            std::array<scalar, 7>{-4.16873e-11, -1.82239e-11, -1.74032e-11, -1.24994e-11, -9.2016e-12, -8.13104e-12, -6.14e-12},
            std::array<scalar, 7>{23.1584, 33.2732, 39.1961, 43.2469, 49.5738, 11.278, -52.6184},
            std::array<scalar, 7>{-0.0802147, -0.111099, -0.12352, -0.126973, -0.138613, 0.00143478, 0.214689},
            std::array<scalar, 7>{0.000105824, 0.000141421, 0.000149015, 0.000142637, 0.000147851, -3.69846e-5, -0.000294882},
            std::array<scalar, 7>{-6.15036e-8, -7.94952e-8, -7.9705e-8, -7.09985e-8, -6.96361e-8, 3.58318e-8, 1.71171e-7},
            std::array<scalar, 7>{1.32453e-11, 1.65836e-11, 1.58772e-11, 1.31646e-11, 1.21595e-11, -9.91225e-12, -3.60582e-11}};

    constexpr std::array<std::array<scalar, 7>, 10> ci = {
            std::array<scalar, 7>{-1.04825, -0.93106, -0.820867, -0.744047, -0.722471, -0.687482, -0.739984},
            std::array<scalar, 7>{0.0166305, 0.0141537, 0.0119916, 0.0104743, 0.00980317, 0.00916594, 0.00952854},
            std::array<scalar, 7>{-9.24263e-5, -7.29862e-5, -5.79835e-5, -4.78544e-5, -4.25245e-5, -3.80932e-5, -3.62727e-5},
            std::array<scalar, 7>{2.72382e-7, 2.00294e-7, 1.50707e-7, 1.18513e-7, 9.95544e-8, 8.51275e-8, 7.3887e-8},
            std::array<scalar, 7>{-2.41355e-10, -1.62006e-10, -1.13026e-10, -8.31498e-11, -6.55175e-11, -5.29972e-11,
                                  -4.23907e-11},
            std::array<scalar, 7>{50.5034, 61.624, 53.2623, 18.2236, -31.8442, -48.7208,
                                  -147.829},  // ci[5][6] исправлено с -147.859 на -147.829
            std::array<scalar, 7>{-0.170541, -0.192967, -0.144342, -0.00840024, 0.168327, 0.222996, 0.531652},
            std::array<scalar, 7>{2.17232e-4, 2.28061e-4, 1.4659e-4, -3.88e-5, -2.62603e-4, -3.21884e-4, -6.71937e-4},
            std::array<scalar, 7>{-1.21902e-7, -1.18715e-7, -6.46443e-8, 4.313846e-8, 1.65454e-7, 1.91495e-7, 3.64787e-7},
            std::array<scalar, 7>{2.54037e-11, 2.296380e-11, 1.04227e-11, -1.23832e-11, -3.69355e-11, -4.08067e-11,
                                  -7.26268e-11}};

    // Эти коэффициенты не зависят от высоты, но для удобства добавляем ещё раз
    constexpr std::array<std::array<scalar, 7>, 10> di = {
            std::array<scalar, 7>{-0.351899, -0.047813, 0.20981, 0.265174, 0.23047, 0.170074, 0.088141},
            std::array<scalar, 7>{0.00577056, 0.00380813, 0.00262881, 0.00275836, 0.00338331, 0.00406131, 0.00468253},
            std::array<scalar, 7>{9.95819e-7, 4.22771e-6, 4.24379e-6, 2.08668e-6, -5.52305e-7, -2.81614e-6,
                                  -4.24609e-6},  // di[3][5] исправлено с -2.82114e-6 на -2.81614e-6
            std::array<scalar, 7>{-7.25324e-9, -8.66826e-9, -6.67328e-9, -3.69543e-9, -8.23607e-10, 1.38369e-9, 2.53509e-9},
            std::array<scalar, 7>{2.9759e-12, 3.06712e-12, 2.13496e-12, 1.11862e-12, 2.21349e-13, -4.27908e-13, -7.29031e-13},

            std::array<scalar, 7>{-0.351899, -0.047813, 0.20981, 0.265174, 0.23047, 0.170074, 0.088141},
            std::array<scalar, 7>{0.00577056, 0.00380813, 0.00262881, 0.00275836, 0.00338331, 0.00406131, 0.00468253},
            std::array<scalar, 7>{9.95819e-7, 4.22771e-6, 4.24379e-6, 2.08668e-6, -5.52305e-7, -2.81614e-6,
                                  -4.24609e-6},  // di[8][5] исправлено с -2.82114e-6 на -2.81614e-6
            std::array<scalar, 7>{-7.25324e-9, -8.66826e-9, -6.67328e-9, -3.69543e-9, -8.23607e-10, 1.38369e-9, 2.53509e-9},
            std::array<scalar, 7>{2.9759e-12, 3.06712e-12, 2.13496e-12, 1.11862e-12, 2.21349e-13, -4.27908e-13, -7.29031e-13}};

    constexpr std::array<std::array<scalar, 7>, 18> ei = {
            std::array<scalar, 7>{-0.731596, -0.752175, -0.570476, -0.949573, -0.967598, -1.02278, -0.757903},
            std::array<scalar, 7>{0.00597345, 0.00565925, 2.95802e-3, 8.13121e-3, 8.41991e-3, 9.23633e-3, 0.00606068},
            std::array<scalar, 7>{-5.82037e-6, 1.8082e-6, 1.68896e-5, -3.87813e-6, -3.585e-6, -6.10128e-6, 7.85296e-6},
            std::array<scalar, 7>{6.84634e-8, 3.33822e-8, -4.7475e-9, 2.37694e-8, 1.74801e-8, 1.78211e-8, -9.74891e-9},
            std::array<scalar, 7>{-9.50483e-11, -5.13965e-11, -1.72711e-11, -2.77469e-11, -1.96221e-11, -1.70073e-11,
                                  1.58377e-12},
            std::array<scalar, 7>{-0.20670, -0.16971, -0.14671, -0.13150, -0.120916, -0.11363, -0.10444},
            std::array<scalar, 7>{9.7533e-2, 7.9830e-2, 6.8808e-2, 6.1603e-2, 5.6538e-2, 5.3178e-2, 4.8551e-2},
            std::array<scalar, 7>{-1.1817e-2, -9.4393e-3, -7.9836e-3, -7.0866e-3, -6.4324e-3, -6.0436e-3, -5.3567e-3},
            std::array<scalar, 7>{1.6145e-3, 1.2622e-3, 1.0535e-3, 9.2813e-4, 8.3723e-4, 7.7982e-4, 6.8809e-4},

            std::array<scalar, 7>{38.6199, 51.249, 68.4746, 58.422, 7.20188, 21.5948, -88.4076},
            std::array<scalar, 7>{-0.132147, -0.167373, -2.15659e-1, -1.66664e-1, 2.16109e-2, -2.02239e-2, 0.338518},
            std::array<scalar, 7>{1.75411e-4, 2.11832e-4, 2.62273e-4, 1.85486e-4, -6.52882e-5, -1.72029e-5, -0.000445581},
            std::array<scalar, 7>{-1.02417e-7, -1.18221e-7, -1.40972e-7, -9.12345e-8, 5.37077e-8, 2.83017e-8, 2.51729e-7},
            std::array<scalar, 7>{2.21446e-11, 2.45055e-11, 2.82285e-11, 1.67118e-11, -1.4095e-11, -8.94486e-12, -5.203e-11},
            std::array<scalar, 7>{-0.20670, -0.16971, -0.14671, -0.13150, -0.120916, -0.11363, -0.10444},
            std::array<scalar, 7>{9.7533e-2, 7.9830e-2, 6.8808e-2, 6.1603e-2, 5.6538e-2, 5.3178e-2, 4.8551e-2},
            std::array<scalar, 7>{-1.1817e-2, -9.4393e-3, -7.9836e-3, -7.0866e-3, -6.4324e-3, -6.0436e-3, -5.3567e-3},
            std::array<scalar, 7>{1.6145e-3, 1.2622e-3, 1.0535e-3, 9.2813e-4, 8.3723e-4, 7.7982e-4, 6.8809e-4}};

    // Эти коэффициенты не зависят от высоты
    constexpr std::array<std::array<scalar, 7>, 4> eti = {
            std::array<scalar, 7>{-0.2061, -0.169279, -0.146377, -0.13121, -0.12067, -0.113399, -0.104243},
            std::array<scalar, 7>{9.4449e-2, 7.7599e-2, 6.7052e-2, 6.0105e-2, 5.5232e-2, 5.1994e-2, 4.7573e-2},
            std::array<scalar, 7>{-8.7953e-3, -7.1375e-3, -6.0951e-3, -5.4388e-3, -4.9580e-3, -4.6876e-3, -4.1711e-3},
            std::array<scalar, 7>{8.8385e-4, 6.9025e-4, 5.7456e-4, 5.0585e-4, 4.5512e-4, 4.2548e-4, 3.7068e-4},
    };

    constexpr std::array<std::array<scalar, 7>, 10> li = {
            std::array<scalar, 7>{-0.407768, -0.902739, -0.733037, -1.31444, -1.20026, -1.52158, -1.67664},
            std::array<scalar, 7>{0.00148506, 0.00826803, 0.00523396, 0.0133124, 0.0114087, 0.015704, 1.77194e-2},
            std::array<scalar, 7>{1.25357e-5, -1.25448e-5, 6.35667e-6, -2.55585e-5, -1.47324e-5, -3.02859e-5, -3.69498e-5},
            std::array<scalar, 7>{3.77311e-8, 6.12853e-8, 1.09065e-8, 5.43981e-8, 2.7804e-8, 4.57668e-8, 5.09134e-8},
            std::array<scalar, 7>{-7.78953e-11, -7.07966e-11, -2.61427e-11, -4.33784e-11, -2.2632e-11, -2.82926e-11,
                                  -2.82878e-11},
            std::array<scalar, 7>{48.6536, 54.4867, 60.1267, 47.0996, 50.6174, 8.01942, -15.5728},
            std::array<scalar, 7>{-0.170291, -0.178298, -0.183144, -0.12526, -0.129047, 0.0185302, 9.36704e-2},
            std::array<scalar, 7>{2.26242e-4, 2.22725e-4, 2.12481e-4, 1.26352e-4, 1.24842e-4, -6.14733e-5, -1.49036e-4},
            std::array<scalar, 7>{-1.32032e-7, -1.227e-7, -1.08497e-7, -5.51584e-8, -5.24993e-8, 4.97674e-8, 9.42151e-8},
            std::array<scalar, 7>{2.85193e-11, 2.51316e-11, 2.0571e-11, 8.75272e-12, 8.08272e-12, -1.26162e-11, -2.0961e-11}};

    // Эти коэффициенты не зависят от уровня солнечной активности или высоты
    constexpr std::array<scalar, 3> ni = {2.058, 5.887e-3, -4.012e-6};

    // Этот коэффициент не зависит от высоты
    constexpr std::array<scalar, 7> phi1 = {0.5411, 0.5515, 0.5585, 0.5585, 0.5585, 0.5585, 0.5585};

    /** Данные из таблицы А.1 для перевода индекса Ap в индекс Kp **/
    // размерность = 2 нанотесла
    const std::array<scalar, 28> GeomagneticConversionArray = {
            0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32, 39, 48, 56, 67, 80, 94, 111, 132, 154, 179, 207, 236, 300, 400};
    // Измеряется в баллах
    const std::array<scalar, 28> Kps = {0.,       1. / 3.,  2. / 3.,  1.,       4. / 3.,  5. / 3.,  2.,
                                        7. / 3.,  8. / 3.,  3.,       10. / 3., 11. / 3., 4.,       13. / 3.,
                                        14. / 3., 5.,       16. / 3., 17. / 3., 6.,       19. / 3., 20. / 3.,
                                        7.,       22. / 3., 23. / 3., 8.,       25. / 3., 26. / 3., 9.};

    /** Константы для расчёта плотности при высоте < 120 км*/

    // Границы применимости к-тов, м
    constexpr std::array<scalar, 5> ConstsLowAltitudeThresholds = {0, 20, 60, 100, 120};

    // Сами коэффициенты

    // Размерность = кг/м^3
    constexpr std::array<scalar, 4> a0i = {1.228, 9.013e-2, 3.104e-4, 3.66e-7};
    // Размерность = 1/км
    constexpr std::array<scalar, 4> k1i = {-9.0764e-2, -0.16739, -0.137, -0.18553};
    // Размерность = 1/км^2
    constexpr std::array<scalar, 4> k2i = {-2.0452e-3, 6.2669e-4, -7.8653e-4, 1.5397e-3};
}

#endif //PROJECT_BALLISTICS_ATMOSPHERECONSTANTS_HPP
