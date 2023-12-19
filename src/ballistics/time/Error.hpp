//
// Created by Арсений Плахотнюк on 01.10.2022.
//

#ifndef BALLISTICS_ERROR_HPP
#define BALLISTICS_ERROR_HPP
enum Error {
    UNACCEPTABLE_DATE = -1, //iauUt1utc, iauTaiutc, iauUtctai, iauUtcut1
    STATUS_OK = 0,
    DUBIOUS_YEAR = 1, //iauUt1utc, iauTaiutc, iauUtctai, iauUtcut1
    DUT_ERROR = 2,
    FILE_PARSER_ERROR = 3,
    INDEX_OUT_OF_RANGE = 4
};
#endif //BALLISTICS_ERROR_HPP
