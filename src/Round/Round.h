#pragma once

#include <cmath>
#include <iostream>

inline double RoundWithPrecision(double numberToRound, double precision) {
    return std::round((numberToRound) / precision) * precision;
}

inline double CeilWithPrecision(double numberToRound, double precision) {
    return std::ceil((numberToRound) / precision) * precision;
}

inline double FloorWithPrecision(double numberToRound, double precision) {
    return std::floor((numberToRound) / precision) * precision;
}

