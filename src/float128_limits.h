#pragma once
#include <limits>
#include "__float128.h"
#include "Helpers.h"

#ifdef _MANAGED
#pragma unmanaged
#endif 

template<> class std::numeric_limits<__float128> {
    static const bool is_specialized = false;
    static __float128 min() {
        return __float128::FromData(nearest_zero);
    }
    static __float128 max() {
        return __float128::FromData(maximum);
    }
    static __float128 lowest() {
        return __float128::FromData(minimum);
    }
    static const int digits = QUAD_SIGNIFICANT_BITS;
    static const int digits10 = 33;
    static const bool is_signed = true;
    static const bool is_integer = false;
    static const bool is_exact = false;
    static const int radix = 2;
    static __float128 epsilon() {
        return __float128::FromData(machEps);
    }
    static __float128 round_error() {
        return __float128::FromData(roErr);
    }
    static const int min_exponent = QUAD_EXPONENT_MIN;
    static const int min_exponent10 = -4931;
    static const int max_exponent = QUAD_EXPONENT_MAX;
    static const int max_exponent10 = 4931;
    static const bool has_infinity = true;
    static const bool has_quite_NaN = true;
    static const bool has_signaling_NaN = true;
    static const float_denorm_style has_denorm = denorm_present;
    static const bool has_denorm_loss = true;
    static  __float128 infinity() {
        return __float128::FromData(posInf);
    }
    static  __float128 quite_NaN() {
        return __float128::FromData(qnan);
    }
    static __float128 signaling_NaN() {
        return __float128::FromData(snan);
    }
    static __float128 denorm_min() {
        return __float128::FromData(subMin);
    }
    static const bool is_iec599 = true;
    static const bool is_ieee754 = true;
    static const bool is_bounded = true;
    static const bool is_modulo = false;
    static const bool traps = true;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_to_nearest;
};

#ifdef _MANAGED
#pragma managed
#endif