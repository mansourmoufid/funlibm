/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

/* Compile with -DNDEBUG, except for debugging. */

#pragma once

#include "constants.h" // pi, etc.
#include "polynomial.h" // eval_polynomial

static inline float _cos_pi_2_float(float x);

static inline double _cos_pi_2_double(double x);

#define _cos_0_pi_2(x) \
    _Generic((x), \
        float: _cos_0_pi_2_float, \
        double: _cos_0_pi_2_double \
    )(x)

static inline float
__attribute__((always_inline))
__attribute__((const))
_cos_0_pi_2_float(float x)
{
    assert(x >= 0.0f && x <= pi_2(float));
    const float as[22 + 1] = {
        /* x^0 */ 1.0f,
        /* x^1 */ 6.759776657698502e-13f,
        /* x^2 */ -0.5f,
        /* x^3 */ -1.9783012727980775e-10f,
        /* x^4 */ 0.0416666679084301f,
        /* x^5 */ 1.084509992921312e-08f,
        /* x^6 */ -0.0013890363043174148f,
        /* x^7 */ 4.3343567313058884e-07f,
        /* x^8 */ 2.706554369069636e-05f,
        /* x^9 */ -2.5627790819271468e-05f,
        /* x^10 */ 0.00011755149898817763f,
        /* x^11 */ -0.00034594471799209714f,
        /* x^12 */ 0.0007237704703584313f,
        /* x^13 */ -0.0011300748446956277f,
        /* x^14 */ 0.0013468097895383835f,
        /* x^15 */ -0.0012371520278975368f,
        /* x^16 */ 0.000876404985319823f,
        /* x^17 */ -0.00047500297660008073f,
        /* x^18 */ 0.00019345934560988098f,
        /* x^19 */ -5.732126737711951e-05f,
        /* x^20 */ 1.1666309546853881e-05f,
        /* x^21 */ -1.4586045153919258e-06f,
        /* x^22 */ 8.446725274779965e-08f,
    };
    return eval_polynomial(as, sizeof as / sizeof as[0], x);
}

static inline double
__attribute__((always_inline))
__attribute__((const))
_cos_0_pi_2_double(double x)
{
    assert(x >= 0.0 && x <= pi_2(double));
    const double as[22 + 1] = {
        /* x^0 */ 1.0,
        /* x^1 */ -1.1399247575174607e-21,
        /* x^2 */ -0.5,
        /* x^3 */ 3.7724774820191193e-19,
        /* x^4 */ 0.041666666666666664,
        /* x^5 */ -2.3875436436542907e-17,
        /* x^6 */ -0.0013888888888885747,
        /* x^7 */ -9.678259018392463e-16,
        /* x^8 */ 2.4801587297486778e-05,
        /* x^9 */ 4.924141452473386e-14,
        /* x^10 */ -2.755734187130068e-07,
        /* x^11 */ 6.595923551445847e-13,
        /* x^12 */ 2.0863109038375228e-09,
        /* x^13 */ 2.1045225666463804e-12,
        /* x^14 */ -1.3945869728049475e-11,
        /* x^15 */ 2.242711520404441e-12,
        /* x^16 */ -1.5190267372248051e-12,
        /* x^17 */ 8.374108806832156e-13,
        /* x^18 */ -3.364846344502883e-13,
        /* x^19 */ 9.827771930895603e-14,
        /* x^20 */ -1.9728105301385762e-14,
        /* x^21 */ 2.4332891569700873e-15,
        /* x^22 */ -1.3903679550728462e-16,
    };
    return eval_polynomial(as, sizeof as / sizeof as[0], x);
}
