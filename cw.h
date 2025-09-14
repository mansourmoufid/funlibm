/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

/* Compile with -DNDEBUG, except for debugging. */

#pragma once

#include <assert.h>
#include <math.h> // fmaf, fma
#include <stdint.h> // int32_t, int64_t

#include "arithmetic.h" // add, mul
#include "constants.h" // pi, etc.
#include "types.h" // rem_result_float, rem_result_double

/*
 * Variant of the Cody-Waite argument reduction algorithm.
 *
 * W. J. Cody and W. Waite, Software manual for elementary functions.
 * Prentice Hall, 1980.
 */

rem_result_float cw_rem_pi_float(float x);
rem_result_float cw_rem_2pi_float(float x);
rem_result_float cw_rem_pi_2_float(float x);

rem_result_double cw_rem_pi_double(double x);
rem_result_double cw_rem_2pi_double(double x);
rem_result_double cw_rem_pi_2_double(double x);

typedef struct {
    float C;
    float C1;
    float C2;
    float C3;
    float R;
    float R1;
    float R2;
    float R3;
} cw_parameters_float;

typedef struct {
    double C;
    double C1;
    double C2;
    double C3;
    double R;
    double R1;
    double R2;
    double R3;
} cw_parameters_double;

static const cw_parameters_float cw_parameters_pi_float = {
    .C = pi(float),
    .C1 = pi_hi(float),
    .C2 = pi_lo1(float),
    .C3 = pi_lo2(float),
    .R = inv_pi(float),
    .R1 = inv_pi_hi(float),
    .R2 = inv_pi_lo1(float),
    .R3 = inv_pi_lo2(float),
};

static const cw_parameters_double cw_parameters_pi_double = {
    .C = pi(double),
    .C1 = pi_hi(double),
    .C2 = pi_lo1(double),
    .C3 = pi_lo2(double),
    .R = inv_pi(double),
    .R1 = inv_pi_hi(double),
    .R2 = inv_pi_lo1(double),
    .R3 = inv_pi_lo2(double),
};

static const cw_parameters_float cw_parameters_2pi_float = {
    .C = twopi(float),
    .C1 = twopi_hi(float),
    .C2 = twopi_lo1(float),
    .C3 = twopi_lo2(float),
    .R = inv_2pi(float),
    .R1 = inv_2pi_hi(float),
    .R2 = inv_2pi_lo1(float),
    .R3 = inv_2pi_lo2(float),
};

static const cw_parameters_double cw_parameters_2pi_double = {
    .C = twopi(double),
    .C1 = twopi_hi(double),
    .C2 = twopi_lo1(double),
    .C3 = twopi_lo2(double),
    .R = inv_2pi(double),
    .R1 = inv_2pi_hi(double),
    .R2 = inv_2pi_lo1(double),
    .R3 = inv_2pi_lo2(double),
};

static const cw_parameters_float cw_parameters_pi_2_float = {
    .C = pi_2(float),
    .C1 = pi_2_hi(float),
    .C2 = pi_2_lo1(float),
    .C3 = pi_2_lo2(float),
    .R = inv_pi_2(float),
    .R1 = inv_pi_2_hi(float),
    .R2 = inv_pi_2_lo1(float),
    .R3 = inv_pi_2_lo2(float),
};

static const cw_parameters_double cw_parameters_pi_2_double = {
    .C = pi_2(double),
    .C1 = pi_2_hi(double),
    .C2 = pi_2_lo1(double),
    .C3 = pi_2_lo2(double),
    .R = inv_pi_2(double),
    .R1 = inv_pi_2_hi(double),
    .R2 = inv_pi_2_lo1(double),
    .R3 = inv_pi_2_lo2(double),
};

static inline rem_result_float
__attribute__((always_inline))
__attribute__((const))
_cw_correction_float(float x, const cw_parameters_float parameters, rem_result_float result)
{
    const float C = parameters.C;
    const float C1 = parameters.C1;
    const float C2 = parameters.C2;
    int32_t z = result.z;
    float v1 = result.v1;
    float v2 = result.v2;
    float r = v1 + v2;
    if (x >= 0.0f) {
        if (r < 0.0f) {
            // r += C;
            v1 += C1;
            v2 += C2;
            z -= 1;
            r = v1 + v2;
        }
        assert(r >= 0.0f && r <= C);
    } else {
        if (r > 0.0f) {
            // r -= C;
            v1 -= C1;
            v2 -= C2;
            z += 1;
            r = v1 + v2;
        }
        assert(r >= -C && r <= 0.0f);
    }
    return (rem_result_float) {z, v1, v2};
}

static inline rem_result_double
__attribute__((always_inline))
__attribute__((const))
_cw_correction_double(double x, const cw_parameters_double parameters, rem_result_double result)
{
    const double C = parameters.C;
    const double C1 = parameters.C1;
    const double C2 = parameters.C2;
    int64_t z = result.z;
    double v1 = result.v1;
    double v2 = result.v2;
    double r = v1 + v2;
    if (x >= 0.0) {
        if (r < 0.0) {
            // r += C;
            v1 += C1;
            v2 += C2;
            z -= 1;
            r = v1 + v2;
        }
        assert(r >= 0.0 && r <= C);
    } else {
        if (r > 0.0) {
            // r -= C;
            v1 -= C1;
            v2 -= C2;
            z += 1;
            r = v1 + v2;
        }
        assert(r >= -C && r <= 0.0);
    }
    return (rem_result_double) {z, v1, v2};
}

#define _cw_correction(x, parameters, result) \
    _Generic((x), \
        float: _cw_correction_float, \
        double: _cw_correction_double \
    )(x, parameters, result)

static inline tuple_float
__attribute__((always_inline))
__attribute__((const))
cw_reduce_float(const cw_parameters_float parameters, float x, int32_t q)
{
    const float C1 = parameters.C1;
    const float C2 = parameters.C2;
    const float C3 = parameters.C3;
    float r, e;

    // r = x - (q * C1) - (q * C2);
    // r = fmaf(-q, C2, fmaf(-q, C1, x));
    // e = 0.0f;

    // r = x - (q * C1) - (q * C2) - (q * C3);
    // r = fmaf(-q, C3, fmaf(-q, C2, fmaf(-q, C1, x)));
    // e = 0.0f;

    float qf = (float) q;
    // r1 = x - q * C1
    tuple_float x1 = xfma(-qf, C1, x);
    float r1 = x1.a;
    float e1 = x1.b;
    // r2 = r1 - q * C2
    tuple_float x2 = xfma(-qf, C2, r1);
    float r2 = x2.a;
    float e2 = x2.b;
    // r3 = r2 - q * C3
    tuple_float x3 = xfma(-qf, C3, r2);
    float r3 = x3.a;
    float e3 = x3.b;
    r = r3;
    e = e1 + e2 + e3;

    return (tuple_float) {r, e};
}

static inline rem_result_float
__attribute__((always_inline))
__attribute__((const))
_cw_rem_float(const cw_parameters_float parameters, float x)
{
    const float R1 = parameters.R1;
    const float R2 = parameters.R2;
    const float R3 = parameters.R3;

    // const float C = parameters.C;
    // const float R = parameters.R;
    // float t = x * R;
    // int32_t q = (int32_t) t;
    // float r = x - q * C;
    // float r = fmaf(-q, C, x);

    // float t = (x * R1) + (x * R2);
    // float t = fmaf(x, R1, x * R2);

    // float t = (x * R1) + (x * R2) + (x * R3);
    float t = fmaf(x, R1, fmaf(x, R2, x * R3));

    int32_t q = (int32_t) t;

    tuple_float rem = cw_reduce_float(parameters, x, q);
    float r = rem.a;
    float e = rem.b;

    rem_result_float result = {q, r, e};
    // return result;
    return _cw_correction(x, parameters, result);
}

static inline rem_result_double
__attribute__((always_inline))
__attribute__((const))
_cw_rem_double(const cw_parameters_double parameters, double x)
{
    const double C1 = parameters.C1;
    const double C2 = parameters.C2;
    const double C3 = parameters.C3;
    const double R1 = parameters.R1;
    const double R2 = parameters.R2;
    const double R3 = parameters.R3;

    // const double C = parameters.C;
    // const double R = parameters.R;
    // double t = x * R;
    // int64_t q = (int64_t) t;
    // double r = x - q * C;
    // double r = fma(-q, C, x);

    // double t = (x * R1) + (x * R2);
    // double t = fma(x, R1, x * R2);

    // double t = (x * R1) + (x * R2) + (x * R3);
    double t = fma(x, R1, fma(x, R2, x * R3));

    int64_t q = (int64_t) t;

    // double r = x - (q * C1) - (q * C2);
    // double r = fma(-q, C2, fma(-q, C1, x));

    // double r = x - (q * C1) - (q * C2) - (q * C3);
    // double r = fma(-q, C3, fma(-q, C2, fma(-q, C1, x)));
    // double e = 0.0;

    double qd = (double) q;
    // r1 = x - q * C1
    tuple_double x1 = xfma(-qd, C1, x);
    double r1 = x1.a;
    double e1 = x1.b;
    // r2 = r1 - q * C2
    tuple_double x2 = xfma(-qd, C2, r1);
    double r2 = x2.a;
    double e2 = x2.b;
    // r3 = r2 - q * C3
    tuple_double x3 = xfma(-qd, C3, r2);
    double r3 = x3.a;
    double e3 = x3.b;
    double r = r3;
    double e = e1 + e2 + e3;

    rem_result_double result = {q, r, e};
    // return result;
    return _cw_correction(x, parameters, result);
}

#define _cw_rem(parameters, x) \
    _Generic((parameters), \
        cw_parameters_float: _cw_rem_float, \
        cw_parameters_double: _cw_rem_double \
    )(parameters, x)

static inline rem_result_float
__attribute__((always_inline))
__attribute__((const))
_cw_rem_pi_float(float x)
{
    return _cw_rem(cw_parameters_pi_float, x);
}

static inline rem_result_double
__attribute__((always_inline))
__attribute__((const))
_cw_rem_pi_double(double x)
{
    return _cw_rem(cw_parameters_pi_double, x);
}

#define _cw_rem_pi(x) \
    _Generic((x), \
        float: _cw_rem_pi_float, \
        double: _cw_rem_pi_double \
    )(x)

static inline rem_result_float
__attribute__((always_inline))
__attribute__((const))
_cw_rem_2pi_float(float x)
{
    return _cw_rem(cw_parameters_2pi_float, x);
}

static inline rem_result_double
__attribute__((always_inline))
__attribute__((const))
_cw_rem_2pi_double(double x)
{
    return _cw_rem(cw_parameters_2pi_double, x);
}

#define _cw_rem_2pi(x) \
    _Generic((x), \
        float: _cw_rem_2pi_float, \
        double: _cw_rem_2pi_double \
    )(x)

static inline rem_result_float
__attribute__((always_inline))
__attribute__((const))
_cw_rem_pi_2_float(float x)
{
    return _cw_rem(cw_parameters_pi_2_float, x);
}

static inline rem_result_double
__attribute__((always_inline))
__attribute__((const))
_cw_rem_pi_2_double(double x)
{
    return _cw_rem(cw_parameters_pi_2_double, x);
}

#define _cw_rem_pi_2(x) \
    _Generic((x), \
        float: _cw_rem_pi_2_float, \
        double: _cw_rem_pi_2_double \
    )(x)
