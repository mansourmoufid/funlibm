/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

/* Compile with -DNDEBUG, except for debugging. */

#pragma once

#include <assert.h>
#include <math.h> // fabs, fma
#include <stdint.h> // int32_t, int64_t

#include "arithmetic.h" // add, mul
#include "constants.h" // pi, etc.
#include "types.h" // rem_result_double, tuple_double

/*
 * Naive argument reduction (for comparison).
 */

static inline rem_result_float
_naive_rem_pi_float(float x)
{
    int32_t q = (int32_t) (x / pi(float));
    return (rem_result_float) {q, x - (float) q * pi(float), 0.0f};
}

static inline rem_result_double
_naive_rem_pi_double(double x)
{
    int64_t q = (int64_t) (x / pi_double);
    return (rem_result_double) {q, x - (double) q * pi_double, 0.0};
}

#define _naive_rem_pi(x) \
    _Generic((x), \
        float: _naive_rem_pi_float, \
        double: _naive_rem_pi_double \
    )(x)

static inline rem_result_float
_naive_rem_2pi_float(float x)
{
    int32_t q = (int32_t) (x / twopi(float));
    return (rem_result_float) {q, x - (float) q * twopi(float), 0.0f};
}

static inline rem_result_double
_naive_rem_2pi_double(double x)
{
    int64_t q = (int64_t) (x / twopi_double);
    return (rem_result_double) {q, x - (double) q * twopi_double, 0.0};
}

#define _naive_rem_2pi(x) \
    _Generic((x), \
        float: _naive_rem_2pi_float, \
        double: _naive_rem_2pi_double \
    )(x)

static inline rem_result_float
_naive_rem_pi_2_float(float x)
{
    int32_t q = (int32_t) (x / pi_2(float));
    return (rem_result_float) {q, x - (float) q * pi_2(float), 0.0f};
}

static inline rem_result_double
_naive_rem_pi_2_double(double x)
{
    int64_t q = (int64_t) (x / pi_2_double);
    return (rem_result_double) {q, x - (double) q * pi_2_double, 0.0};
}

#define _naive_rem_pi_2(x) \
    _Generic((x), \
        float: _naive_rem_pi_2_float, \
        double: _naive_rem_pi_2_double \
    )(x)

/*
 * The Boldo-Daumas-Li exact argument reduction algorithm.
 *
 * [1] Sylvie Boldo, Marc Daumas, and Ren-Cang Li. "Formally verified argument
 *     reduction with a fused multiply-add."
 *     IEEE Transactions on Computers 58, no. 8 (2008): 1139-1145.
 *     https://arxiv.org/pdf/0708.3722
 */

/*
 * Table I [1]
 */

#define R_PI_FLOAT          (0xa2f983p-25f) // 10680707⋅2⁻²⁵
#define C1_PI_FLOAT         (0xc90fdcp-22f) // 13176796⋅2⁻²²
#define C2_PI_FLOAT         (-0xaeef48p-45f) // -11464520⋅2⁻⁴⁵

#define R_2PI_FLOAT         (R_PI_FLOAT / 2.0f)
#define C1_2PI_FLOAT        (C1_PI_FLOAT * 2.0f)
#define C2_2PI_FLOAT        (C2_PI_FLOAT * 2.0f)

#define R_PI_2_FLOAT        (R_PI_FLOAT * 2.0f)
#define C1_PI_2_FLOAT       (C1_PI_FLOAT / 2.0f)
#define C2_PI_2_FLOAT       (C2_PI_FLOAT / 2.0f)

#define R_PI_DOUBLE         (0x145f306dc9c883p-54) // 5734161139222659⋅2⁻⁵⁴
#define C1_PI_DOUBLE        (0x1921fb54442d18p-51) // 7074237752028440⋅2⁻⁵¹
#define C2_PI_DOUBLE        (0x11a62633145c00p-105) // 4967757600021504⋅2⁻¹⁰⁵

#define R_2PI_DOUBLE        (R_PI_DOUBLE / 2.0)
#define C1_2PI_DOUBLE       (C1_PI_DOUBLE * 2.0)
#define C2_2PI_DOUBLE       (C2_PI_DOUBLE * 2.0)

#define R_PI_2_DOUBLE       (R_PI_DOUBLE * 2.0)
#define C1_PI_2_DOUBLE      (C1_PI_DOUBLE / 2.0)
#define C2_PI_2_DOUBLE      (C2_PI_DOUBLE / 2.0)

/*
 * bias = 3⋅2^(p - N - 2)
 * p: bits of mantissa
 * N: Table I
 */

// p - N - 2 = 24 - 25 - 2 = -3
#define BIAS_FLOAT          (0x3p-3f) // 3⋅2⁻³

// p - N - 2 = 53 - 54 - 2 = -3
#define BIAS_DOUBLE         (0x3p-3) // 3⋅2⁻³

/*
 * See §III. [1]
 */

static inline float
bias_float(int N)
{
    // σ = 3·2^(p − N − 2)
    const int p = 24; // 24 bit mantissa
    return ldexpf(3.0f, p - N - 2);
}

static inline double
bias_double(int N)
{
    // σ = 3·2^(p − N − 2)
    const int p = 53; // 53 bit mantissa
    return ldexp(3.0, p - N - 2);
}

typedef struct {
    float C;
    float C1;
    float C2;
    float R;
    float bias;
} bdl_parameters_float;

typedef struct {
    double C;
    double C1;
    double C2;
    double R;
    double bias;
} bdl_parameters_double;

static const bdl_parameters_float bdl_parameters_pi_float = {
    .C = pi(float),
    .C1 = C1_PI_FLOAT,
    .C2 = C2_PI_FLOAT,
    .R = R_PI_FLOAT,
    .bias = BIAS_FLOAT,
};

static const bdl_parameters_float bdl_parameters_2pi_float = {
    .C = twopi(float),
    .C1 = C1_2PI_FLOAT,
    .C2 = C2_2PI_FLOAT,
    .R = R_2PI_FLOAT,
    .bias = BIAS_FLOAT,
};

static const bdl_parameters_float bdl_parameters_pi_2_float = {
    .C = pi_2(float),
    .C1 = C1_PI_2_FLOAT,
    .C2 = C2_PI_2_FLOAT,
    .R = R_PI_2_FLOAT,
    .bias = BIAS_FLOAT,
};

static const bdl_parameters_double bdl_parameters_pi_double = {
    .C = pi_double,
    .C1 = C1_PI_DOUBLE,
    .C2 = C2_PI_DOUBLE,
    .R = R_PI_DOUBLE,
    .bias = BIAS_DOUBLE,
};

static const bdl_parameters_double bdl_parameters_2pi_double = {
    .C = twopi_double,
    .C1 = C1_2PI_DOUBLE,
    .C2 = C2_2PI_DOUBLE,
    .R = R_2PI_DOUBLE,
    .bias = BIAS_DOUBLE,
};

static const bdl_parameters_double bdl_parameters_pi_2_double = {
    .C = pi_2_double,
    .C1 = C1_PI_2_DOUBLE,
    .C2 = C2_PI_2_DOUBLE,
    .R = R_PI_2_DOUBLE,
    .bias = BIAS_DOUBLE,
};

/*
 * Returns the quotient of x/R.
 * See §III. [1]
 */

static inline int32_t
_bdl_quotient_float(float x, float R, float bias)
{
    // z = fma(x⋅R + σ) − σ
    return (int32_t) (fmaf(x, R, bias) - bias);
}

static inline int64_t
_bdl_quotient_double(double x, double R, double bias)
{
    // z = fma(x⋅R + σ) − σ
    return (int64_t) (fma(x, R, bias) - bias);
}

static inline int32_t
bdl_quotient_float(float x, float R)
{
    assert(bias_float(25) == BIAS_FLOAT);
    // return _bdl_quotient_float(x, R, bias_float(25));
    return _bdl_quotient_float(x, R, BIAS_FLOAT);
}

static inline int64_t
bdl_quotient_double(double x, double R)
{
    assert(bias_double(54) == BIAS_DOUBLE);
    // return _bdl_quotient_double(x, R, bias_double(54));
    return _bdl_quotient_double(x, R, BIAS_DOUBLE);
}

#define bdl_quotient(x, R) \
    _Generic((x), \
        float: bdl_quotient_float, \
        double: bdl_quotient_double \
    )(x, R)

static inline int32_t
bdl_quotient_pi_float(float x)
{
    return bdl_quotient(x, R_PI_FLOAT);
}

static inline int64_t
bdl_quotient_pi_double(double x)
{
    return bdl_quotient(x, R_PI_DOUBLE);
}

#define bdl_quotient_pi(x) \
    _Generic((x), \
        float: bdl_quotient_pi_float, \
        double: bdl_quotient_pi_double \
    )(x)

static inline int32_t
bdl_quotient_2pi_float(float x)
{
    return bdl_quotient(x, R_2PI_FLOAT);
}

static inline int64_t
bdl_quotient_2pi_double(double x)
{
    return bdl_quotient(x, R_2PI_DOUBLE);
}

#define bdl_quotient_2pi(x) \
    _Generic((x), \
        float: bdl_quotient_2pi_float, \
        double: bdl_quotient_2pi_double \
    )(x)

static inline int32_t
bdl_quotient_pi_2_float(float x)
{
    return bdl_quotient(x, R_PI_2_FLOAT);
}

static inline int64_t
bdl_quotient_pi_2_double(double x)
{
    return bdl_quotient(x, R_PI_2_DOUBLE);
}

#define bdl_quotient_pi_2(x) \
    _Generic((x), \
        float: bdl_quotient_pi_2_float, \
        double: bdl_quotient_pi_2_double \
    )(x)

static inline rem_result_float
_bdl_correction_float(bdl_parameters_float parameters, float x, rem_result_float result)
{
    float C = parameters.C;
    float C1 = parameters.C1;
    float C2 = parameters.C2;
    int32_t z = result.z;
    float v1 = result.v1;
    float v2 = result.v2;
    float r = v1 + v2;
    if (x < 0.0f) {
        if (r < -C) {
            // v1 += C;
            v1 += C1;
            v2 += C2;
            z--;
            r = v1 + v2;
        }
        if (r > 0.0f) {
            // v1 -= C;
            v1 -= C1;
            v2 -= C2;
            z++;
            r = v1 + v2;
        }
    } else {
        if (r < 0.0f) {
            // v1 += C;
            v1 += C1;
            v2 += C2;
            z--;
            r = v1 + v2;
        }
        if (r > C) {
            // v1 -= C;
            v1 -= C1;
            v2 -= C2;
            z++;
            r = v1 + v2;
        }
    }
    assert(r >= -C);
    assert(r <= C);
    return (rem_result_float) {z, v1, v2};
}

static inline rem_result_double
_bdl_correction_double(bdl_parameters_double parameters, double x, rem_result_double result)
{
    double C = parameters.C;
    double C1 = parameters.C1;
    double C2 = parameters.C2;
    int64_t z = result.z;
    double v1 = result.v1;
    double v2 = result.v2;
    double r = v1 + v2;
    if (x < 0.0) {
        if (r < -C) {
            // v1 += C;
            v1 += C1;
            v2 += C2;
            z--;
            r = v1 + v2;
        }
        if (r > 0.0) {
            // v1 -= C;
            v1 -= C1;
            v2 -= C2;
            z++;
            r = v1 + v2;
        }
    } else {
        if (r < 0.0) {
            // v1 += C;
            v1 += C1;
            v2 += C2;
            z--;
            r = v1 + v2;
        }
        if (r > C) {
            // v1 -= C;
            v1 -= C1;
            v2 -= C2;
            z++;
            r = v1 + v2;
        }
    }
    assert(r >= -C);
    assert(r <= C);
    return (rem_result_double) {z, v1, v2};
}

#define _bdl_correction(parameters, x, result) \
    _Generic((parameters), \
        bdl_parameters_float: _bdl_correction_float, \
        bdl_parameters_double: _bdl_correction_double \
    )(parameters, x, result)

/*
 * Reduce an argument x to the range [0, C], given C = C₁ + C₂, R = 1∕C,
 * C₁ and C₂ meet the requirement of Theorem 6, using Algorithm 5.1. [1]
 * Returns (z, v₁, v₂) such that x = z⋅(C₁ + C₂) + v₁ + v₂.
 */

static inline rem_result_float
bdl_reduce_float(bdl_parameters_float parameters, float x)
{
    float C1 = parameters.C1;
    float C2 = parameters.C2;
    float R = parameters.R;

    int32_t z = bdl_quotient(x, R);
    float u = fmaf(-z, C1, x);
    float v1 = fmaf(-z, C2, u);
    tuple_float p = mul((float) z, C2);
    tuple_float t = add(u, -p.a);
    float v2 = ((t.a - v1) + t.b) - p.b;

    rem_result_float result = {z, v1, v2};
    /*
     * Algorithm 5.1 guarantees v₁+v₂ = x - z⋅(C₁ + C₂) but
     * z may be off by one or two.
     * Correct {z, v1, v2} such that 0 ≤ |v₁+v₂| ≤ C.
     */
    result = _bdl_correction(parameters, x, result);
    return result;
}

static inline rem_result_double
bdl_reduce_double(bdl_parameters_double parameters, double x)
{
    double C1 = parameters.C1;
    double C2 = parameters.C2;
    double R = parameters.R;

    int64_t z = bdl_quotient(x, R);
    double u = fma(-z, C1, x);
    double v1 = fma(-z, C2, u);
    tuple_double p = mul((double) z, C2);
    tuple_double t = add(u, -p.a);
    double v2 = ((t.a - v1) + t.b) - p.b;

    rem_result_double result = {z, v1, v2};
    /*
     * Algorithm 5.1 guarantees v₁+v₂ = x - z⋅(C₁ + C₂) but
     * z may be off by one or two.
     * Correct {z, v1, v2} such that 0 ≤ |v₁+v₂| ≤ C.
     */
    if (x <= -0x1.0p24 || x >= 0x1.0p24)
        result = _bdl_correction(parameters, x, result);
    return result;
}

#define bdl_reduce(parameters, x) \
    _Generic((parameters), \
        bdl_parameters_float: bdl_reduce_float, \
        bdl_parameters_double: bdl_reduce_double \
    )(parameters, x)

/*
 * Return the quotient z and remainder (v₁+v₂) of x∕C
 * as {z, v₁, v₂}, such that x = z⋅C + v₁+v₂ and |v₁+v₂| ≤ C.
 * If x < 0, then z < 0 and v₁+v₂ < 0.
 *
 * If x > 2⁵³, the identity x = z⋅C + v₁+v₂ is guaranteed,
 * but not |v₁+v₂| ≤ C.
 *
 * If x > 2⁶³, this method fails.
 */

static inline rem_result_float
_bdl_rem_pi_float(float x)
{
    return bdl_reduce(bdl_parameters_pi_float, x);
}

static inline rem_result_double
_bdl_rem_pi_double(double x)
{
    return bdl_reduce(bdl_parameters_pi_double, x);
}

#define _bdl_rem_pi(x) \
    _Generic((x), \
        float: _bdl_rem_pi_float, \
        double: _bdl_rem_pi_double \
    )(x)

static inline rem_result_float
_bdl_rem_2pi_float(float x)
{
    return bdl_reduce(bdl_parameters_2pi_float, x);
}

static inline rem_result_double
_bdl_rem_2pi_double(double x)
{
    return bdl_reduce(bdl_parameters_2pi_double, x);
}

#define _bdl_rem_2pi(x) \
    _Generic((x), \
        float: _bdl_rem_2pi_float, \
        double: _bdl_rem_2pi_double \
    )(x)

static inline rem_result_float
_bdl_rem_pi_2_float(float x)
{
    return bdl_reduce(bdl_parameters_pi_2_float, x);
}

static inline rem_result_double
_bdl_rem_pi_2_double(double x)
{
    return bdl_reduce(bdl_parameters_pi_2_double, x);
}

#define _bdl_rem_pi_2(x) \
    _Generic((x), \
        float: _bdl_rem_pi_2_float, \
        double: _bdl_rem_pi_2_double \
    )(x)

rem_result_float naive_rem_pi_float(float x);
rem_result_float naive_rem_2pi_float(float x);
rem_result_float naive_rem_pi_2_float(float x);
rem_result_double naive_rem_pi_double(double x);
rem_result_double naive_rem_2pi_double(double x);
rem_result_double naive_rem_pi_2_double(double x);

#define naive_rem_pi(x) \
    _Generic((x), \
        float: naive_rem_pi_float, \
        double: naive_rem_pi_double \
    )(x)
#define naive_rem_2pi(x) \
    _Generic((x), \
        float: naive_rem_2pi_float, \
        double: naive_rem_2pi_double \
    )(x)
#define naive_rem_pi_2(x) \
    _Generic((x), \
        float: naive_rem_pi_2_float, \
        double: naive_rem_pi_2_double \
    )(x)

rem_result_float bdl_rem_pi_float(float x);
rem_result_float bdl_rem_2pi_float(float x);
rem_result_float bdl_rem_pi_2_float(float x);
rem_result_double bdl_rem_pi_double(double x);
rem_result_double bdl_rem_2pi_double(double x);
rem_result_double bdl_rem_pi_2_double(double x);

#define bdl_rem_pi(x) \
    _Generic((x), \
        float: bdl_rem_pi_float, \
        double: bdl_rem_pi_double \
    )(x)
#define bdl_rem_2pi(x) \
    _Generic((x), \
        float: bdl_rem_2pi_float, \
        double: bdl_rem_2pi_double \
    )(x)
#define bdl_rem_pi_2(x) \
    _Generic((x), \
        float: bdl_rem_pi_2_float, \
        double: bdl_rem_pi_2_double \
    )(x)
