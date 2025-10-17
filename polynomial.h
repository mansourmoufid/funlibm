/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#pragma once

#include <math.h> // fma
#include <stddef.h> // size_t

#include "arithmetic.h" // add, mul
#include "types.h" // tuple_double, tuple_float

static inline float
__attribute__((always_inline))
__attribute__((const))
eval_polynomial_float1(const float as[], size_t n, float x)
{
    float r = as[n - 1];
    #pragma unroll
    for (size_t i = 1; i < n; i++)
        r = fmaf(r, x, as[n - 1 - i]);
    return r;
}

/*
 * Implementation of the Graillat–Langlois–Louvet error-free polynomial
 * evaluation algorithm. [1]
 *
 * [1] S. Graillat, P. Langlois, and N. Louvet. Algorithms for accurate,
 * validated and fast computations with polynomials. Japan Journal of
 * Industrial and Applied Mathematics, Special issue on Verified Numerical
 * Computation, 2009.
 */

static inline float
__attribute__((always_inline))
__attribute__((const))
eval_polynomial_float(const float as[], size_t n, float x)
{
    float r; // result
    tuple_float p; // product
    tuple_float s; // sum
    float pe; // product error
    float se; // sum error
    float e = 0.0; // total error
    r = as[n - 1];
    #pragma unroll
    for (size_t i = 1; i < n; i++) {
        p = mul(r, x);
        r = p.a;
        pe = p.b;
        s = add(r, as[n - 1 - i]);
        r = s.a;
        se = s.b;
        e = fmaf(e, x, pe + se);
    }
    return r + e;
}

static inline double
__attribute__((always_inline))
__attribute__((const))
eval_polynomial_double(const double as[], size_t n, double x)
{
    double r; // result
    tuple_double p; // product
    tuple_double s; // sum
    double pe; // product error
    double se; // sum error
    double e = 0.0; // total error
    r = as[n - 1];
    #pragma unroll
    for (size_t i = 1; i < n; i++) {
        p = mul(r, x);
        r = p.a;
        pe = p.b;
        s = add(r, as[n - 1 - i]);
        r = s.a;
        se = s.b;
        e = fma(e, x, pe + se);
    }
    return r + e;
}

#define eval_polynomial(as, n, x) \
    _Generic((as)[0], \
        float: eval_polynomial_float, \
        double: eval_polynomial_double \
    )(as, n, x)
