/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

/* Compile with -DNDEBUG, except for debugging. */

#pragma once

#include <assert.h>
#include <math.h> // fma

#include "types.h" // tuple_double

/*
 * [1] Marc Daumas, Laurence Rideau, Laurent Thery. A Generic Library for
 *     Floating-Point Numbers and Its Application to Exact Computing.
 *     Theorem Proving in Higher Order Logics, 2001, Edinburgh, United Kingdom.
 *     pp.169-184. https://hal.science/hal-00157285
 *
 * [2] Alan H. Karp and Peter Markstein. 1997. High-precision division and
 *     square root. ACM Trans. Math. Softw. 23, 4 (Dec. 1997), 561–589.
 *     https://dl.acm.org/doi/pdf/10.1145/279232.279237
 */

/*
 * Return the sum and its error. See TwoSum_k, page 178. [1]
 */

static inline tuple_float
__attribute__((always_inline))
__attribute__((const))
add_float(float a, float b)
{
    float x = a + b;
    float e = (a - (x - b)) + (b - (x - (x - b)));
    return (tuple_float) {x, e};
}

static inline tuple_double
__attribute__((always_inline))
__attribute__((const))
add_double(double a, double b)
{
    double x = a + b;
    double e = (a - (x - b)) + (b - (x - (x - b)));
    return (tuple_double) {x, e};
}

#define add(a, b) \
    _Generic((a), \
        float: add_float, \
        double: add_double \
    )(a, b)

/*
 * Return the product and its error. See Figure 2, page 566. [2]
 */

static inline tuple_float
__attribute__((always_inline))
__attribute__((const))
mul_float(float a, float b)
{
    float x = a * b;
    return (tuple_float) {x, fmaf(a, b, -x)};
}

static inline tuple_double
__attribute__((always_inline))
__attribute__((const))
mul_double(double a, double b)
{
    double x = a * b;
    return (tuple_double) {x, fma(a, b, -x)};
}

#define mul(a, b) \
    _Generic((a), \
        float: mul_float, \
        double: mul_double \
    )(a, b)

/*
 * Return a * b + c.
 */

static inline tuple_float
__attribute__((always_inline))
__attribute__((const))
xfma_float(float a, float b, float c)
{
    tuple_float p = mul(a, b);
    tuple_float s1 = add(c, p.a);
    tuple_float s2  = add(s1.a, p.b);
    return (tuple_float) {s2.a, s1.b + s2.b};
}

static inline tuple_double
__attribute__((always_inline))
__attribute__((const))
xfma_double(double a, double b, double c)
{
    tuple_double p = mul(a, b);
    tuple_double s1 = add(c, p.a);
    tuple_double s2  = add(s1.a, p.b);
    return (tuple_double) {s2.a, s1.b + s2.b};
}

#define xfma(a, b, c) \
    _Generic((a), \
        float: xfma_float, \
        double: xfma_double \
    )(a, b, c)
