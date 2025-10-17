/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#include <assert.h>
#include <math.h> // fmaf
#include <stddef.h> // size_t
#include <stdint.h> // int64_t
#include <stdio.h>

#include "constants.h" // pi, etc.
#include "cw.h" // _cw_rem_pi_2, etc.
#include "reduce.h" // _bdl_rem_pi_2, etc.
#include "sincos.h"
#include "types.h" // rem_result_double

#include "cos.h"
#include "sin.h"

#define rem_pi_2 _cw_rem_pi_2

static inline float
__attribute__((always_inline))
__attribute__((const))
_sin_pi_2_float(float x)
{
    float sign = 1.0f;
    if (x < 0.0f) { // sin(-x) = -sin(x)
        x = -x;
        sign = -1.0f;
    }
    assert(x >= 0.0f);
    if (x <= 2.7e-4f)
        return sign * x;
    int32_t q = 0;
    float r = x;
    float v1 = r;
    float v2 = 0.0f;
    if (x > pi_2(float)) {
        rem_result_float result = rem_pi_2(x);
        q = result.z % 4;
        v1 = result.v1;
        v2 = result.v2;
        r = result.v1 + result.v2;
    }
    assert(q >= 0);
    assert(q < 4);
    assert(r >= 0.0f);
    assert(r <= pi_2(float));
    float s;
    if (q <= 0) { // '<=' rather than '==' otherwise llvm won't vectorize
        s = _sin_0_pi_2(v1);
        if (v2 != 0.0f) {
            // Newton-Raphson
            // sin(x + 𝜹) ≅ sin(x) + sinʹ(x)⋅𝜹 = sin(x) + cos(x)⋅𝜹
            // sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v₂
            s = s + _cos_0_pi_2(v1) * v2;
        }
    } else if (q <= 1) {
        // sin(π∕2 + x) = cos(x)
        s = _cos_0_pi_2(v1);
        if (v2 != 0.0f) {
            // Newton-Raphson
            // cos(x + 𝜹) ≅ cos(x) + cosʹ(x)⋅𝜹 = cos(x) - sin(x)⋅𝜹
            // cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            // sin(π∕2 + v₁ + v₂) = cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            s = s - _sin_0_pi_2(v1) * v2;
        }
    } else if (q <= 2) {
        // sin(π + x) = -sin(x)
        s = -_sin_0_pi_2(v1);
        if (v2 != 0.0f) {
            // Newton-Raphson
            // sin(x + 𝜹) ≅ sin(x) + sinʹ(x)⋅𝜹 = sin(x) + cos(x)⋅𝜹
            // sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v₂
            // sin(π + v₁ + v₂) = -sin(v₁ + v₂) ≅ -sin(v₁) - cos(v₁)⋅v₂
            s = s - _cos_0_pi_2(v1) * v2;
        }
    } else {
        // sin(3π∕2 + x) = -cos(x)
        s = -_cos_0_pi_2(v1);
        if (v2 != 0.0f) {
            // Newton-Raphson
            // cos(x + 𝜹) ≅ cos(x) + cosʹ(x)⋅𝜹 = cos(x) - sin(x)⋅𝜹
            // cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            // sin(3π∕2 + v₁ + v₂) = -cos(v₁ + v₂) ≅ -cos(v₁) + sin(v₁)⋅v2
            s = s + _sin_0_pi_2(v1) * v2;
        }
    }
    return sign * s;
}

static inline double
__attribute__((always_inline))
__attribute__((const))
_sin_pi_2_double(double x)
{
    double sign = 1.0;
    if (x < 0.0) { // sin(-x) = -sin(x)
        x = -x;
        sign = -1.0;
    }
    assert(x >= 0.0);
    if (x <= 2.1e-8f)
        return sign * x;
    int64_t q = 0;
    double r = x;
    double v1 = r;
    double v2 = 0.0;
    if (x > pi_2(double)) {
        rem_result_double result = rem_pi_2(x);
        q = result.z % 4;
        v1 = result.v1;
        v2 = result.v2;
        r = result.v1 + result.v2;
    }
    assert(q >= 0);
    assert(q < 4);
    assert(r >= 0.0);
    assert(r <= pi_2(double));
    double s;
    if (q <= 0) { // '<=' rather than '==' otherwise llvm won't vectorize
        s = _sin_0_pi_2(v1);
        if (v2 != 0.0) {
            // Newton-Raphson
            // sin(x + 𝜹) ≅ sin(x) + sinʹ(x)⋅𝜹 = sin(x) + cos(x)⋅𝜹
            // sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v₂
            s = s + _cos_0_pi_2(v1) * v2;
        }
    } else if (q <= 1) {
        // sin(π∕2 + x) = cos(x)
        s = _cos_0_pi_2(v1);
        if (v2 != 0.0) {
            // Newton-Raphson
            // cos(x + 𝜹) ≅ cos(x) + cosʹ(x)⋅𝜹 = cos(x) - sin(x)⋅𝜹
            // cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            // sin(π∕2 + v₁ + v₂) = cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            s = s - _sin_0_pi_2(v1) * v2;
        }
    } else if (q <= 2) {
        // sin(π + x) = -sin(x)
        s = -_sin_0_pi_2(v1);
        if (v2 != 0.0) {
            // Newton-Raphson
            // sin(x + 𝜹) ≅ sin(x) + sinʹ(x)⋅𝜹 = sin(x) + cos(x)⋅𝜹
            // sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v₂
            // sin(π + v₁ + v₂) = -sin(v₁ + v₂) ≅ -sin(v₁) - cos(v₁)⋅v₂
            s = s - _cos_0_pi_2(v1) * v2;
        }
    } else {
        // sin(3π∕2 + x) = -cos(x)
        s = -_cos_0_pi_2(v1);
        if (v2 != 0.0) {
            // Newton-Raphson
            // cos(x + 𝜹) ≅ cos(x) + cosʹ(x)⋅𝜹 = cos(x) - sin(x)⋅𝜹
            // cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            // sin(3π∕2 + v₁ + v₂) = -cos(v₁ + v₂) ≅ -cos(v₁) + sin(v₁)⋅v2
            s = s + _sin_0_pi_2(v1) * v2;
        }
    }
    return sign * s;
}

#define _sin_pi_2(x) \
    _Generic((x), \
        float: _sin_pi_2_float, \
        double: _sin_pi_2_double \
    )(x)

float
_sin_float(float x)
{
    return _sin_pi_2(x);
}

double
_sin_double(double x)
{
    return _sin_pi_2(x);
}

void
_sin_array_float(float xs[], float ys[], size_t n)
{
    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < n; i++)
        ys[i] = _sin_pi_2(xs[i]);
}

void
_sin_array_double(double xs[], double ys[], size_t n)
{
    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < n; i++)
        ys[i] = _sin_pi_2(xs[i]);
}

static inline float
__attribute__((always_inline))
__attribute__((const))
_cos_pi_2_float(float x)
{
    if (x < 0.0f) // cos(-x) = cos(x)
        x = -x;
    assert(x >= 0.0);
    int32_t q = 0;
    float r = x;
    float v1 = r;
    float v2 = 0.0;
    if (x > pi_2(float)) {
        rem_result_float result = rem_pi_2(x);
        q = result.z % 4;
        v1 = result.v1;
        v2 = result.v2;
        r = result.v1 + result.v2;
    }
    assert(q >= 0);
    assert(q < 4);
    assert(r >= 0.0f);
    assert(r <= pi_2(float));
    float c;
    if (q <= 0) { // '<=' rather than '==' otherwise llvm won't vectorize
        c = _cos_0_pi_2(v1);
        if (v2 != 0.0f) {
            // Newton-Raphson
            // cos(x + 𝜹) ≅ cos(x) + cosʹ(x)⋅𝜹 = cos(x) - sin(x)⋅𝜹
            // cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            c = c - _sin_0_pi_2(v1) * v2;
        }
    } else if (q <= 1) {
        // cos(π∕2 + x) = -sin(x)
        c = -_sin_0_pi_2(v1);
        if (v2 != 0.0f) {
            // Newton-Raphson
            // sin(x + 𝜹) ≅ sin(x) + sinʹ(x)⋅𝜹 = sin(x) + cos(x)⋅𝜹
            // sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v₂
            // cos(π∕2 + v₁ + v₂) = -sin(v₁ + v₂) ≅ -sin(v₁) - cos(v₁)⋅v₂
            c = c - _cos_0_pi_2(v1) * v2;
        }
    } else if (q <= 2) {
        // cos(π + x) = -cos(x)
        c = -_cos_0_pi_2(v1);
        if (v2 != 0.0f) {
            // Newton-Raphson
            // cos(x + 𝜹) ≅ cos(x) + cosʹ(x)⋅𝜹 = cos(x) - sin(x)⋅𝜹
            // cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            // cos(π + v₁ + v₂) = -cos(v₁ + v₂) ≅ -cos(v₁) + sin(v₁)⋅v₂
            c = c + _sin_0_pi_2(v1) * v2;
        }
    } else {
        // cos(3π∕2 + x) = sin(x)
        c = _sin_0_pi_2(v1);
        if (v2 != 0.0f) {
            // Newton-Raphson
            // sin(x + 𝜹) ≅ sin(x) + sinʹ(x)⋅𝜹 = sin(x) + cos(x)⋅𝜹
            // sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v₂
            // cos(3π∕2 + v₁ + v₂) = sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v2
            c = c + _cos_0_pi_2(v1) * v2;
        }
    }
    return c;
}

static inline double
__attribute__((always_inline))
__attribute__((const))
_cos_pi_2_double(double x)
{
    if (x < 0.0) // cos(-x) = cos(x)
        x = -x;
    assert(x >= 0.0);
    int64_t q = 0;
    double r = x;
    double v1 = r;
    double v2 = 0.0;
    if (x > pi_2(double)) {
        rem_result_double result = rem_pi_2(x);
        q = result.z % 4;
        v1 = result.v1;
        v2 = result.v2;
        r = result.v1 + result.v2;
    }
    assert(q >= 0);
    assert(q < 4);
    assert(r >= 0.0);
    assert(r <= pi_2(double));
    double c;
    if (q <= 0) { // '<=' rather than '==' otherwise llvm won't vectorize
        c = _cos_0_pi_2(v1);
        if (v2 != 0.0) {
            // Newton-Raphson
            // cos(x + 𝜹) ≅ cos(x) + cosʹ(x)⋅𝜹 = cos(x) - sin(x)⋅𝜹
            // cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            c = c - _sin_0_pi_2(v1) * v2;
        }
    } else if (q <= 1) {
        // cos(π∕2 + x) = -sin(x)
        c = -_sin_0_pi_2(v1);
        if (v2 != 0.0) {
            // Newton-Raphson
            // sin(x + 𝜹) ≅ sin(x) + sinʹ(x)⋅𝜹 = sin(x) + cos(x)⋅𝜹
            // sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v₂
            // cos(π∕2 + v₁ + v₂) = -sin(v₁ + v₂) ≅ -sin(v₁) - cos(v₁)⋅v₂
            c = c - _cos_0_pi_2(v1) * v2;
        }
    } else if (q <= 2) {
        // cos(π + x) = -cos(x)
        c = -_cos_0_pi_2(v1);
        if (v2 != 0.0) {
            // Newton-Raphson
            // cos(x + 𝜹) ≅ cos(x) + cosʹ(x)⋅𝜹 = cos(x) - sin(x)⋅𝜹
            // cos(v₁ + v₂) ≅ cos(v₁) - sin(v₁)⋅v₂
            // cos(π + v₁ + v₂) = -cos(v₁ + v₂) ≅ -cos(v₁) + sin(v₁)⋅v₂
            c = c + _sin_0_pi_2(v1) * v2;
        }
    } else {
        // cos(3π∕2 + x) = sin(x)
        c = _sin_0_pi_2(v1);
        if (v2 != 0.0) {
            // Newton-Raphson
            // sin(x + 𝜹) ≅ sin(x) + sinʹ(x)⋅𝜹 = sin(x) + cos(x)⋅𝜹
            // sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v₂
            // cos(3π∕2 + v₁ + v₂) = sin(v₁ + v₂) ≅ sin(v₁) + cos(v₁)⋅v2
            c = c + _cos_0_pi_2(v1) * v2;
        }
    }
    return c;
}

#define _cos_pi_2(x) \
    _Generic((x), \
        float: _cos_pi_2_float, \
        double: _cos_pi_2_double \
    )(x)

float
_cos_float(float x)
{
    return _cos_pi_2(x);
}

double
_cos_double(double x)
{
    return _cos_pi_2(x);
}

void
_cos_array_float(float xs[], float ys[], size_t n)
{
    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < n; i++)
        ys[i] = _cos_pi_2(xs[i]);
}

void
_cos_array_double(double xs[], double ys[], size_t n)
{
    #pragma clang loop vectorize(enable)
    for (size_t i = 0; i < n; i++)
        ys[i] = _cos_pi_2(xs[i]);
}
