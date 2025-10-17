/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#pragma once

#include <stddef.h> // size_t

float _sin_float(float x);
double _sin_double(double x);
#define _sin(x) \
    _Generic((x), \
        float: _sin_float, \
        double: _sin_double \
    )(x)
void _sin_array_float(float xs[], float ys[], size_t n);
void _sin_array_double(double xs[], double ys[], size_t n);
#define _sin_array(xs, ys, n) \
    _Generic((xs[0]), \
        float: _sin_array_float, \
        double: _sin_array_double \
    )(xs, ys, n)

float _cos_float(float x);
double _cos_double(double x);
#define _cos(x) \
    _Generic((x), \
        float: _cos_float, \
        double: _cos_double \
    )(x)
void _cos_array_float(float xs[], float ys[], size_t n);
void _cos_array_double(double xs[], double ys[], size_t n);
#define _cos_array(xs, ys, n) \
    _Generic((xs[0]), \
        float: _cos_array_float, \
        double: _cos_array_double \
    )(xs, ys, n)
