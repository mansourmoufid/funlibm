/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#pragma once

#include <stdint.h> // int32_t, int64_t

typedef struct {
    float a;
    float b;
} tuple_float;

typedef struct {
    double a;
    double b;
} tuple_double;

typedef struct {
    int32_t z;
    float v1;
    float v2;
} rem_result_float;

typedef struct {
    int64_t z;
    double v1;
    double v2;
} rem_result_double;
