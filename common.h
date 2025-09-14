/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#pragma once

#include <assert.h>
#include <math.h> // fabs, frexp, ldexp
#include <stddef.h> // NULL
#include <sys/time.h> // struct timeval, gettimeofday

static float
ulp_float(float x)
{
    int exp = 0;
    frexpf(fabsf(x), &exp);
    return ldexpf(1.0f, exp - 24);
}

static double
ulp_double(double x)
{
    int exp = 0;
    frexp(fabs(x), &exp);
    return ldexp(1.0, exp - 53);
}

#define ulp(x) \
    _Generic((x), \
        float: ulp_float, \
        double: ulp_double \
    )(x)

static inline int
significant_digits_float(float x)
{
    // at least 2 significant digits of the ulp
    return 2 - (int) log10f(ulp(x));
}

static inline int
significant_digits_double(double x)
{
    // at least 2 significant digits of the ulp
    return 2 - (int) log10(ulp(x));
}

#define significant_digits(x) \
    _Generic((x), \
        float: significant_digits_float, \
        double: significant_digits_double \
    )(x)

static long int
xtime(void)
{
    struct timeval now = {0};
    (void) gettimeofday(&now, NULL);
    return now.tv_sec * 1000000L + now.tv_usec; // microseconds
}
