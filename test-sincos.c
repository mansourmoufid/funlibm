/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#include <assert.h>
#include <math.h> // M_PI
#include <stddef.h> // size_t, NULL
#include <stdlib.h> // calloc, drand48
#include <stdio.h> // printf
#include <time.h>

#include <mpfr.h>

#include "common.h" // significant_digits, ulp, xtime
#include "sincos.h"
#include "types.h" // tuple_double

static const mpfr_prec_t mp_precision = 128;

static int indent = 0;

// The size of the precision test array.
static const size_t n = 10000000;

static void
test_float(
    const char *function_name,
    float (*function)(float),
    int (*mpfr_function)(mpfr_t rop, const mpfr_t op, mpfr_rnd_t rnd)
) {
    printf("testing float %s ...\n\n", function_name);
    indent += 4;

    float *xs = calloc(n, sizeof (float)); // array of input values x
    float *ys = calloc(n, sizeof (float)); // array of function(x)
    float *zs = calloc(n, sizeof (float)); // array of function(x) using mpfr
    assert(xs != NULL);
    assert(ys != NULL);
    assert(zs != NULL);

    srand48(time(NULL));

    tuple_float D; // test domain
    D = (tuple_float) {-2.0 * M_PI, 2.0 * M_PI};
    for (size_t i = 0; i < n / 2; i++) xs[i] = drand48() * (D.b - D.a) + D.a;
    // D = (tuple_float) {-0x1.0p23, 0x1.0p23};
    for (size_t i = n / 2; i < n; i++) xs[i] = drand48() * (D.b - D.a) + D.a;

    for (size_t i = 0; i < n; i++)
        ys[i] = (*function)(xs[i]);

    mpfr_t mp_x;
    mpfr_t mp_result;
    mpfr_t mp_error;
    mpfr_init2(mp_x, mp_precision);
    mpfr_init2(mp_result, mp_precision);
    mpfr_init2(mp_error, mp_precision);
    float max_abs_error = 0.0;
    float max_rel_error = 0.0;
    int error_dist[4] = {0};
    for (size_t i = 0; i < n; i++) {
        mpfr_set_flt(mp_x, xs[i], MPFR_RNDN);
        (*mpfr_function)(mp_result, mp_x, MPFR_RNDN);
        zs[i] = mpfr_get_flt(mp_result, MPFR_RNDN);

        // error = |mpfr_function(x) - function(x)|
        mpfr_sub_d(mp_error, mp_result, ys[i], MPFR_RNDN);
        mpfr_abs(mp_error, mp_error, MPFR_RNDN);
        float abs_error = mpfr_get_flt(mp_error, MPFR_RNDN);
        float rel_error = abs_error / ulp(zs[i]);
        if (rel_error >= 1.0f) {
            int N = significant_digits(zs[i]);
            fprintf(stderr, "%*sx = %+.*f\n", indent, "", 20, xs[i]);
            fprintf(stderr, "%*s    expected %s(x) = %+.*f\n", indent, "", function_name, N, zs[i]);
            fprintf(stderr, "%*s         got %s(x) = %+.*f\n", indent, "", function_name, N, ys[i]);
            fprintf(stderr, "%*s              error = %+.*f = %e\n", indent, "", N, abs_error, abs_error);
            fprintf(stderr, "%*s        ulp(%s(x)) = %+.*f = %e\n", indent, "", function_name, N, ulp(zs[i]), ulp(zs[i]));
            fprintf(stderr, "\n");
        }
        if (rel_error >= 3.0f)
            error_dist[3]++;
        else if (rel_error >= 2.0f)
            error_dist[2]++;
        else if (rel_error >= 1.0f)
            error_dist[1]++;
        else
            error_dist[0]++;
        if (abs_error > max_abs_error)
            max_abs_error = abs_error;
        if (rel_error > max_rel_error)
            max_rel_error = rel_error;
    }
    assert(error_dist[0] + error_dist[1] + error_dist[2] + error_dist[3] == n);
    printf("%*serror distribution:\n", indent, "");
    printf("%*s 0 ulp %i (%.2f%%)\n", indent, "", error_dist[0], (float) error_dist[0] / n * 100.0f);
    printf("%*s 1 ulp %i (%.2f%%)\n", indent, "", error_dist[1], (float) error_dist[1] / n * 100.0f);
    printf("%*s 2 ulp %i (%.2f%%)\n", indent, "", error_dist[2], (float) error_dist[2] / n * 100.0f);
    printf("%*s≥3 ulp %i (%.2f%%)\n", indent, "", error_dist[3], (float) error_dist[3] / n * 100.0f);
    printf("\n");

    mpfr_clear(mp_error);
    mpfr_clear(mp_result);
    mpfr_clear(mp_x);

    free(xs);
    free(ys);
    free(zs);

    indent -= 4;
}

static void
test_double(
    const char *function_name,
    double (*function)(double),
    int (*mpfr_function)(mpfr_t rop, const mpfr_t op, mpfr_rnd_t rnd)
) {
    printf("testing double %s ...\n\n", function_name);
    indent += 4;

    double *xs = calloc(n, sizeof (double)); // array of input values x
    double *ys = calloc(n, sizeof (double)); // array of function(x)
    double *zs = calloc(n, sizeof (double)); // array of function(x) using mpfr
    assert(xs != NULL);
    assert(ys != NULL);
    assert(zs != NULL);

    srand48(time(NULL));

    tuple_double D; // test domain
    D = (tuple_double) {-2.0 * M_PI, 2.0 * M_PI};
    for (size_t i = 0; i < n / 2; i++) xs[i] = drand48() * (D.b - D.a) + D.a;
    // D = (tuple_double) {-0x1.0p52, 0x1.0p52};
    for (size_t i = n / 2; i < n; i++) xs[i] = drand48() * (D.b - D.a) + D.a;

    for (size_t i = 0; i < n; i++)
        ys[i] = (*function)(xs[i]);

    mpfr_t mp_x;
    mpfr_t mp_result;
    mpfr_t mp_error;
    mpfr_init2(mp_x, mp_precision);
    mpfr_init2(mp_result, mp_precision);
    mpfr_init2(mp_error, mp_precision);
    double max_abs_error = 0.0;
    double max_rel_error = 0.0;
    int error_dist[4] = {0};
    for (size_t i = 0; i < n; i++) {
        mpfr_set_d(mp_x, xs[i], MPFR_RNDN);
        (*mpfr_function)(mp_result, mp_x, MPFR_RNDN);
        zs[i] = mpfr_get_d(mp_result, MPFR_RNDN);

        // error = |mpfr_function(x) - function(x)|
        mpfr_sub_d(mp_error, mp_result, ys[i], MPFR_RNDN);
        mpfr_abs(mp_error, mp_error, MPFR_RNDN);
        double abs_error = mpfr_get_d(mp_error, MPFR_RNDN);
        double rel_error = abs_error / ulp(zs[i]);
        if (rel_error >= 1.0) {
            int N = significant_digits(zs[i]);
            fprintf(stderr, "%*sx = %+.*f\n", indent, "", 20, xs[i]);
            fprintf(stderr, "%*s    expected %s(x) = %+.*f\n", indent, "", function_name, N, zs[i]);
            fprintf(stderr, "%*s         got %s(x) = %+.*f\n", indent, "", function_name, N, ys[i]);
            fprintf(stderr, "%*s              error = %+.*f = %e\n", indent, "", N, abs_error, abs_error);
            fprintf(stderr, "%*s        ulp(%s(x)) = %+.*f = %e\n", indent, "", function_name, N, ulp(zs[i]), ulp(zs[i]));
            fprintf(stderr, "\n");
        }
        if (rel_error >= 3.0)
            error_dist[3]++;
        else if (rel_error >= 2.0)
            error_dist[2]++;
        else if (rel_error >= 1.0)
            error_dist[1]++;
        else
            error_dist[0]++;
        if (abs_error > max_abs_error)
            max_abs_error = abs_error;
        if (rel_error > max_rel_error)
            max_rel_error = rel_error;
    }
    assert(error_dist[0] + error_dist[1] + error_dist[2] + error_dist[3] == n);
    printf("%*serror distribution:\n", indent, "");
    printf("%*s 0 ulp %i (%.2f%%)\n", indent, "", error_dist[0], (double) error_dist[0] / n * 100.0);
    printf("%*s 1 ulp %i (%.2f%%)\n", indent, "", error_dist[1], (double) error_dist[1] / n * 100.0);
    printf("%*s 2 ulp %i (%.2f%%)\n", indent, "", error_dist[2], (double) error_dist[2] / n * 100.0);
    printf("%*s≥3 ulp %i (%.2f%%)\n", indent, "", error_dist[3], (double) error_dist[3] / n * 100.0);
    printf("\n");

    mpfr_clear(mp_error);
    mpfr_clear(mp_result);
    mpfr_clear(mp_x);

    free(xs);
    free(ys);
    free(zs);

    indent -= 4;
}

int
main(void)
{
    // test_float("libm sinf", &sinf, &mpfr_sin);
    // test_float("libm cosf", &cosf, &mpfr_cos);
    test_float("sin", &_sin_float, &mpfr_sin);
    test_float("cos", &_cos_float, &mpfr_cos);
    test_double("sin", &_sin_double, &mpfr_sin);
    test_double("cos", &_cos_double, &mpfr_cos);
    printf("\n");
    return 0;
}
