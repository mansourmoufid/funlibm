/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#include <assert.h>
#include <math.h> // fabs, frexp, ldexp
#include <stdio.h>
#include <stddef.h> // size_t
#include <stdint.h> // int64_t
#include <stdlib.h> // drand48
#include <string.h> // strcmp

#include <mpfr.h>

#include "common.h" // ulp
#include "reduce.h"
#include "types.h" // rem_result_double, rem_result_float

static const mpfr_prec_t mp_precision = 2 * (sizeof(double) * 8);

static const size_t n = 1000000;

static int indent = 0;

static void
test_q_and_r_float(rem_result_float (*rem_function)(float), mpfr_t mp_divisor)
{
    size_t error_dist[4] = {0};
    const float magnitudes[] = {
        0x1.0p+00f,
        0x1.0p+01f,
        0x1.0p+02f,
        0x1.0p+04f,
        0x1.0p+08f,
        0x1.0p+16f,
        0x1.0p+20f,
        0x1.0p+23f,
    };
    float max_abs_error_by_magnitude[24] = {0.0f};
    float max_rel_error_by_magnitude[24] = {0.0f};
    for (size_t j = 0; j < sizeof magnitudes / sizeof magnitudes[0]; j++) {
        float max_abs_error = 0.0f;
        float max_rel_error = 0.0f;
        for (size_t i = 0; i < n; i++) {
            float x = (drand48() - 0.5) * 2.0f * magnitudes[j];
            mpfr_t mp_x;
            mpfr_init2(mp_x, mp_precision);
            mpfr_set_flt(mp_x, x, MPFR_RNDN);

            rem_result_float rem = (*rem_function)(x);
            int32_t q = rem.z;
            float r = rem.v1 + rem.v2;
            float abs_error = 0.0f;
            float rel_error = 0.0f;

            // mp_q = trunc(mp_x / mp_divisor)
            mpfr_t mp_q;
            mpfr_init2(mp_q, mp_precision);
            mpfr_div(mp_q, mp_x, mp_divisor, MPFR_RNDN);
            mpfr_trunc(mp_q, mp_q);
            long int mp_q_int = mpfr_get_si(mp_q, MPFR_RNDN);
            if (q != mp_q_int) {
                fprintf(stderr, "%*sx = %.12f\n", indent, "", x);
                fprintf(
                    stderr,
                    "%*squotient of %.12f ÷ %.12f = %.0f\n",
                    indent, "",
                    mpfr_get_flt(mp_x, MPFR_RNDN),
                    mpfr_get_flt(mp_divisor, MPFR_RNDN),
                    mpfr_get_flt(mp_q, MPFR_RNDN)
                );
                fprintf(stderr, "%*s    expected q = %+li\n", indent, "", mp_q_int);
                fprintf(stderr, "%*s         got q = %+i\n", indent, "", q);
            }
            assert(q == mp_q_int);

            // mp_r = remainder(mp_x / mp_divisor)
            mpfr_t mp_r;
            mpfr_init2(mp_r, mp_precision);
            mpfr_fmod(mp_r, mp_x, mp_divisor, MPFR_RNDN);
            float mp_r_float = mpfr_get_flt(mp_r, MPFR_RNDN);
            abs_error = fabsf(r - mp_r_float);
            rel_error = abs_error / ulp(mp_r_float);
            if (rel_error >= 3.0f) {
                fprintf(stderr, "%*sx = %+.12f\n", indent, "", x);
                int N = significant_digits(mp_r_float);
                fprintf(
                    stderr,
                    "%*sremainder of %+.12f ÷ %+.12f = %+.*f\n",
                    indent, "",
                    mpfr_get_flt(mp_x, MPFR_RNDN),
                    mpfr_get_flt(mp_divisor, MPFR_RNDN),
                    N, mp_r_float
                );
                fprintf(stderr, "%*s    expected r = %+.*f\n", indent, "", N, mp_r_float);
                fprintf(stderr, "%*s         got r = %+.*f\n", indent, "", N, r);
                fprintf(stderr, "%*s         error = %+.*f\n", indent, "", N, abs_error);
                fprintf(stderr, "%*s        ulp(r) = %+.*f\n", indent, "", N, ulp(mp_r_float));
            }
            if (abs_error > max_abs_error)
                max_abs_error = abs_error;
            if (rel_error > max_rel_error)
                max_rel_error = rel_error;

            mpfr_clear(mp_r);
            mpfr_clear(mp_q);
            mpfr_clear(mp_x);

            if (rel_error >= 3.0f)
                error_dist[3]++;
            else if (rel_error >= 2.0f)
                error_dist[2]++;
            else if (rel_error >= 1.0f)
                error_dist[1]++;
            else
                error_dist[0]++;
        }
        max_abs_error_by_magnitude[j] = max_abs_error;
        max_rel_error_by_magnitude[j] = max_rel_error;
    }
    size_t m = n * sizeof magnitudes / sizeof magnitudes[0];
    assert(error_dist[0] + error_dist[1] + error_dist[2] + error_dist[3] == m);
    printf("\n");
    printf("%*serror distribution:\n", indent, "");
    printf("%*s 0 ulp %zu (%.2f%%)\n", indent, "", error_dist[0], (float) error_dist[0] / m * 100.0f);
    printf("%*s 1 ulp %zu (%.2f%%)\n", indent, "", error_dist[1], (float) error_dist[1] / m * 100.0f);
    printf("%*s 2 ulp %zu (%.2f%%)\n", indent, "", error_dist[2], (float) error_dist[2] / m * 100.0f);
    printf("%*s≥3 ulp %zu (%.2f%%)\n", indent, "", error_dist[3], (float) error_dist[3] / m * 100.0f);
    printf("\n");

    printf("%*smaximum error by magnitude:\n", indent, "");
    printf(
        "%*s%30s  %30s  %30s\n",
        indent, "",
        "order of magnitude of x",
        "absolute error",
        "relative error (ulp)"
    );
    printf("%*s┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈  ", indent, "");
    printf("┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈  ");
    printf("┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈\n");
    for (size_t j = 0; j < sizeof magnitudes / sizeof magnitudes[0]; j++) {
        printf(
            "%*s%30a  %30.12e  %30.1f\n",
            indent, "",
            magnitudes[j],
            max_abs_error_by_magnitude[j],
            max_rel_error_by_magnitude[j]
        );
    }
    printf("\n");
}

static void
test_q_and_r_double(rem_result_double (*rem_function)(double), mpfr_t mp_divisor)
{
    size_t error_dist[4] = {0};
    const double magnitudes[] = {
        0x1.0p+00,
        0x1.0p+01,
        0x1.0p+02,
        0x1.0p+04,
        0x1.0p+08,
        0x1.0p+16,
        0x1.0p+24,
        0x1.0p+32,
        0x1.0p+40,
        0x1.0p+48,
        0x1.0p+52,
    };
    double max_abs_error_by_magnitude[24] = {0.0};
    double max_rel_error_by_magnitude[24] = {0.0};
    for (size_t j = 0; j < sizeof magnitudes / sizeof magnitudes[0]; j++) {
        double max_abs_error = 0.0;
        double max_rel_error = 0.0;
        for (size_t i = 0; i < n; i++) {
            double x = (drand48() - 0.5) * 2.0 * magnitudes[j];
            mpfr_t mp_x;
            mpfr_init2(mp_x, mp_precision);
            mpfr_set_d(mp_x, x, MPFR_RNDN);

            rem_result_double rem = (*rem_function)(x);
            int64_t q = rem.z;
            double r = rem.v1 + rem.v2;
            double abs_error = 0.0;
            double rel_error = 0.0;

            // mp_q = trunc(mp_x / mp_divisor)
            mpfr_t mp_q;
            mpfr_init2(mp_q, mp_precision);
            mpfr_div(mp_q, mp_x, mp_divisor, MPFR_RNDN);
            mpfr_trunc(mp_q, mp_q);
            long int mp_q_int = mpfr_get_si(mp_q, MPFR_RNDN);
            if (q != mp_q_int) {
                fprintf(stderr, "%*sx = %.20f\n", indent, "", x);
                fprintf(
                    stderr,
                    "%*squotient of %.20f ÷ %.20f = %.0f\n",
                    indent, "",
                    mpfr_get_d(mp_x, MPFR_RNDN),
                    mpfr_get_d(mp_divisor, MPFR_RNDN),
                    mpfr_get_d(mp_q, MPFR_RNDN)
                );
                fprintf(stderr, "%*s    expected q = %+li\n", indent, "", mp_q_int);
                fprintf(stderr, "%*s         got q = %+li\n", indent, "", q);
            }
            assert(q == mp_q_int);

            // mp_r = remainder(mp_x / mp_divisor)
            mpfr_t mp_r;
            mpfr_init2(mp_r, mp_precision);
            mpfr_fmod(mp_r, mp_x, mp_divisor, MPFR_RNDN);
            double mp_r_double = mpfr_get_d(mp_r, MPFR_RNDN);
            abs_error = fabs(r - mp_r_double);
            rel_error = abs_error / ulp(mp_r_double);
            if (rel_error >= 3.0) {
                fprintf(stderr, "%*sx = %+.20f\n", indent, "", x);
                int N = significant_digits(mp_r_double);
                fprintf(
                    stderr,
                    "%*sremainder of %+.20f ÷ %+.20f = %+.*f\n",
                    indent, "",
                    mpfr_get_d(mp_x, MPFR_RNDN),
                    mpfr_get_d(mp_divisor, MPFR_RNDN),
                    N, mp_r_double
                );
                fprintf(stderr, "%*s    expected r = %+.*f\n", indent, "", N, mp_r_double);
                fprintf(stderr, "%*s         got r = %+.*f\n", indent, "", N, r);
                fprintf(stderr, "%*s         error = %+.*f\n", indent, "", N, abs_error);
                fprintf(stderr, "%*s        ulp(r) = %+.*f\n", indent, "", N, ulp(mp_r_double));
            }
            // assert(fabs(r - mp_r_double) <= ulp(mp_r_double));
            if (abs_error > max_abs_error)
                max_abs_error = abs_error;
            if (rel_error > max_rel_error)
                max_rel_error = rel_error;

            mpfr_clear(mp_r);
            mpfr_clear(mp_q);
            mpfr_clear(mp_x);

            if (rel_error >= 3.0)
                error_dist[3]++;
            else if (rel_error >= 2.0)
                error_dist[2]++;
            else if (rel_error >= 1.0)
                error_dist[1]++;
            else
                error_dist[0]++;
        }
        max_abs_error_by_magnitude[j] = max_abs_error;
        max_rel_error_by_magnitude[j] = max_rel_error;
    }
    size_t m = n * sizeof magnitudes / sizeof magnitudes[0];
    assert(error_dist[0] + error_dist[1] + error_dist[2] + error_dist[3] == m);
    printf("\n");
    printf("%*serror distribution:\n", indent, "");
    printf("%*s 0 ulp %zu (%.2f%%)\n", indent, "", error_dist[0], (double) error_dist[0] / m * 100.0);
    printf("%*s 1 ulp %zu (%.2f%%)\n", indent, "", error_dist[1], (double) error_dist[1] / m * 100.0);
    printf("%*s 2 ulp %zu (%.2f%%)\n", indent, "", error_dist[2], (double) error_dist[2] / m * 100.0);
    printf("%*s≥3 ulp %zu (%.2f%%)\n", indent, "", error_dist[3], (double) error_dist[3] / m * 100.0);
    printf("\n");

    printf("%*smaximum error by magnitude:\n", indent, "");
    printf(
        "%*s%30s  %30s  %30s\n",
        indent, "",
        "order of magnitude of x",
        "absolute error",
        "relative error (ulp)"
    );
    printf("%*s┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈  ", indent, "");
    printf("┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈  ");
    printf("┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈\n");
    for (size_t j = 0; j < sizeof magnitudes / sizeof magnitudes[0]; j++) {
        printf(
            "%*s%30a  %30.18e  %30.1f\n",
            indent, "",
            magnitudes[j],
            max_abs_error_by_magnitude[j],
            max_rel_error_by_magnitude[j]
        );
    }
    printf("\n");
}

int
main(int argc, char *argv[])
{
    // π
    mpfr_t mp_pi;
    mpfr_init2(mp_pi, mp_precision);
    mpfr_const_pi(mp_pi, MPFR_RNDN);
    {
        printf("%*stesting float [0, π] ...\n", indent, "");
        indent += 4;
        rem_result_float (*rem_pi)(float) = NULL;
        if (argc == 2 && strcmp(argv[1], "naive") == 0)
            rem_pi = &naive_rem_pi_float;
        else if (argc == 2 && strcmp(argv[1], "cw") == 0)
            rem_pi = &cw_rem_pi_float;
        else
            rem_pi = &bdl_rem_pi_float;
        assert(rem_pi != NULL);
        test_q_and_r_float(rem_pi, mp_pi);
        indent -= 4;
    }
    {
        printf("%*stesting double [0, π] ...\n", indent, "");
        indent += 4;
        rem_result_double (*rem_pi)(double) = NULL;
        if (argc == 2 && strcmp(argv[1], "naive") == 0)
            rem_pi = &naive_rem_pi_double;
        else if (argc == 2 && strcmp(argv[1], "cw") == 0)
            rem_pi = &cw_rem_pi_double;
        else
            rem_pi = &bdl_rem_pi_double;
        assert(rem_pi != NULL);
        test_q_and_r_double(rem_pi, mp_pi);
        indent -= 4;
    }

    // 2π
    mpfr_t mp_2pi;
    mpfr_init2(mp_2pi, mp_precision);
    mpfr_mul_d(mp_2pi, mp_pi, 2.0, MPFR_RNDN);
    {
        printf("%*stesting float [0, 2π] ...\n", indent, "");
        indent += 4;
        rem_result_float (*rem_2pi)(float) = NULL;
        if (argc == 2 && strcmp(argv[1], "naive") == 0)
            rem_2pi = &naive_rem_2pi_float;
        else if (argc == 2 && strcmp(argv[1], "cw") == 0)
            rem_2pi = &cw_rem_2pi_float;
        else
            rem_2pi = &bdl_rem_2pi_float;
        assert(rem_2pi != NULL);
        test_q_and_r_float(rem_2pi, mp_2pi);
        indent -= 4;
    }
    {
        printf("%*stesting double [0, 2π] ...\n", indent, "");
        indent += 4;
        rem_result_double (*rem_2pi)(double) = NULL;
        if (argc == 2 && strcmp(argv[1], "naive") == 0)
            rem_2pi = &naive_rem_2pi_double;
        else if (argc == 2 && strcmp(argv[1], "cw") == 0)
            rem_2pi = &cw_rem_2pi_double;
        else
            rem_2pi = &bdl_rem_2pi_double;
        assert(rem_2pi != NULL);
        test_q_and_r_double(rem_2pi, mp_2pi);
        indent -= 4;
    }

    // π∕2
    mpfr_t mp_pi_2;
    mpfr_init2(mp_pi_2, mp_precision);
    mpfr_div_d(mp_pi_2, mp_pi, 2.0, MPFR_RNDN);
    {
        printf("%*stesting float [0, π∕2] ...\n", indent, "");
        indent += 4;
        rem_result_float (*rem_pi_2)(float) = NULL;
        if (argc == 2 && strcmp(argv[1], "naive") == 0)
            rem_pi_2 = &naive_rem_pi_2_float;
        else if (argc == 2 && strcmp(argv[1], "cw") == 0)
            rem_pi_2 = &cw_rem_pi_2_float;
        else
            rem_pi_2 = &bdl_rem_pi_2_float;
        assert(rem_pi_2 != NULL);
        test_q_and_r_float(rem_pi_2, mp_pi_2);
        indent -= 4;
    }
    {
        printf("%*stesting double [0, π∕2] ...\n", indent, "");
        indent += 4;
        rem_result_double (*rem_pi_2)(double) = NULL;
        if (argc == 2 && strcmp(argv[1], "naive") == 0)
            rem_pi_2 = &naive_rem_pi_2_double;
        else if (argc == 2 && strcmp(argv[1], "cw") == 0)
            rem_pi_2 = &cw_rem_pi_2_double;
        else
            rem_pi_2 = &bdl_rem_pi_2_double;
        assert(rem_pi_2 != NULL);
        test_q_and_r_double(rem_pi_2, mp_pi_2);
        indent -= 4;
    }

    mpfr_clear(mp_pi_2);
    mpfr_clear(mp_2pi);
    mpfr_clear(mp_pi);

    return 0;
}
