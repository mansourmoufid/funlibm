/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

/* Compile with -DNDEBUG, except for debugging. */

#pragma once

#include "constants.h" // pi, etc.
#include "polynomial.h" // eval_polynomial

static inline float _sin_pi_2_float(float x);

static inline double _sin_pi_2_double(double x);

#define _sin_0_pi_2(x) \
    _Generic((x), \
        float: _sin_0_pi_2_float, \
        double: _sin_0_pi_2_double \
    )(x)

static inline float
__attribute__((always_inline))
__attribute__((const))
_sin_0_pi_2_float(float x)
{
    assert(x >= 0.0f && x <= pi_2(float));
    const float as[22 + 1] = {
        /* x^0 */ -3.1513280585027375e-15f,
        /* x^1 */ 1.0f,
        /* x^2 */ 1.9950128354873087e-11f,
        /* x^3 */ -0.1666666716337204f,
        /* x^4 */ 2.4357765582294633e-08f,
        /* x^5 */ 0.008334207348525524f,
        /* x^6 */ -9.616296665626578e-06f,
        /* x^7 */ -0.00019375116971787065f,
        /* x^8 */ 0.0005109433550387621f,
        /* x^9 */ -0.004382471088320017f,
        /* x^10 */ 0.02023240737617016f,
        /* x^11 */ -0.06232820823788643f,
        /* x^12 */ 0.13840991258621216f,
        /* x^13 */ -0.2299337536096573f,
        /* x^14 */ 0.29131436347961426f,
        /* x^15 */ -0.28387677669525146f,
        /* x^16 */ 0.21279233694076538f,
        /* x^17 */ -0.12171296030282974f,
        /* x^18 */ 0.05217766389250755f,
        /* x^19 */ -0.01623234525322914f,
        /* x^20 */ 0.0034605867695063353f,
        /* x^21 */ -0.00045222308835946023f,
        /* x^22 */ 2.7316056730342098e-05f,
    };
    return eval_polynomial(as, sizeof as / sizeof as[0], x);
}

static inline double
__attribute__((always_inline))
__attribute__((const))
_sin_0_pi_2_double(double x)
{
    assert(x >= 0.0 && x <= pi_2(double));
    const double as[22 + 1] = {
        /* x^0 */ 6.043787009651245e-24,
        /* x^1 */ 1.0,
        /* x^2 */ -3.8476051657630597e-20,
        /* x^3 */ -0.16666666666666666,
        /* x^4 */ -4.560220434818188e-17,
        /* x^5 */ 0.008333333333331708,
        /* x^6 */ 1.7932328717674466e-14,
        /* x^7 */ -0.00019841269842175437,
        /* x^8 */ -9.496433393698673e-13,
        /* x^9 */ 2.75574008558096e-06,
        /* x^10 */ -3.7687566443477193e-11,
        /* x^11 */ -2.49359603476905e-08,
        /* x^12 */ -2.5800492027669335e-10,
        /* x^13 */ 5.893141042766951e-10,
        /* x^14 */ -5.433041056837119e-10,
        /* x^15 */ 5.287967272010351e-10,
        /* x^16 */ -3.970551101455e-10,
        /* x^17 */ 2.271697426732926e-10,
        /* x^18 */ -9.741207870774355e-11,
        /* x^19 */ 3.0313510689886855e-11,
        /* x^20 */ -6.464558935029198e-12,
        /* x^21 */ 8.45051197920576e-13,
        /* x^22 */ -5.10618043509081e-14,
    };
    return eval_polynomial(as, sizeof as / sizeof as[0], x);
}
