/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#include <fenv.h> // fegetround, FE_TONEAREST

#include "types.h" // rem_result_float, rem_result_double

#include "reduce.h"

/*
 * π
 */

rem_result_float
naive_rem_pi_float(float x)
{
    return _naive_rem_pi(x);
}

rem_result_double
naive_rem_pi_double(double x)
{
    return _naive_rem_pi(x);
}

rem_result_float
bdl_rem_pi_float(float x)
{
    return _bdl_rem_pi(x);
}

rem_result_double
bdl_rem_pi_double(double x)
{
    return _bdl_rem_pi(x);
}

/*
 * 2π
 */

rem_result_float
naive_rem_2pi_float(float x)
{
    return _naive_rem_2pi(x);
}

rem_result_double
naive_rem_2pi_double(double x)
{
    return _naive_rem_2pi(x);
}

rem_result_float
bdl_rem_2pi_float(float x)
{
    return _bdl_rem_2pi(x);
}

rem_result_double
bdl_rem_2pi_double(double x)
{
    return _bdl_rem_2pi(x);
}

/*
 * π∕2
 */

rem_result_float
naive_rem_pi_2_float(float x)
{
    return _naive_rem_pi_2(x);
}

rem_result_double
naive_rem_pi_2_double(double x)
{
    return _naive_rem_pi_2(x);
}

rem_result_float
bdl_rem_pi_2_float(float x)
{
    return _bdl_rem_pi_2(x);
}

rem_result_double
bdl_rem_pi_2_double(double x)
{
    return _bdl_rem_pi_2(x);
}

#undef NDEBUG
static void
__attribute__((constructor))
init(void)
{
    assert(fegetround() == FE_TONEAREST); // very important
}
