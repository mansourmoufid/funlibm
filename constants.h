/* Copyright 2025, Mansour Moufid <mansourmoufid@gmail.com> */

#pragma once

// π
static const float pi_float = 3.1415927410125732421875f ;
static const float pi_hi_float = 3.141592502593994140625f ;
static const float pi_lo_float = 1.509958025280866422690451145172119140625e-7f ;
static const double pi_double = 3.141592653589793115997963468544185161590576171875 ;
static const double pi_hi_double = 3.141592653589793115997963468544185161590576171875 ;
static const double pi_lo_double = 1.2246467991473532071737640294583966046256921246776e-16 ;

#define pi(T) \
    _Generic((T){0}, \
        float: pi_float, \
        double: pi_double \
    )
#define pi_hi(T) \
    _Generic((T){0}, \
        float: pi_hi_float, \
        double: pi_hi_double \
    )
#define pi_lo(T) \
    _Generic((T){0}, \
        float: pi_lo_float, \
        double: pi_lo_double \
    )

// 2π
static const float twopi_float = 6.283185482025146484375f ;
static const float twopi_hi_float = 6.28318500518798828125f ;
static const float twopi_lo_float = 3.01991605056173284538090229034423828125e-7f ;
static const double twopi_double = 6.28318530717958623199592693708837032318115234375 ;
static const double twopi_hi_double = 6.28318530717958623199592693708837032318115234375 ;
static const double twopi_lo_double = 2.4492935982947064143475280589167932092513842493552e-16 ;

#define twopi(T) \
    _Generic((T){0}, \
        float: twopi_float, \
        double: twopi_double \
    )
#define twopi_hi(T) \
    _Generic((T){0}, \
        float: twopi_hi_float, \
        double: twopi_hi_double \
    )
#define twopi_lo(T) \
    _Generic((T){0}, \
        float: twopi_lo_float, \
        double: twopi_lo_double \
    )

// π∕2
static const float pi_2_float = 1.57079637050628662109375f ;
static const float pi_2_hi_float = 1.5707962512969970703125f ;
static const float pi_2_lo_float = 7.549790126404332113452255725860595703125e-8f ;
static const double pi_2_double = 1.5707963267948965579989817342720925807952880859375 ;
static const double pi_2_hi_double = 1.5707963267948965579989817342720925807952880859375 ;
static const double pi_2_lo_double = 6.123233995736766035868820147291983023128460623388e-17 ;

#define pi_2(T) \
    _Generic((T){0}, \
        float: pi_2_float, \
        double: pi_2_double \
    )
#define pi_2_hi(T) \
    _Generic((T){0}, \
        float: pi_2_hi_float, \
        double: pi_2_hi_double \
    )
#define pi_2_lo(T) \
    _Generic((T){0}, \
        float: pi_2_lo_float, \
        double: pi_2_lo_double \
    )

// π∕n
static const float pi_3_float = 1.0471975803375244140625f ;
static const float pi_4_float = 0.785398185253143310546875f ;
static const float pi_6_float = 0.52359879016876220703125f ;
static const float pi_8_float = 0.3926990926265716552734375f ;
static const float pi_16_float = 0.19634954631328582763671875f ;
static const double pi_3_double = 1.04719755119659785336239110620226711034774780273437 ;
static const double pi_4_double = 0.78539816339744827899949086713604629039764404296875 ;
static const double pi_6_double = 0.52359877559829892668119555310113355517387390136719 ;
static const double pi_8_double = 0.39269908169872413949974543356802314519882202148437 ;
static const double pi_16_double = 0.19634954084936206974987271678401157259941101074219 ;

// 1∕π
static const float inv_pi_float = 0.3183098733425140380859375f ;
static const float inv_pi_hi_float = 0.3183098733425140380859375f ;
static const float inv_pi_lo_float = 1.284127648659705300815403461456298828125e-8f ;
static const double inv_pi_double = 0.31830988618379069121644420192751567810773849487305 ;
static const double inv_pi_hi_double = 0.31830988618379063570529297066968865692615509033203 ;
static const double inv_pi_lo_double = 3.583247455607534113928639249044920144085258958814e-17 ;

// 1∕(2π)
static const float inv_2pi_float = 0.15915493667125701904296875f ;
static const float inv_2pi_hi_float = 0.15915493667125701904296875f ;
static const float inv_2pi_lo_float = 6.420638243298526504077017307281494140625e-9f ;
static const double inv_2pi_double = 0.159154943091895345608222100963757839053869247436523 ;
static const double inv_2pi_hi_double = 0.159154943091895317852646485334844328463077545166016 ;
static const double inv_2pi_lo_double = 1.791623727803767056964319624522460072042629479407e-17 ;

// 2∕π
static const float inv_pi_2_float = 0.636619746685028076171875f ;
static const float inv_pi_2_hi_float = 0.636619746685028076171875f ;
static const float inv_pi_2_lo_float = 2.56825529731941060163080692291259765625e-8f ;
static const double inv_pi_2_double = 0.6366197723675813824328884038550313562154769897461 ;
static const double inv_pi_2_hi_double = 0.63661977236758127141058594133937731385231018066406 ;
static const double inv_pi_2_lo_double = 7.166494911215068227857278498089840288170517917628e-17 ;

// k⋅2π ≤ 2²⁴, k = 2670176
static const float twopi_2_24_float = 16777211.0f ;
static const float twopi_2_24_hi_float = 16777210.0f ;
static const float twopi_2_24_lo_float = 0.61078357696533203125f ;
static const double twopi_2_24_double = 1.677721061078356020152568817138671875e7 ;
static const double twopi_2_24_hi_double = 1.67772106107835583388805389404296875e7 ;
static const double twopi_2_24_lo_double = 1.1617299155767537633817251742796711200256254414853e-9 ;
