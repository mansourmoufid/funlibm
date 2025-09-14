A fun experimental math library.


# Features

## Exact floating-point arithmetic

Implements exact addition[^1] and exact multiplication.[^2]

See the file named `arithmetic.h`.

## Argument reduction

Implements the Cody-Waite argument reduction algorithm.[^3]
See the files named `cw.h` and [`test-reduce-cw.txt`](test-reduce-cw.txt).

Implements the Boldo-Daumas-Li exact argument reduction algorithm.[^4]
See the files named `reduce.h` and [`test-reduce-bdl.txt`](test-reduce-bdl.txt).


# Requirements

 - [clang] or [GCC];
 - [GNU Make];
 - [GNU MPFR];
 - [pkg-config];

To install the requirements on Ubuntu:

    sudo apt install clang gcc libmpfr-dev make pkg-config


# Build

To build and test the argument reduction algorithm:

    make test-reduce
    ./test-reduce cw 2>/dev/null | tee test-reduce-cw.txt
    ./test-reduce bdl 2>/dev/null | tee test-reduce-bdl.txt


# References

[^1]: Marc Daumas, Laurence Rideau, Laurent Thery. A Generic Library for
    Floating-Point Numbers and Its Application to Exact Computing.
    Theorem Proving in Higher Order Logics, 2001, Edinburgh, United Kingdom.
    pp.169-184. https://hal.science/hal-00157285

[^2]: Alan H. Karp and Peter Markstein. 1997. High-precision division and
    square root. ACM Trans. Math. Softw. 23, 4 (Dec. 1997), 561–589.
    https://dl.acm.org/doi/pdf/10.1145/279232.279237

[^3]: W. J. Cody and W. Waite, Software manual for elementary functions.
    Prentice Hall, 1980.

[^4]: Sylvie Boldo, Marc Daumas, and Ren-Cang Li. "Formally verified argument
    reduction with a fused multiply-add."
    IEEE Transactions on Computers 58, no. 8 (2008): 1139-1145.
    https://arxiv.org/pdf/0708.3722


[clang]: https://clang.llvm.org/
[GCC]: https://gcc.gnu.org/
[GNU Make]: https://www.gnu.org/software/make/
[GNU MPFR]: https://www.mpfr.org/
[pkg-config]: https://www.freedesktop.org/wiki/Software/pkg-config/
