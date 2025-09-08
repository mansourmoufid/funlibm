CC=			clang
CXX=		clang++
LD=			$(CC)

CPPFLAGS+=	-I.
CPPFLAGS+=	-D_XOPEN_SOURCE=600 # drand48
ifeq ("$(DEBUG)","0")
CPPFLAGS+=	-DNDEBUG
endif

CFLAGS+=	-march=native
CFLAGS+=	-mtune=generic

# CFLAGS+=	-Wall -Wextra
CFLAGS+=	-Weverything
CFLAGS+=	-Wno-declaration-after-statement
CFLAGS+=	-pedantic

# debugging
CFLAGS+=	-fno-omit-frame-pointer
CFLAGS+=	-g

CFLAGS+=	-fno-fast-math # no cheating
CFLAGS+=	-fvectorize
CFLAGS+=	-Rpass=loop-vectorize
CFLAGS+=	-Rpass-analysis=loop-vectorize

CFLAGS+=	-std=c11

ifeq ("$(PROFILE)","1")
CFLAGS+=	-fcoverage-mapping
CFLAGS+=	-fprofile-instr-generate
LDFLAGS+=	-fprofile-instr-generate
endif

# GNU MPFR
ifeq ("$(shell command -v pkg-config)", "")
$(error "need pkg-config to find mpfr")
else
CPPFLAGS+=	$(shell pkg-config --cflags mpfr)
LDFLAGS+=	$(shell pkg-config --libs mpfr)
LDFLAGS+=	-rpath $(shell pkg-config --variable=prefix mpfr)
endif

LDFLAGS+=	-lm

%.s: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -O2 -S $< -o $@

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -O2 -c $< -o $@

reduce.h: arithmetic.h constants.h types.h
reduce.c: reduce.h types.h
test-reduce.c: common.h reduce.h types.h

test-reduce: test-reduce.c reduce.o
	$(CC) $(CPPFLAGS) -UNDEBUG $(CFLAGS) -O0 -c test-reduce.c -o test-reduce.o
	$(LD) reduce.o test-reduce.o -o test-reduce $(LDFLAGS)

.PHONY: clean
clean:
	rm -f *.s
	rm -f *.o
	rm -f test-reduce
