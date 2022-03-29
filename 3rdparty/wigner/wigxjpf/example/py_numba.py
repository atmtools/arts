from pywigxjpf_ffi import ffi, lib
import pywigxjpf_ffi
import time

import numba
from numba import cffi_support
cffi_support.register_module(pywigxjpf_ffi)
nb_wig3jj = pywigxjpf_ffi.lib.wig3jj

lib.wig_table_init(100,9);
lib.wig_temp_init(100);

print(lib.wig3jj(5,6,7,3,4,-7))
print(lib.wig6jj(6,6,6,6,6,6))
print(lib.wig9jj(6,6,6,7,7,6,7,7,8))
print(lib.wig9jj(6,6,6,7,7,6,7,7,8))
print(lib.wig6jj(6,6,6,6,6,int(4+0.5*4)))

# Benchmark 1, make calls for a lot of trivially-0 symbols.
@numba.jit(nopython=True)
def benchmark(jjmax):
    sum = 0.0
    calls = 0
    for jj1 in range(0, jjmax + 1):
        for jj2 in range(0, jjmax + 1):
            for jj3 in range(0, jjmax + 1):
                for mm1 in range(-jjmax, jjmax + 1):
                    for mm2 in range(-jjmax, jjmax + 1):
                        for mm3 in range(-jjmax, jjmax + 1):
                            w = nb_wig3jj(jj1, jj2, jj3, \
                                           mm1, mm2, mm3)
                            # print((w)
                            sum = sum + w
                            calls = calls+1
    return (sum,calls)

# Benchmark 2, avoiding trivial 0 by triangle rules.
@numba.jit(nopython=True)
def benchmark_opt(jjmax):
    sum = 0.0
    calls = 0
    for jj1 in range(0, jjmax + 1, 1):
        for jj2 in range(0, jjmax + 1, 1):
            jj3_min = abs(jj1-jj2)
            jj3_max = jj1+jj2
            if (jj3_max > jjmax):
                jj3_max = jjmax
            for jj3 in range(jj3_min, jj3_max + 1, 2):
                for mm1 in range(-jj1, jj1 + 1, 2):
                    for mm2 in range(-jj2, jj2 + 1, 2):
                        #for m3 in range(-j3, j3 + 1):
                        mm3 = -mm1-mm2
                        if (abs(mm3) <= jjmax):
                            w = nb_wig3jj(jj1, jj2, jj3, \
                                           mm1, mm2, mm3)
                            sum = sum + w
                            calls = calls+1
    return (sum,calls)

jjmax=10

for i in range(0,5):
    start_time = time.time()
    wigsum, total_calls = benchmark(jjmax)
    total_time = time.time()-start_time
    print("Benchmark 1 for jjmax=%d, sum=%.10f, time=%.5fs, "
          "time/call=%4.0fns [%d calls]" %
          (jjmax, wigsum, total_time, total_time/total_calls*1e9, total_calls))

for i in range(0,5):
    start_time = time.time()
    wigsum, total_calls = benchmark_opt(jjmax)
    total_time = time.time()-start_time
    print("Benchmark 2 for jjmax=%d, sum=%.10f, time=%.5fs, "
          "time/call=%4.0fns [%d calls]" %
          (jjmax, wigsum, total_time, total_time/total_calls*1e9, total_calls))

lib.wig_temp_free();
lib.wig_table_free();

print('Done')

# i7-8550U:   (with jjmax=20)

#                  incl. trivial-0    excl. trivial-0
#
#  Total:          638277381 calls        357599 calls
#
#  cffi   numba         16 ns/call        273 ns/call
#  C      native        10 ns/call        263 ns/call

# E3-1240 v3:  (with jjmax=10)

#                  incl. trivial-0    excl. trivial-0
#
#  Total:           12326391 calls        15475 calls
#
#  cffi   direct       188 ns/call        486 ns/call
#  cffi   numba         16 ns/call        185 ns/call
#  C      native        10 ns/call        175 ns/call
