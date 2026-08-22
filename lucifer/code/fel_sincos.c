/* fel_sincos_c: the (sin, cos) pair in ONE libm call. The FEL step's ODE evaluates
   both per particle per RK stage; gfortran emits two separate calls where one
   paired call costs about one. VALUE-PRESERVING, verified: on macOS __sincos is
   bitwise-identical to the separate sin/cos pair over a 44M-point sweep of the
   physical theta domain (the perf goal's lever-b test), so the benchmark keystone
   stays bit-for-bit. */

#include <math.h>

#ifdef __APPLE__
extern void __sincos(double, double*, double*);
void fel_sincos_c(double th, double* s, double* c) { __sincos(th, s, c); }
#else
extern void sincos(double, double*, double*);
void fel_sincos_c(double th, double* s, double* c) { sincos(th, s, c); }
#endif
