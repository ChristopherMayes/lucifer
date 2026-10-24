/* fel_sincos_c: the (sin, cos) pair in one libm call. The FEL step's ODE needs both
   per particle per RK stage. gfortran emits two calls, and one paired call costs about
   as much as one of them.
   On macOS __sincos agrees bit for bit with the separate libm sin and cos. It differs
   from gfortran's sin and cos intrinsics by one ulp on about 2e-6 of arguments of the
   theta domain (doc/validation.md). */

#include <math.h>

#ifdef __APPLE__
extern void __sincos(double, double*, double*);
void fel_sincos_c(double th, double* s, double* c) { __sincos(th, s, c); }
#else
extern void sincos(double, double*, double*);
void fel_sincos_c(double th, double* s, double* c) { sincos(th, s, c); }
#endif
