/* Why code/fel_sincos.c still calls libm.

   Replacing that call with arithmetic is the next thing anyone reaches for when
   speeding up the FEL step's particle path, and this is the audit that refuses it. The candidate here is the
   standard one, a Cody-Waite reduction with a two-term pi/2 and the fdlibm kernel
   polynomials, which is what a hand-written double-precision sincos looks like. The
   sweep reports how far it lands from libm over the argument range the step actually
   uses, and over ranges a longer line or a shorter wavelength would reach.

   The second question decides it. A two-term reduction degrades with |theta|, the
   caller has no way to notice, and theta is not bounded by anything the code owns: it
   carries the common phase accumulating along the line. Measured worst case runs 6 ulp
   to 1e4, 10 at 1e5, 24 at 1e6 and 2681 at 1e8.

   The shim this would replace was admitted at one ulp on 2e-6 of arguments. Anything
   proposing to replace it is measured against that number, not against the idea that
   arithmetic beats a library call. FINDINGS 7.38 records the decision, and
   doc/performance.md carries the table.

   Not built by the library. Build it when the question comes up again:

     cc -O2 -o sincos_audit sincos_audit.c -lm   */

#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

static const double TWO_OVER_PI = 6.36619772367581382433e-01;
static const double PIO2_1  = 1.57079632673412561417e+00;
static const double PIO2_1T = 6.07710050650619224932e-11;

static const double S1 = -1.66666666666666324348e-01, S2 =  8.33333333332248946124e-03;
static const double S3 = -1.98412698298579493134e-04, S4 =  2.75573137070700676789e-06;
static const double S5 = -2.50507477628578072866e-08, S6 =  1.58962301576546568060e-10;

static const double C1 =  4.16666666666666019037e-02, C2 = -1.38888888888741095749e-03;
static const double C3 =  2.48015872894767294178e-05, C4 = -2.75573143513906633035e-07;
static const double C5 =  2.08757232129817482790e-09, C6 = -1.13596475577881948265e-11;

static void poly_sincos(double th, double* s, double* c) {
  double fk = nearbyint(th * TWO_OVER_PI);
  double r  = (th - fk * PIO2_1) - fk * PIO2_1T;
  double z  = r * r;
  double sr = r + r * z * (S1 + z * (S2 + z * (S3 + z * (S4 + z * (S5 + z * S6)))));
  double cr = 1.0 - 0.5 * z + z * z * (C1 + z * (C2 + z * (C3 + z * (C4 + z * (C5 + z * C6)))));
  long k = ((long) fk) & 3L;
  switch (k) {
    case 0:  *s =  sr;  *c =  cr;  break;
    case 1:  *s =  cr;  *c = -sr;  break;
    case 2:  *s = -sr;  *c = -cr;  break;
    default: *s = -cr;  *c =  sr;  break;
  }
}

/* Distance in units in the last place, on the reference value's own exponent. */
static double ulp_diff(double got, double ref) {
  if (got == ref) return 0.0;
  double u = nextafter(fabs(ref), INFINITY) - fabs(ref);
  if (u == 0.0) u = 4.9406564584124654e-324;
  return fabs(got - ref) / u;
}

static void sweep(const char* label, double lo, double hi, long n) {
  double worst_s = 0, worst_c = 0, at_s = 0, at_c = 0;
  long ne_s = 0, ne_c = 0;
  for (long i = 0; i <= n; i++) {
    double th = lo + (hi - lo) * ((double) i / (double) n);
    double sp, cp, sr, cr;
    poly_sincos(th, &sp, &cp);
    sr = sin(th);  cr = cos(th);
    double us = ulp_diff(sp, sr), uc = ulp_diff(cp, cr);
    if (sp != sr) ne_s++;
    if (cp != cr) ne_c++;
    if (us > worst_s) { worst_s = us; at_s = th; }
    if (uc > worst_c) { worst_c = uc; at_c = th; }
  }
  printf("  %-26s  sin: worst %8.2f ulp at %+.4e, %6.2f%% differ | "
         "cos: worst %8.2f ulp at %+.4e, %6.2f%% differ\n",
         label, worst_s, at_s, 100.0 * ne_s / (n + 1),
         worst_c, at_c, 100.0 * ne_c / (n + 1));
}

int main(void) {
  printf("Cody-Waite (two term) + fdlibm kernels against libm sin/cos\n");
  printf("-----------------------------------------------------------------------------\n");
  sweep("measured range, averaged", -159.83, 19.03, 20000000L);
  sweep("measured range, unavg",    -127.41, 17.15, 20000000L);
  printf("\nHow it degrades as |theta| grows (the safety question):\n");
  sweep("|theta| < 1e2",  -1e2, 1e2, 5000000L);
  sweep("|theta| < 1e3",  -1e3, 1e3, 5000000L);
  sweep("|theta| < 1e4",  -1e4, 1e4, 5000000L);
  sweep("|theta| < 1e5",  -1e5, 1e5, 5000000L);
  sweep("|theta| < 1e6",  -1e6, 1e6, 5000000L);
  sweep("|theta| < 1e8",  -1e8, 1e8, 5000000L);
  return 0;
}
