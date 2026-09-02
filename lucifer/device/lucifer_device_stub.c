/* The device seam with no backend behind it: every platform without one builds
 * this file, so the same tree compiles everywhere and a deck asking for a device
 * is refused by name rather than failing at link time or, worse, quietly taking
 * the CPU path. Only luc_dev_available can legitimately be called on a stub
 * build; everything else exists to satisfy the linker and would be a caller bug.
 */

#include "lucifer_device.h"

#include <stdio.h>
#include <string.h>

static void put_str (char *dst, int len, const char *src)
{
  if (dst == NULL || len < 1) return;
  snprintf (dst, (size_t) len, "%s", src);
}

int luc_dev_available (char *name, int name_len, char *reason, int reason_len)
{
  put_str (name, name_len, "none");
  put_str (reason, reason_len,
           "this build carries no device backend (the Metal backend builds on Darwin only)");
  return 0;
}

int luc_dev_init (int nslice, int npart, int ngrid, char *reason, int reason_len)
{
  (void) nslice; (void) npart; (void) ngrid;
  put_str (reason, reason_len, "no device backend in this build");
  return 1;
}

void luc_dev_close (void) {}

void luc_dev_upload_slice (int is, int n, const float *x, const float *px,
                           const float *y, const float *py, const float *goff,
                           const int64_t *uphase, const float *w)
{
  (void) is; (void) n; (void) x; (void) px; (void) y; (void) py;
  (void) goff; (void) uphase; (void) w;
}

void luc_dev_download_slice (int is, int n, float *x, float *px, float *y,
                             float *py, float *goff, int64_t *uphase)
{
  (void) is; (void) n; (void) x; (void) px; (void) y; (void) py;
  (void) goff; (void) uphase;
}

void luc_dev_upload_field_slice (int ifld, const float *e) { (void) ifld; (void) e; }
void luc_dev_download_field_slice (int ifld, float *e) { (void) ifld; (void) e; }
void luc_dev_zero_field_slice (int ifld) { (void) ifld; }
void luc_dev_download_source_slice (int ifld, float *s) { (void) ifld; (void) s; }
void luc_dev_set_kernel (const float *expk) { (void) expk; }
void luc_dev_set_slice_phases (const float *base, const float *base_dep)
{
  (void) base; (void) base_dep;
}

int luc_dev_step (const luc_dev_step_par *par, int64_t cret_ticks,
                  char *reason, int reason_len)
{
  (void) par; (void) cret_ticks;
  put_str (reason, reason_len, "no device backend in this build");
  return 1;
}

void luc_dev_sync (void) {}
int luc_dev_wrap_check (int64_t bucket_ticks) { (void) bucket_ticks; return 1; }
double luc_dev_seconds (void) { return 0; }
int64_t luc_dev_bytes (void) { return 0; }
