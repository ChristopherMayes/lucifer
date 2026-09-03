/* The device seam: the C interface fel_device_mod talks to through iso_c_binding.
 *
 * One implementation file sits behind this header per build: lucifer_metal.mm where
 * the toolchain can carry it (macOS with a Clang-family Objective-C++ compiler, since
 * the backend is ARC-managed Objective-C++ against the Metal framework) and
 * lucifer_device_stub.c everywhere else, which refuses by name. No device type or
 * call appears outside that one file, so a second backend is a new implementation of
 * these functions plus one CMake branch. The shape follows GPUEngine.h of this project's
 * own prior work, the Genesis 1.3 v4 backends on branch gpu/metal-engine (commit
 * 4919b01, unmerged upstream): residency for the whole element, refusal by name,
 * and single precision arranged for rather than accepted.
 *
 * Units and charts are the caller's business. Everything crossing this seam is
 * already in the device representation: transverse coordinates and the energy
 * offset as 32-bit floats, the longitudinal state as a 64-bit fixed-point phase
 * (ticks of 2 pi / 2^32 off the slice reference), weights as floats, field slices
 * as interleaved re/im float pairs. fel_device_mod owns every conversion, beside
 * the FP32 reformulations it mirrors from fel_fp32_mod.
 *
 * Slice and field indices at this seam are 0-based. Error reporting: functions
 * returning int give 0 on success and nonzero on refusal, with 'reason' filled as
 * a NUL-terminated string naming what was refused.
 */

#ifndef LUCIFER_DEVICE_H
#define LUCIFER_DEVICE_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* One integration step's constants. Mirrored field for field by
 * fel_device_par_struct in fel_device_mod.f90: editing one side alone skews the
 * struct layout silently, the hazard MetalEngine.mm records for its own shader
 * structs. All doubles first, then 32-bit ints, so both mirrors agree on layout
 * without padding surprises. cret_ticks travels separately (luc_dev_step's own
 * argument) to keep this struct free of 64-bit members. */
typedef struct {
  double dz;             /* full step length [m] */
  double ks, ku;         /* radiation and undulator wavenumbers [1/m] */
  double aw;             /* rms undulator parameter */
  double qres;           /* 2 ku / ks, the detuning difference's resonance */
  double e0;             /* p0_mc^2 / gamma0: pz <-> goff conversion */
  double gam0;           /* reference gamma */
  double p0_mc;          /* reference momentum / m_e c */
  double rtmp;           /* energy-exchange coupling fc / (sqrt(2) m_e) */
  double kx, ky;         /* natural-focusing roll-off [1/m^2] */
  double ax, ay;         /* undulator field offset [m] */
  double cos_t, sin_t;   /* wiggle-plane tilt */
  double k1x, k1y;       /* Bmad transverse map k1 locals (pre 1/rel_p^2) */
  double gridmax, dgrid; /* grid half width and spacing [m] */
  double scl_w;          /* deposit scale, fel_field_step's scl_w */
  int32_t first;         /* field ring offset, Genesis's Field::first */
  int32_t helical;       /* octupole kick shape (1 = both planes) */
  int32_t mutate;        /* falsifiability hook: perturb the kernel's detuning */
  int32_t pad;
} luc_dev_step_par;

/* Backend presence and the device it would run on. Returns 1 when a backend is
 * compiled in and a usable device exists, else 0 with 'reason' naming what is
 * missing. Safe to call on any machine; allocates nothing. */
int luc_dev_available (char *name, int name_len, char *reason, int reason_len);

/* Allocate the resident buffers and compile the kernels for this run shape.
 * Refuses by name (nonzero return, reason filled): no device, a grid the
 * transform does not handle (powers of two 64 to 1024, the message names the
 * nearest supported size), or a failed allocation with the wanted and free
 * bytes as the device reports them. */
int luc_dev_init (int nslice, int npart, int ngrid, char *reason, int reason_len);

void luc_dev_close (void);

/* Whole-slice transfers between the caller's staging arrays and the resident
 * buffers. The upload carries the weights; they are constant while resident
 * (migration is refused), so downloads do not return them.
 *
 * Concurrency contract: every transfer drains the device before touching a
 * buffer, and that drain is not thread-safe, so concurrent transfers are
 * illegal in general. After one serial luc_dev_sync with no further encoding,
 * however, a transfer is a pure copy of a caller-chosen region of the
 * shared-storage buffers, and transfers of DISJOINT regions (distinct slice
 * indices) may then run from concurrent threads. The parallel readback is
 * built on exactly that: one drain, then one download per slice per thread. */
void luc_dev_upload_slice (int is, int n, const float *x, const float *px,
                           const float *y, const float *py, const float *goff,
                           const int64_t *uphase, const float *w);
void luc_dev_download_slice (int is, int n, float *x, float *px, float *y,
                             float *py, float *goff, int64_t *uphase);

/* Field-slice transfers, 2*ngrid*ngrid floats interleaved re/im, plus the zero
 * fill slippage needs and the source-grid readback the instrument's rows read. */
void luc_dev_upload_field_slice (int ifld, const float *e);
void luc_dev_download_field_slice (int ifld, float *e);
void luc_dev_zero_field_slice (int ifld);
void luc_dev_download_source_slice (int ifld, float *s);

/* The step propagator exp(K2 dz), 2*ngrid*ngrid floats, FFT order. The caller
 * keys rebuilds (fel_device_mod mirrors fp32_kernel_cache's key). */
void luc_dev_set_kernel (const float *expk);

/* Per-slice phase rotators e^{-i(phi0 + ks z_ref)}, 2*nslice floats each: the
 * push works against the step's entry phase (base) and the deposit against its
 * exit phase (base_dep, read as i times the rotator), the two epochs
 * fel_fp32_mod's twin uses. Uploaded per step: phi0 advances, and in lockstep
 * z_ref moves. */
void luc_dev_set_slice_phases (const float *base, const float *base_dep);

/* Encode one integration step: transverse half step, longitudinal push in the
 * (goff, phase-tick) chart, transverse half step, source deposit, four-pass FFT
 * field solve. One command buffer for the whole step; nothing is waited on here.
 * cret_ticks is this step's phi0 advance in ticks, subtracted exactly. */
int luc_dev_step (const luc_dev_step_par *par, int64_t cret_ticks,
                  char *reason, int reason_len);

/* Drain the device: everything encoded completes before this returns. Every
 * transfer above syncs itself; this exists for the caller's own sequencing. */
void luc_dev_sync (void);

/* The exact-wrap assertion, run on the device itself: a probe set of phase
 * accumulators is shifted by whole buckets and back through device arithmetic.
 * Returns 0 only if the round trip is bit-exact and the extracted phase is
 * bit-identical under bucket shifts (wraps are modular arithmetic, so this
 * asserts exactly rather than to a tolerance). */
int luc_dev_wrap_check (int64_t bucket_ticks);

/* Seconds the device spent executing since init, from command-buffer
 * timestamps, and the resident footprint in bytes. */
double luc_dev_seconds (void);
int64_t luc_dev_bytes (void);

#ifdef __cplusplus
}
#endif

#endif
