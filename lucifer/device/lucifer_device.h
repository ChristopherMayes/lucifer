/* The device seam: the C interface fel_device_mod talks to through iso_c_binding.
 *
 * One implementation file sits behind this header per build: lucifer_metal.mm where
 * the toolchain can carry it (macOS with a Clang-family Objective-C++ compiler, since
 * the backend is ARC-managed Objective-C++ against the Metal framework) and
 * lucifer_device_stub.c everywhere else, which refuses and says why. No device type or
 * call appears outside that one file, so a second backend is a new implementation of
 * these functions plus one CMake branch. The shape follows GPUEngine.h of this project's
 * own prior work, the Genesis 1.3 v4 backends on branch gpu/metal-engine (commit
 * 4919b01, unmerged upstream): residency for the whole element, refusal with a stated reason,
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

/* The field set's size bound, wavefront_init%harmonics(9)'s. Mirrored by
 * fel_dev_max_field$ in fel_device_mod.f90. */
#define LUC_DEV_MAX_FIELD 9

/* One integration step's constants. Mirrored field for field by
 * fel_device_par_struct in fel_device_mod.f90: editing one side alone skews the
 * struct layout silently, the hazard MetalEngine.mm records for its own shader
 * structs. All doubles first, then 32-bit ints, so both mirrors agree on layout
 * without padding surprises. cret_ticks travels separately (luc_dev_step's own
 * argument) to keep this struct free of 64-bit members.
 *
 * The field set: nfield members, each a harmonic of the fundamental with its own
 * coupling fc(h) and deposit scale, and npol planes per member (Ex alone, or the
 * (Ex, Ey) pair when the run carries two live polarizations). The push gathers
 * every member at h*theta and the deposit writes every member's source with
 * exp(-i h theta). pol is the element's polarization 2-vector, read as conj(pol).E
 * in the gather and written as pol*src at the field add, fel_advance's and
 * fel_field_step's own conventions. With one member and one plane, every kernel
 * reduces to the single-field arithmetic it had before the set. */
typedef struct {
  double dz;             /* full step length [m] */
  double ks, ku;         /* fundamental radiation and undulator wavenumbers [1/m] */
  double aw;             /* rms undulator parameter */
  double qres;           /* 2 ku / ks, the detuning difference's resonance */
  double e0;             /* p0_mc^2 / gamma0: pz <-> goff conversion */
  double gam0;           /* reference gamma */
  double p0_mc;          /* reference momentum / m_e c */
  double kx, ky;         /* natural-focusing roll-off [1/m^2] */
  double ax, ay;         /* undulator field offset [m] */
  double cos_t, sin_t;   /* wiggle-plane tilt */
  double k1x, k1y;       /* Bmad transverse map k1 locals (pre 1/rel_p^2) */
  double gridmax, dgrid; /* grid half width and spacing [m], one grid for the set */
  double harm[LUC_DEV_MAX_FIELD];   /* per member: the harmonic number, as a real factor */
  double rtmp[LUC_DEV_MAX_FIELD];   /* per member: energy-exchange coupling fc(h) / (sqrt(2) m_e) */
  double scl_w[LUC_DEV_MAX_FIELD];  /* per member: deposit scale, fel_field_step's scl_w at h */
  double pol_re[2], pol_im[2];      /* the element's polarization pair on (Ex, Ey) */
  int32_t first;         /* field ring offset, Genesis's Field::first, one for the set */
  int32_t helical;       /* octupole kick shape (1 = both planes) */
  int32_t mutate;        /* falsifiability hook: perturb the kernel's detuning */
  int32_t nfield;        /* members in use, 1 to LUC_DEV_MAX_FIELD */
  int32_t npol;          /* planes per member, 1 or 2 */
  int32_t pad;
} luc_dev_step_par;

/* Backend presence and the device it would run on. Returns 1 when a backend is
 * compiled in and a usable device exists, else 0 with 'reason' naming what is
 * missing. Safe to call on any machine; allocates nothing. */
int luc_dev_available (char *name, int name_len, char *reason, int reason_len);

/* Allocate the resident buffers and compile the kernels for this run shape: nslice
 * slices of npart particles, and a field set of nfield members with npol planes
 * each on one ngrid grid. Refuses (nonzero return, reason filled): no device, a grid
 * the transform does not handle (powers of two 64 to 1024, the message names the
 * nearest supported size), a set outside 1..LUC_DEV_MAX_FIELD members or 1..2
 * planes, or a failed allocation with the wanted bytes as the device reports them. */
int luc_dev_init (int nslice, int npart, int ngrid, int nfield, int npol,
                  char *reason, int reason_len);

void luc_dev_close (void);

/* Grow the particle buffers to npart per slice, keeping the grid, the field set and
 * the compiled kernels: the kernels take npart at dispatch, so no recompile. Growth
 * only, since a smaller rectangle is never needed and shrinking would free buffers a
 * later element could ask for again. The buffers' contents are not preserved, and the
 * caller re-uploads every slice afterwards, which is what an element entry does. This
 * is what slice migration needs: fel_migrate_slices grows a slice's arrays on the
 * host, and the resident rectangle was sized at setup to the largest fill then. */
void luc_dev_resize_particles (int npart);

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

/* Field-plane transfers, 2*ngrid*ngrid floats interleaved re/im, addressed by
 * member im (0-based in the set), plane ip (0 = Ex, 1 = Ey) and record slice is,
 * plus the zero fill slippage needs and the source-grid readback the instrument's
 * rows read. The source is per member, not per plane: the polarization factors
 * apply at the field add, as in fel_field_step. */
void luc_dev_upload_field_slice (int im, int ip, int is, const float *e);
void luc_dev_download_field_slice (int im, int ip, int is, float *e);
void luc_dev_zero_field_slice (int im, int ip, int is);
void luc_dev_download_source_slice (int im, int is, float *s);

/* The step propagator exp(K2 dz) of member im, 2*ngrid*ngrid floats, FFT order:
 * each member diffracts at its own wavelength. The caller keys rebuilds
 * (fel_device_mod mirrors fp32_kernel_cache's key, per member). */
void luc_dev_set_kernel (int im, const float *expk);

/* Per-slice phase rotators e^{-i h (phi0 + ks z_ref)} per member, member-major
 * (index im*nslice + is), 2*nfield*nslice floats each: the push works against the
 * step's entry phase (base) and the deposit against its exit phase (base_dep, read
 * as i times the rotator), the two epochs fel_fp32_mod's twin uses. Uploaded per
 * step: phi0 advances, and in lockstep z_ref moves. */
void luc_dev_set_slice_phases (const float *base, const float *base_dep);

/* Encode one integration step: transverse half step, longitudinal push in the
 * (goff, phase-tick) chart gathering every member, transverse half step, source
 * deposit into every member, four-pass FFT field solve of every plane with its
 * member's propagator. One command buffer for the whole step; nothing is waited
 * on here. cret_ticks is this step's phi0 advance in ticks, subtracted exactly. */
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
