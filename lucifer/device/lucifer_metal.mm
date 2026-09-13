// The Metal implementation of the device seam (lucifer_device.h): resident FP32
// state on Apple Silicon, one command buffer per integration step.
//
// Provenance. The transform kernels, the compare-and-swap-free device atomics,
// the command-buffer discipline and the buffer-sync rule are transcribed from
// this project's own prior work, the Genesis 1.3 v4 Metal backend:
// src/Core/MetalEngine.mm on branch gpu/metal-engine, commit 4919b01, unmerged
// upstream (cited by branch and commit since no release carries it). The physics
// kernels transcribe Lucifer's own sources instead: the push mirrors
// fel_fp32_mod's rk32/ode32 (the priced FP32 reformulations), the deposit
// mirrors its dep32, and the transverse map mirrors fel_transverse_track_bmad
// plus quad_mat2_calc (bmad/low_level/quad_mat2_calc.f90).
//
// One deliberate divergence from the reference backend, stated here where the
// citation is: the longitudinal state is not an absolute FP32 theta but a
// 64-bit fixed-point phase accumulator, in ticks of 2 pi / 2^32 off a static
// FP64 per-slice reference held by the host. The low 32 bits are the phase
// modulo one radiation period (uniform resolution 1.5e-9 rad, finer than FP32
// spacing anywhere past |theta| ~ 0.01), the high bits count whole periods, and
// a bucket crossing is exact integer arithmetic -- the migration operation a
// later landing needs, asserted exactly by luc_dev_wrap_check. The RK4 stages
// still run in FP32 radians on the small extracted angle, so the arithmetic
// floor stays the one the lockstep instrument priced.
//
// The field set. The resident field is nfield members (harmonics of the
// fundamental, each with its own coupling fc(h), deposit scale and propagator)
// of npol planes each (Ex, or the (Ex, Ey) pair when two polarizations are
// live), stored member-major: plane (m, p, is) sits at ((m*npol + p)*nslice +
// is). The source is per member, since the polarization factors apply at the
// field add exactly as fel_field_step writes pol*src. The push gathers every
// member once per particle and sums them into every RK stage at h*theta,
// fel_ode_multi's structure in fel_fp32_mod's reformulation. The deposit writes
// every member's source at exp(-i h theta), and the transform runs over every
// plane with its member's propagator. With one member and one plane every kernel is
// the single-field arithmetic it was before the set, term for term.
//
// The unaveraged mode. A record step there is nsub Strang substeps of half
// magnetic push, radiation kick with its deposit, half push, and the slice's own
// diffract and source add, so three kernels are that mode's own (unavg_push,
// unavg_kick and unavg_spont) and the rest of a substep is the kernels above:
// the four-pass solve, the accumulator's clear and the fixed-point scatter. The
// state rides the same seven particle buffers in that mode's chart, which the
// caller owns and lucifer_device.h states. What sets that mode's recorded levels
// is the transform pair rather than its own kernels: an FP32 pair loses about
// 1.3e-7 of the field's energy every time it runs, and this mode runs it nsub
// times a record step where the averaged mode runs it once (FINDINGS 7.72).
//
// The shader structs below are mirrored in host C++ immediately after the MSL
// string, and luc_dev_step_par is mirrored again in fel_device_mod.f90. Editing
// any one alone skews a buffer layout silently; keep all mirrors together.

#include "lucifer_device.h"

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

static const MTLResourceOptions kShared = MTLResourceStorageModeShared;

static NSString *kMSL = @R"MSL(
#include <metal_stdlib>
using namespace metal;

constant uint N = NG;

// The spontaneous reduction's shape: threads a group, and groups a slice.
#define SPONT_T 256u

// Ticks of phase: 2^32 per radiation period. The extraction takes the low 32
// bits as a SIGNED fraction, so the working angle sits in [-pi, pi) where FP32
// spacing is finest, and is exactly invariant under whole-period shifts.
constant float TICKS_PER_RAD = 683565275.576f;   // 2^32 / 2 pi
constant float RAD_PER_TICK  = 1.46291808e-9f;   // 2 pi / 2^32

inline float phase_of (long u){
    return float(as_type<int>(uint(u & 0xffffffffL))) * RAD_PER_TICK;
}

inline float2 cmul (float2 a, float2 b){ return float2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }

// The check's own failure hook: quantize to 256 ulps of the value, which is
// fel_fp32_mod's anint(x / (256 spacing(x))) * 256 spacing(x). The step has to come
// from the exponent alone. Taking it as a fixed fraction of |x| makes x/step the same
// integer for every x and the hook returns its argument unchanged, which is what this
// once did (FINDINGS 7.71).
inline float coarsen256 (float x){
    if (!(fabs(x) > 0.0f)) return x;
    int ex;
    frexp(x, ex);                          // |x| in [2^(ex-1), 2^ex)
    float sp = ldexp(1.0f, ex - 16);       // 256 ulps of x
    return rint(x / sp) * sp;
}

// ---------------- the source accumulator ----------------
// The deposit accumulates in fixed point so that the answer does not depend on the
// order threads reach a cell. Integer addition is associative and commutative where
// float addition is neither, so any arrival order gives one bit pattern.
//
// Sixty-four bits are carried as two 32-bit words because this hardware has no 64-bit
// atomic of any kind: atomic_ulong exists as a type and store, load, fetch_add,
// fetch_min, fetch_or and compare_exchange are all invalid for it. Every low addend is
// unsigned, so the number of carries is the number of times the low word passes 2^32,
// which is a property of the sum and not of the order, and both words are modular adds.
// The pair is therefore the exact 64-bit modular sum, and reinterpreting it as signed
// recovers the true sum whenever that sum fits, which the host's scale guarantees. An
// intermediate word may wrap freely: nothing reads one, and modular addition composes.
//
// Thirty-two bits would not do. A scale that keeps the worst case in an int leaves a
// quantum coarse enough to lose to the float accumulation it replaces, measured at 7.3
// times its error, where sixty-four bits beat it by the same factor.
//
// The scale is a power of two, so it and its reciprocal are exact and the conversion
// costs one rounding rather than the n the float accumulation paid.

inline void acc_fixed (device atomic_uint* SI, ulong ci, float v, float scale){
    ulong u = (ulong) (long) rint(v * scale);
    uint qlo = (uint) (u & 0xFFFFFFFFul);
    uint qhi = (uint) (u >> 32);
    uint prev = atomic_fetch_add_explicit(&SI[2u*ci], qlo, memory_order_relaxed);
    uint add_hi = qhi + ((prev + qlo < prev) ? 1u : 0u);
    if (add_hi != 0u) atomic_fetch_add_explicit(&SI[2u*ci + 1u], add_hi, memory_order_relaxed);
}

// One complex source element, decoded where its first consumer reads it.
inline float2 dec_fixed (const device uint* SI, ulong ie, float sinv){
    ulong b = 4u*ie;
    ulong ur = ((ulong) SI[b + 1u] << 32) | (ulong) SI[b];
    ulong ui = ((ulong) SI[b + 3u] << 32) | (ulong) SI[b + 2u];
    return float2(float((long) ur) * sinv, float((long) ui) * sinv);
}

// ---------------- transverse half step ----------------
// fel_transverse_track_bmad flattened per particle: optional tilt rotation in,
// the leading/trailing octupole-like kick, one quad_mat2_calc map per plane
// with the 1/rel_p^2 chromatic scaling, tilt rotation out. rel_p comes from
// the stored energy offset through the chart's exact relation pz = goff/e0.

struct TrkPar {
    float delz, inv_e0, k1x, k1y, kz, cos_t, sin_t;
    uint  npart, helical, leading;
};

inline void quad_map (float k1, float length, float rel_p, thread float* v, thread float* vp){
    // quad_mat2_calc: k1 > 0 defocuses. The small-argument branch keeps the
    // transcription complete; at FP32 it is unreachable for any real k1.
    float sqrt_k = sqrt(fabs(k1));
    float sk_l = sqrt_k * length;
    float cx, sx;
    if (fabs(sk_l) < 1e-10f) {
        float k_l2 = k1 * length * length;
        cx = 1.0f + k_l2 * 0.5f;
        sx = (1.0f + k_l2 / 6.0f) * length;
    } else if (k1 < 0.0f) {
        cx = cos(sk_l);
        sx = sin(sk_l) / sqrt_k;
    } else {
        cx = cosh(sk_l);
        sx = sinh(sk_l) / sqrt_k;
    }
    float v1 = *v, v2 = *vp;
    *v  = cx * v1 + (sx / rel_p) * v2;
    *vp = (k1 * sx * rel_p) * v1 + cx * v2;
}

kernel void trk_half (device float* X [[buffer(0)]], device float* PX [[buffer(1)]],
                      device float* Y [[buffer(2)]], device float* PY [[buffer(3)]],
                      const device float* G [[buffer(4)]],
                      constant TrkPar& P [[buffer(5)]],
                      uint gid [[thread_position_in_grid]]){
    float x = X[gid], px = PX[gid], y = Y[gid], py = PY[gid];
    float rel_p = 1.0f + G[gid] * P.inv_e0;
    float k1xx = P.k1x / (rel_p * rel_p);
    float k1yy = P.k1y / (rel_p * rel_p);
    float k3l = 2.0f * P.delz * k1yy;

    if (P.sin_t != 0.0f) {
        float t;
        t = P.cos_t * x  + P.sin_t * y;   y  = -P.sin_t * x  + P.cos_t * y;   x  = t;
        t = P.cos_t * px + P.sin_t * py;  py = -P.sin_t * px + P.cos_t * py;  px = t;
    }
    if (P.leading != 0u) {
        py += k3l * rel_p * P.kz * P.kz * y * y * y / 3.0f;
        if (P.helical != 0u) px += k3l * rel_p * P.kz * P.kz * x * x * x / 3.0f;
    }
    quad_map (k1xx, P.delz, rel_p, &x, &px);
    quad_map (k1yy, P.delz, rel_p, &y, &py);
    if (P.leading == 0u) {
        py += k3l * rel_p * P.kz * P.kz * y * y * y / 3.0f;
        if (P.helical != 0u) px += k3l * rel_p * P.kz * P.kz * x * x * x / 3.0f;
    }
    if (P.sin_t != 0.0f) {
        float t;
        t = P.cos_t * x  - P.sin_t * y;   y  =  P.sin_t * x  + P.cos_t * y;   x  = t;
        t = P.cos_t * px - P.sin_t * py;  py =  P.sin_t * px + P.cos_t * py;  px = t;
    }
    X[gid] = x; PX[gid] = px;
    Y[gid] = y; PY[gid] = py;
}

// ---------------- longitudinal push ----------------
// fel_fp32_mod's rk32/ode32 in the (goff, delta) chart: gamma as an offset,
// the detuning as the difference 0.5*ks*(qres - q), the phase factor through
// the per-slice FP64-seeded base rotator times the small-angle FP32 pair.
// The stage algebra is fel_runge_kutta's, verbatim. qmut is the check's own
// falsifiability hook: a perturbed kernel constant that must move the recorded
// theta level, or the ceilings prove nothing.

constant uint MAXF = 9;   // LUC_DEV_MAX_FIELD

struct PushPar {
    float dz, ks, qres, gam0, p0_mc, kx, ky, ax, ay, cos_t, sin_t;
    float gridmax, dgrid, aw, qmut;
    uint  ngrid, npart, nslice, first, nf, npol;
};

// Per member of the field set: the harmonic number as a factor, the coupling
// fc(h)/(sqrt(2) m_e) and the deposit scale. The polarization pair enters the
// gather as conj(pol) and the field add as pol.
struct FieldPar { float h, rtmp, scl, pad; };
struct SetPar {
    FieldPar f[MAXF];
    float2 polc[2];
    float2 pol[2];
};

struct Acc { float gg, pp; };

// fel_ode_multi's sum in the FP32 reformulation: every member's gathered phasor
// times its own rotator e^{-i h (phi0 + ks z_ref)} at its own phase h*d. One
// member with h = 1 is the single-field stage, term for term.
inline Acc ode (float g, float d, float btpar, thread const float2* rp,
                const device float2* BASE, uint is,
                constant PushPar& P, constant SetPar& S, Acc k){
    float2 ctmp = float2(0.0f);
    for (uint m = 0; m < P.nf; m++) {
        float dm = S.f[m].h * d;
        float s_t, c_t;
        s_t = sin(dm);  c_t = cos(dm);
        float2 rot = cmul(BASE[m * P.nslice + is], float2(c_t, -s_t));
        ctmp += cmul(rp[m], rot);
    }
    float gam_l = P.gam0 + g;
    float btper = btpar + (-2.0f / P.ks) * ctmp.x;
    float q_l = btper / (gam_l * gam_l);
    k.pp += 0.5f * P.ks * (P.qres * (1.0f + P.qmut) - q_l);
    k.gg += ctmp.y / gam_l;
    return k;
}

kernel void push (device float* G [[buffer(0)]], device long* U [[buffer(1)]],
                  const device float* X [[buffer(2)]], const device float* Y [[buffer(3)]],
                  const device float* PX [[buffer(4)]], const device float* PY [[buffer(5)]],
                  const device float2* F [[buffer(6)]],
                  const device float2* BASE [[buffer(7)]],
                  constant PushPar& P [[buffer(8)]],
                  constant long& cret [[buffer(9)]],
                  constant SetPar& S [[buffer(10)]],
                  uint gid [[thread_position_in_grid]]){
    uint is = gid / P.npart;
    float x = X[gid], y = Y[gid], px = PX[gid], py = PY[gid];

    // faw, first-order roll-off in the wiggle frame (fel_fp32_mod's faw32).
    float ddx = x - P.ax, ddy = y - P.ay;
    if (P.sin_t != 0.0f) {
        float t = P.cos_t * ddx + P.sin_t * ddy;
        ddy = -P.sin_t * ddx + P.cos_t * ddy;
        ddx = t;
    }
    float awloc = 1.0f + 0.5f * (P.kx * ddx * ddx + P.ky * ddy * ddy);
    float px_g = px * P.p0_mc, py_g = py * P.p0_mc;
    float btpar = 1.0f + px_g * px_g + py_g * py_g + P.aw * P.aw * awloc * awloc;

    // Bilinear gather from every member's FP32 planes, ring-rotated slice,
    // off-grid = dark (fel_fp32_mod's gather32, including the FP32 edge clamp).
    // Two live planes couple as conj(pol).E, fel_advance's own read. One plane
    // takes the scalar field untouched.
    float2 rp[MAXF];
    for (uint m = 0; m < MAXF; m++) rp[m] = float2(0.0f);
    if (x > -P.gridmax && x < P.gridmax && y > -P.gridmax && y < P.gridmax) {
        float wx = (x + P.gridmax) / P.dgrid;
        float wy = (y + P.gridmax) / P.dgrid;
        float fx = floor(wx), fy = floor(wy);
        wx = 1.0f + fx - wx;
        wy = 1.0f + fy - wy;
        int jx = int(fx), jy = int(fy);
        if (jx >= 0 && jy >= 0 && jx + 1 < int(P.ngrid) && jy + 1 < int(P.ngrid)) {
            uint fs = (is + P.first) % P.nslice;
            uint nn = P.ngrid * P.ngrid;
            uint cell = uint(jy) * P.ngrid + uint(jx);
            for (uint m = 0; m < P.nf; m++) {
                float2 cp[2] = { float2(0.0f), float2(0.0f) };
                for (uint q = 0; q < P.npol; q++) {
                    uint b = ((m * P.npol + q) * P.nslice + fs) * nn + cell;
                    cp[q] = F[b] * (wx * wy)
                          + F[b + 1u] * ((1.0f - wx) * wy)
                          + F[b + P.ngrid] * (wx * (1.0f - wy))
                          + F[b + P.ngrid + 1u] * ((1.0f - wx) * (1.0f - wy));
                }
                float2 c = (P.npol == 2u) ? (cmul(S.polc[0], cp[0]) + cmul(S.polc[1], cp[1]))
                                          : cp[0];
                float s = S.f[m].rtmp * awloc;
                rp[m] = float2(s * c.x, -s * c.y);
            }
        }
    }

    float g = G[gid];
    float d0 = phase_of(U[gid]);
    float d = d0;

    // fel_runge_kutta's stage bookkeeping, verbatim.
    Acc k2 = ode(g, d, btpar, rp, BASE, is, P, S, Acc{0.0f, 0.0f});
    float stpz = 0.5f * P.dz;
    g += stpz * k2.gg;  d += stpz * k2.pp;
    Acc k3 = k2;
    k2 = ode(g, d, btpar, rp, BASE, is, P, S, Acc{0.0f, 0.0f});
    g += stpz * (k2.gg - k3.gg);  d += stpz * (k2.pp - k3.pp);
    k3.gg /= 6.0f;  k3.pp /= 6.0f;
    k2.gg *= -0.5f; k2.pp *= -0.5f;
    k2 = ode(g, d, btpar, rp, BASE, is, P, S, k2);
    stpz = P.dz;
    g += stpz * k2.gg;  d += stpz * k2.pp;
    k3.gg -= k2.gg;  k3.pp -= k2.pp;
    k2.gg *= 2.0f;   k2.pp *= 2.0f;
    k2 = ode(g, d, btpar, rp, BASE, is, P, S, k2);
    g += stpz * (k3.gg + k2.gg / 6.0f);
    d += stpz * (k3.pp + k2.pp / 6.0f);

    // The step's phase advance, rounded once into ticks; the common phi0
    // advance is subtracted exactly (whole ticks, computed FP64 on the host).
    long dt = long(rint((d - d0) * TICKS_PER_RAD));
    U[gid] += dt - cret;
    G[gid] = g;
}

// ---------------- source deposit ----------------
// fel_fp32_mod's dep32: faw2 (no half -- Genesis's own roll-off, transcribed),
// part = sqrt(faw2)*scl*w/gamma, cpart = i e^{-i theta} * part through the
// per-slice base rotator, bilinear scatter. The accumulation is in fixed point
// through acc_fixed, so the order threads reach a cell does not reach the answer:
// both reference backends add floats there and two of their runs differ in the
// last bit or two of the source (manual/GPU.md records it). Every member of the
// set is written from the one particle at its own phase h*theta and scale,
// fel_field_step's harm*theta, into its own source plane.

struct DepPar {
    float gridmax, dgrid, kx, ky, ax, ay, cos_t, sin_t, gam0, dscale, gam_floor;
    uint  ngrid, npart, nslice, first, mutate, nf;
};

// The accumulator's clear, four words a thread. It is a memory fill and nothing else, so
// what it costs is the width of a store and the number of threads issued, not arithmetic.
// One word a thread left it the largest single pass on a many-slice window.
kernel void zero_src (device uint4* s [[buffer(0)]], uint gid [[thread_position_in_grid]]){
    s[gid] = uint4(0u);
}

kernel void deposit (device atomic_uint* S [[buffer(0)]],
                     const device float* X [[buffer(1)]], const device float* Y [[buffer(2)]],
                     const device float* G [[buffer(3)]], const device long* U [[buffer(4)]],
                     const device float* W [[buffer(5)]],
                     const device float2* BASE [[buffer(6)]],
                     constant DepPar& P [[buffer(7)]],
                     constant SetPar& SP [[buffer(8)]],
                     device atomic_uint* FAULT [[buffer(9)]],
                     uint gid [[thread_position_in_grid]]){
    uint is = gid / P.npart;
    float x = X[gid], y = Y[gid];
    if (!(x > -P.gridmax && x < P.gridmax && y > -P.gridmax && y < P.gridmax)) return;

    float wx = (x + P.gridmax) / P.dgrid;
    float wy = (y + P.gridmax) / P.dgrid;
    float fx = floor(wx), fy = floor(wy);
    wx = 1.0f + fx - wx;
    wy = 1.0f + fy - wy;
    int jx = int(fx), jy = int(fy);
    if (jx < 0 || jy < 0 || jx + 1 >= int(P.ngrid) || jy + 1 >= int(P.ngrid)) return;

    float ddx = x - P.ax, ddy = y - P.ay;
    if (P.sin_t != 0.0f) {
        float t = P.cos_t * ddx + P.sin_t * ddy;
        ddy = -P.sin_t * ddx + P.cos_t * ddy;
        ddx = t;
    }
    float gam = P.gam0 + G[gid];

    // The scale the host chose assumes no particle deposits below this gamma, a
    // contribution carrying w/gamma. One that does is not converted and not accumulated,
    // and the fault is recorded for the host to refuse the run at the next readback,
    // before anything derived from this step is written. The test is written so that a
    // gamma that is not a number fails it. Checking here rather than on the host closes
    // the interval between two readbacks, which is the only place the bound could have
    // been breached and used.
    if (!(gam >= P.gam_floor)) {
        atomic_fetch_or_explicit(FAULT, 1u, memory_order_relaxed);
        return;
    }

    float sq = sqrt(1.0f + P.kx * ddx * ddx + P.ky * ddy * ddy);

    float d = phase_of(U[gid]);
    // The check's mutation reaches the source row here: coarsen the residual
    // angle by eight mantissa bits, fel_fp32_mod's own hook transcribed.
    if (P.mutate != 0u) d = coarsen256(d);

    uint fs = (is + P.first) % P.nslice;
    uint nn = P.ngrid * P.ngrid;
    uint cell = uint(jy) * P.ngrid + uint(jx);
    for (uint m = 0; m < P.nf; m++) {
        // The member's own scale, in the single-field order of operations.
        float ppart = sq * SP.f[m].scl * W[gid] / gam;
        float dm = SP.f[m].h * d;
        float s_d, c_d;
        s_d = sin(dm);  c_d = cos(dm);

        // cbase = i * base, so (sin(phi+d) + i cos(phi+d)) = cbase * (cos d - i sin d),
        // with base the member's own rotator e^{-i h phi}.
        float2 b = BASE[m * P.nslice + is];
        float2 cp = cmul(float2(-b.y, b.x), float2(c_d, -s_d)) * ppart;

        ulong idx = (ulong) (m * P.nslice + fs) * nn + cell;
        float w;
        ulong dcell;
        w = wx * wy;                     dcell = 2u * idx;
        acc_fixed(S, dcell,      w * cp.x, P.dscale);
        acc_fixed(S, dcell + 1u, w * cp.y, P.dscale);
        w = (1.0f - wx) * wy;            dcell = 2u * (idx + 1u);
        acc_fixed(S, dcell,      w * cp.x, P.dscale);
        acc_fixed(S, dcell + 1u, w * cp.y, P.dscale);
        w = wx * (1.0f - wy);            dcell = 2u * (idx + P.ngrid);
        acc_fixed(S, dcell,      w * cp.x, P.dscale);
        acc_fixed(S, dcell + 1u, w * cp.y, P.dscale);
        w = (1.0f - wx) * (1.0f - wy);   dcell = 2u * (idx + P.ngrid + 1u);
        acc_fixed(S, dcell,      w * cp.x, P.dscale);
        acc_fixed(S, dcell + 1u, w * cp.y, P.dscale);
    }
}


// ---------------- the unaveraged mode ----------------
// The quiver-resolving advance (fel_unaveraged_mod's step header): nsub Strang
// substeps of half magnetic push, radiation kick and source deposit at the
// midpoint, half push, then the slice's own diffract and source add. The four
// passes of that last part are the transform kernels below, unchanged, so what
// is written here is the push, the kick with its deposit, and the ledger's
// spontaneous term. The FP64 routines are fel_unavg_bfield, unavg_ode,
// unavg_push_all and the kick block of fel_unavg_step, and the single-precision
// forms are fel_fp32_mod's fel_fp32_unavg_bfield, _ode and _push, whose
// divergence from FP64 doc/validation.md records.
//
// Three quantities do not depend on the particle and are computed in FP64 on the
// host, once a substep rather than once a particle, which is the hoist FINDINGS
// 7.37 records:
//
//   FQ    the envelope g, its slope gp, and cos(ku s), sin(ku s) at the four RK
//         stage positions of each half push. Eight float4 a substep.
//   CB    e^{i psi_mid} of each substep, per slice, the optical carrier's base.
//         psi_mid reaches some 1700 radians over a segment, which no float can
//         carry, so the host reduces it modulo 2 pi and sends the rotator.
//   the propagator exp(K2 dsub), uploaded through luc_dev_set_kernel as member 0.
//
// The lag. tau is the one state variable of the push that no other line reads:
// dtau/ds depends on gamma, ux and uy alone. The push therefore integrates it
// from zero over the substep and adds the increment to the 64-bit tick
// accumulator, so the lag is carried exactly and the arithmetic never differences
// the increment against a lag that grew. That is what the FP32 twin could not do,
// and its guard on the residual is the measurement of what it cost there
// (doc/validation.md): the residual's own quantum overtakes the per-substep change
// on a long window, where a tick is 2 pi / 2^32 of a radiation period whatever the
// window. The kick then reads the lag as phase_of, an angle in [-pi, pi) where
// FP32 spacing is finest, exactly as the averaged push reads its own.

struct UnavgPar {
    float h, dsub, ks, ku, aw, cos_t, sin_t;
    float gam0, beta0, g0inv2, me;
    float gridmax, dgrid, scl_u, dscale, u_bound;
    uint  ngrid, npart, nslice, first, helical, mutate;
};

struct U5 { float x, y, ux, uy, tau; };

inline U5 u5axpy (U5 a, float s, U5 b){
    U5 r;
    r.x = a.x + s*b.x;  r.y = a.y + s*b.y;  r.ux = a.ux + s*b.ux;
    r.uy = a.uy + s*b.uy;  r.tau = a.tau + s*b.tau;
    return r;
}

// fel_unavg_bfield, term for term, with the four s-dependent factors arriving in fq.
inline void ubfield (constant UnavgPar& P, float x, float y, float4 fq,
                     thread float& bx, thread float& by, thread float& bz){
    float xl = x, yl = y;
    if (P.sin_t != 0.0f) {
        xl =  P.cos_t * x + P.sin_t * y;
        yl = -P.sin_t * x + P.cos_t * y;
    }
    float g = fq.x, gp = fq.y, c_u = fq.z, s_u = fq.w;
    float a0;
    if (P.helical != 0u) {
        a0 = P.aw;
        float fperp = 1.0f + (P.ku*P.ku / 4.0f) * (xl*xl + yl*yl);
        bx = -a0 * (gp * s_u + g * P.ku * c_u) * fperp;
        by =  a0 * (gp * c_u - g * P.ku * s_u) * fperp;
        bz =  a0 * g * (P.ku*P.ku / 2.0f) * (s_u * xl - c_u * yl);
    } else {
        a0 = sqrt(2.0f) * P.aw;
        bx = 0.0f;
        by = a0 * (gp * c_u - g * P.ku * s_u) * cosh(P.ku * yl);
        bz = -a0 * g * c_u * P.ku * sinh(P.ku * yl);
    }
    if (P.sin_t != 0.0f) {
        float bt = P.cos_t * bx - P.sin_t * by;
        by = P.sin_t * bx + P.cos_t * by;
        bx = bt;
    }
}

// unavg_ode in single precision, fel_fp32_unavg_ode's reformulation of the fifth
// line: gamma/u_s - 1/beta0 is a difference of two numbers that are both one to
// within 1.3e-8, so in single precision it is exactly zero and the slippage is
// gone. The identity 1/ra - 1/rb = (a - b)/(ra rb (ra + rb)) moves the
// cancellation into a difference of two small like quantities (FINDINGS 7.69).
inline U5 uode (constant UnavgPar& P, U5 y, float4 fq, float goff){
    float gam = P.gam0 + goff;
    float bx, by, bz;
    ubfield(P, y.x, y.y, fq, bx, by, bz);
    float us = sqrt(gam*gam - 1.0f - y.ux*y.ux - y.uy*y.uy);
    U5 d;
    d.x = y.ux / us;
    d.y = y.uy / us;
    d.ux = by - y.uy * bz / us;
    d.uy = -bx + y.ux * bz / us;
    float aa = (1.0f + y.ux*y.ux + y.uy*y.uy) / (gam*gam);
    float ra = us / gam;
    d.tau = (aa - P.g0inv2) / (ra * P.beta0 * (ra + P.beta0));
    return d;
}

kernel void unavg_push (device float* X [[buffer(0)]], device float* Y [[buffer(1)]],
                        device float* UX [[buffer(2)]], device float* UY [[buffer(3)]],
                        device long* U [[buffer(4)]], const device float* G [[buffer(5)]],
                        const device float4* FQ [[buffer(6)]],
                        constant UnavgPar& P [[buffer(7)]],
                        uint gid [[thread_position_in_grid]]){
    U5 y0;
    y0.x = X[gid];  y0.y = Y[gid];  y0.ux = UX[gid];  y0.uy = UY[gid];  y0.tau = 0.0f;
    float goff = G[gid];
    float h = P.h;

    // unavg_push_all's stage bookkeeping, verbatim. gamma is untouched: B does no work.
    U5 k1 = uode(P, y0, FQ[0], goff);
    U5 k2 = uode(P, u5axpy(y0, 0.5f*h, k1), FQ[1], goff);
    U5 k3 = uode(P, u5axpy(y0, 0.5f*h, k2), FQ[2], goff);
    U5 k4 = uode(P, u5axpy(y0, h, k3), FQ[3], goff);

    float f = h / 6.0f;
    X[gid]  = y0.x  + f * (k1.x  + 2.0f*k2.x  + 2.0f*k3.x  + k4.x);
    Y[gid]  = y0.y  + f * (k1.y  + 2.0f*k2.y  + 2.0f*k3.y  + k4.y);
    UX[gid] = y0.ux + f * (k1.ux + 2.0f*k2.ux + 2.0f*k3.ux + k4.ux);
    UY[gid] = y0.uy + f * (k1.uy + 2.0f*k2.uy + 2.0f*k3.uy + k4.uy);

    // The lag's increment, rounded once into ticks and added exactly.
    float dtau = f * (k1.tau + 2.0f*k2.tau + 2.0f*k3.tau + k4.tau);
    U[gid] += long(rint((dtau * P.ks) * TICKS_PER_RAD));
}

// The radiation kick and the source deposit at the substep midpoint, one kernel
// because the two are exact energy duals only while they share operands: the same
// u_s, the same bilinear weights and the same phase (fel_unaveraged_mod's step
// header). u_s needs no reformulation, which measurement settled: it loses the 1
// and the ux^2 entirely, and they are 6.7e-9 of the result.
kernel void unavg_kick (device atomic_uint* S [[buffer(0)]],
                        device float* G [[buffer(1)]],
                        const device float* X [[buffer(2)]], const device float* Y [[buffer(3)]],
                        const device float* UX [[buffer(4)]], const device float* UY [[buffer(5)]],
                        const device float* W [[buffer(6)]], const device long* U [[buffer(7)]],
                        const device float2* E [[buffer(8)]],
                        const device float2* CB [[buffer(9)]],
                        constant UnavgPar& P [[buffer(10)]],
                        device atomic_uint* FAULT [[buffer(11)]],
                        uint gid [[thread_position_in_grid]]){
    uint is = gid / P.npart;
    float x = X[gid], y = Y[gid];

    // Off the grid is dark and deposits nothing, fel_unavg_step's on_grid branch.
    if (!(x > -P.gridmax && x < P.gridmax && y > -P.gridmax && y < P.gridmax)) return;
    float wx = (x + P.gridmax) / P.dgrid;
    float wy = (y + P.gridmax) / P.dgrid;
    float fx = floor(wx), fy = floor(wy);
    wx = 1.0f + fx - wx;
    wy = 1.0f + fy - wy;
    int jx = int(fx), jy = int(fy);
    if (jx < 0 || jy < 0 || jx + 1 >= int(P.ngrid) || jy + 1 >= int(P.ngrid)) return;

    float gam = P.gam0 + G[gid];
    float ux = UX[gid], uy = UY[gid];
    float us = sqrt(gam*gam - 1.0f - ux*ux - uy*uy);
    float2 jhat = (P.helical != 0u) ? float2(ux, -uy) * (1.0f / sqrt(2.0f))
                                    : float2(ux, 0.0f);

    // A contribution carries |j|/u_s, and the host's scale assumes that ratio stays
    // under u_bound. One that does not is neither converted nor accumulated, and the
    // fault stops the run at the next readback, before anything derived from this step
    // is written. Written so that a value that is not a number fails it. Checking here
    // rather than on the host closes the interval between two readbacks, which is the
    // only place the bound could have been breached and used.
    float ub = P.u_bound * us;
    if (!(jhat.x*jhat.x + jhat.y*jhat.y <= ub * ub)) {
        atomic_fetch_or_explicit(FAULT, 1u, memory_order_relaxed);
        return;
    }

    uint fs = (is + P.first) % P.nslice;
    uint nn = P.ngrid * P.ngrid;
    uint cell = uint(jy) * P.ngrid + uint(jx);
    uint b = fs * nn + cell;
    float2 ehat = E[b] * (wx * wy)
                + E[b + 1u] * ((1.0f - wx) * wy)
                + E[b + P.ngrid] * (wx * (1.0f - wy))
                + E[b + P.ngrid + 1u] * ((1.0f - wx) * (1.0f - wy));

    // The carrier e^{i(psi_mid - ks tau)}: the host's FP64 base rotator times the
    // lag's own small angle, which is the split the averaged push already uses.
    float d = phase_of(U[gid]);
    if (P.mutate != 0u) d = coarsen256(d);
    float s_d, c_d;
    s_d = sin(d);  c_d = cos(d);
    float2 cph = cmul(CB[is], float2(c_d, -s_d));

    // W = -i Ehat e^{i Psi};  dgamma = -dsub Re[W conj(j)] / (u_s m_e).
    float2 t = cmul(ehat, cph);
    float2 wph = float2(t.y, -t.x);
    float dgam = -P.dsub * (wph.x * jhat.x + wph.y * jhat.y) / (us * P.me);
    G[gid] = G[gid] + dgam;

    // src += i e^{-i Psi} j scl_u w / u_s. The /u_s where the averaged deposit has
    // Genesis's /gamma is what makes the pair exact duals (sec-unaveraged).
    float2 cj = cmul(float2(cph.x, -cph.y), jhat);
    float2 cdep = float2(-cj.y, cj.x) * (P.scl_u * W[gid] / us);

    ulong idx = (ulong) b;
    float w;
    ulong dcell;
    w = wx * wy;                     dcell = 2u * idx;
    acc_fixed(S, dcell,      w * cdep.x, P.dscale);
    acc_fixed(S, dcell + 1u, w * cdep.y, P.dscale);
    w = (1.0f - wx) * wy;            dcell = 2u * (idx + 1u);
    acc_fixed(S, dcell,      w * cdep.x, P.dscale);
    acc_fixed(S, dcell + 1u, w * cdep.y, P.dscale);
    w = wx * (1.0f - wy);            dcell = 2u * (idx + P.ngrid);
    acc_fixed(S, dcell,      w * cdep.x, P.dscale);
    acc_fixed(S, dcell + 1u, w * cdep.y, P.dscale);
    w = (1.0f - wx) * (1.0f - wy);   dcell = 2u * (idx + P.ngrid + 1u);
    acc_fixed(S, dcell,      w * cdep.x, P.dscale);
    acc_fixed(S, dcell + 1u, w * cdep.y, P.dscale);
}

// The ledger's spontaneous term, 4 sum|src|^2 over one slice's source grid, summed
// over the substeps of an element. This is the one field-energy increment the
// kick/deposit duality does not charge to the beam, so the time-dependent ledger
// closes only if it is banked (fel_unaveraged_mod's step).
//
// No atomic, and that is the point. One threadgroup owns one chunk of one slice
// and one accumulator slot, reduces its chunk in threadgroup memory over a fixed
// pairing, and adds the chunk's sum into its own slot. Nothing is contended, so
// the arithmetic runs in one order and two runs of a deck give one answer, which
// is what the fixed-point deposit bought for the source and what a float atomic
// would have given back here. The slot is a compensated pair, since a slice's
// grid contributes some five thousand times over an element and a bare float sum
// would carry that accumulation's own error into the ledger.
kernel void unavg_spont (const device uint* SI [[buffer(0)]],
                         device float2* ACC [[buffer(1)]],
                         constant uint& nn [[buffer(2)]],
                         constant float& sinv [[buffer(3)]],
                         constant float& fac [[buffer(4)]],
                         uint2 tg [[threadgroup_position_in_grid]],
                         uint t [[thread_index_in_threadgroup]],
                         uint2 ng [[threadgroups_per_grid]]){
    threadgroup float sh[SPONT_T];
    ulong base = (ulong) tg.y * (ulong) nn;
    float acc = 0.0f;
    for (uint i = tg.x * SPONT_T + t; i < nn; i += ng.x * SPONT_T) {
        float2 v = dec_fixed(SI, base + (ulong) i, sinv);
        acc += v.x * v.x + v.y * v.y;
    }
    sh[t] = acc;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint s = SPONT_T / 2u; s > 0u; s >>= 1u) {
        if (t < s) sh[t] += sh[t + s];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (t != 0u) return;
    uint slot = tg.y * ng.x + tg.x;
    float2 k = ACC[slot];
    float yv = sh[0] * fac - k.y;
    float tt = k.x + yv;
    k.y = (tt - k.x) - yv;
    k.x = tt;
    ACC[slot] = k;
}

// The record step's energy baseline: gamma changes only in the kick, so the step's
// beam-energy change is sum w (gamma_end - gamma_start) m_e, and the host forms it
// in FP64 from the two readbacks rather than reducing a bounded quantity here.
kernel void unavg_gsave (device float* G0 [[buffer(0)]], const device float* G [[buffer(1)]],
                         uint gid [[thread_position_in_grid]]){
    G0[gid] = G[gid];
}

// ---------------- the transform ----------------
// Four-step Cooley-Tukey N = REGS x LANES, transcribed whole from
// MetalEngine.mm (gpu/metal-engine, 4919b01): register DFTs with folded
// permutations, one XOR-swizzled threadgroup exchange, shapes injected as
// preprocessor macros. The solve is the same fused four passes, per plane:
//     field = IFFT(FFT(field) * expK_m) / N^2 + 2 * pol_p * src_m,
// with tg.y the plane, expK_m its member's propagator and src_m its member's
// source. One member and one plane is the scalar solve with pol = 1.

struct AddPar {
    uint   ppm, nslice, npol, pad;   // planes per member = npol * nslice
    float2 pol[2];
};

#define WSTRIDE(P) ((uint)(NG)/(uint)(P))

inline void dft4 (thread float2* a, float s){
    float2 t0 = a[0] + a[2], t1 = a[0] - a[2], t2 = a[1] + a[3], d = a[1] - a[3];
    float2 t3 = float2(s * d.y, -s * d.x);
    a[0] = t0 + t2; a[1] = t1 + t3; a[2] = t0 - t2; a[3] = t1 - t3;
}
inline void dft8 (thread float2* a, const device float2* W, float s){
    for (uint n2 = 0; n2 < 4; n2++){
        float2 u = a[n2], v = a[4 + n2];
        float2 w = W[(WSTRIDE(8) * n2) & (N - 1u)]; w.y *= s;
        a[n2]     = u + v;
        a[4 + n2] = cmul(u - v, w);
    }
    float2 t[4];
    for (uint k1 = 0; k1 < 2; k1++){
        for (uint n2 = 0; n2 < 4; n2++) t[n2] = a[4*k1 + n2];
        dft4(t, s);
        for (uint k2 = 0; k2 < 4; k2++) a[4*k1 + k2] = t[k2];
    }
}
inline void dft8n (thread float2* a, const device float2* W, float s){
    dft8(a, W, s);
    float2 t[8];
    for (uint j = 0; j < 8; j++) t[j] = a[j];
    for (uint j = 0; j < 8; j++) a[(j >> 2) + 2u*(j & 3u)] = t[j];
}
inline void dft16 (thread float2* a, const device float2* W, float s){
    float2 t[4];
    for (uint n2 = 0; n2 < 4; n2++){
        for (uint j = 0; j < 4; j++) t[j] = a[4*j + n2];
        dft4(t, s);
        for (uint k1 = 0; k1 < 4; k1++){
            float2 w = W[(WSTRIDE(16) * n2 * k1) & (N - 1u)]; w.y *= s;
            a[4*k1 + n2] = cmul(t[k1], w);
        }
    }
    for (uint k1 = 0; k1 < 4; k1++){
        for (uint n2 = 0; n2 < 4; n2++) t[n2] = a[4*k1 + n2];
        dft4(t, s);
        for (uint k2 = 0; k2 < 4; k2++) a[4*k1 + k2] = t[k2];
    }
}
inline void dft32 (thread float2* a, const device float2* W, float s){
    float2 t[8];
    for (uint n2 = 0; n2 < 8; n2++){
        for (uint j = 0; j < 4; j++) t[j] = a[8*j + n2];
        dft4(t, s);
        for (uint k1 = 0; k1 < 4; k1++){
            float2 w = W[(WSTRIDE(32) * n2 * k1) & (N - 1u)]; w.y *= s;
            a[8*k1 + n2] = cmul(t[k1], w);
        }
    }
    for (uint k1 = 0; k1 < 4; k1++){
        for (uint n2 = 0; n2 < 8; n2++) t[n2] = a[8*k1 + n2];
        dft8n(t, W, s);
        for (uint k2 = 0; k2 < 8; k2++) a[8*k1 + k2] = t[k2];
    }
}

#if REGS == 8
#  define DFT_REGS(a)  dft8(a, W, sgn)
#  define PERM_REGS(j) (((j) >> 2) + 2u*((j) & 3u))
#elif REGS == 16
#  define DFT_REGS(a)  dft16(a, W, sgn)
#  define PERM_REGS(j) (((j) >> 2) + 4u*((j) & 3u))
#elif REGS == 32
#  define DFT_REGS(a)  dft32(a, W, sgn)
#  define PERM_REGS(j) (((j) >> 3) + 4u*((j) & 7u))
#endif

#if LANES == 8
#  define DFT_LANES(a)  dft8(a, W, sgn)
#  define PERM_LANES(j) (((j) >> 2) + 2u*((j) & 3u))
#elif LANES == 16
#  define DFT_LANES(a)  dft16(a, W, sgn)
#  define PERM_LANES(j) (((j) >> 2) + 4u*((j) & 3u))
#elif LANES == 32
#  define DFT_LANES(a)  dft32(a, W, sgn)
#  define PERM_LANES(j) (((j) >> 3) + 4u*((j) & 7u))
#endif

inline void fftN (thread float2* a, threadgroup float2* s,
                  const device float2* W, uint lane, float sgn){
    DFT_REGS(a);
    for (uint j = 0; j < REGS; j++){
        uint k1 = PERM_REGS(j);
        float2 w = W[(lane * k1) & (N - 1u)]; w.y *= sgn;
        s[k1 * LANES + (lane ^ (k1 & (LANES - 1u)))] = cmul(a[j], w);
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint cc = 0; cc < CHUNK; cc++)
        for (uint n2 = 0; n2 < LANES; n2++)
            a[cc*LANES + n2] = s[(lane + cc*LANES) * LANES + (n2 ^ lane)];
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint cc = 0; cc < CHUNK; cc++) DFT_LANES(a + cc*LANES);
}
#define OUTK(cc, j) ((lane + (cc)*LANES) + REGS*PERM_LANES(j))

kernel void fft_rows (device float2* d [[buffer(0)]], const device float2* W [[buffer(1)]],
                      constant float& sgn [[buffer(2)]],
                      uint2 tg [[threadgroup_position_in_grid]], uint t [[thread_index_in_threadgroup]]){
    threadgroup float2 sh[RF_ROWS * NG];
    uint lane = t % LANES, r = t / LANES;
    device float2* p = d + (ulong)tg.y * ((ulong)N * N) + (ulong)(tg.x * RF_ROWS + r) * N;
    float2 a[REGS];
    for (uint n1 = 0; n1 < REGS; n1++) a[n1] = p[LANES*n1 + lane];
    fftN(a, sh + r*NG, W, lane, sgn);
    for (uint cc = 0; cc < CHUNK; cc++)
        for (uint j = 0; j < LANES; j++) p[OUTK(cc, j)] = a[cc*LANES + j];
}
// The unaveraged solve's first pass: the record with the substep's fixed-point source
// landed on it, E + 2 src, read, converted and transformed in one pass, so the source
// meets the diffraction on the record the kick read (fel_unaveraged_mod). The plane's
// source is found as fft_cols_add_fix finds it.
kernel void fft_rows_add_fix (device float2* d [[buffer(0)]], const device float2* W [[buffer(1)]],
                              constant float& sgn [[buffer(2)]], const device uint* si [[buffer(3)]],
                              constant AddPar& A [[buffer(4)]], constant float& sinv [[buffer(5)]],
                              uint2 tg [[threadgroup_position_in_grid]], uint t [[thread_index_in_threadgroup]]){
    threadgroup float2 sh[RF_ROWS * NG];
    uint lane = t % LANES, r = t / LANES;
    uint m = tg.y / A.ppm, rem = tg.y % A.ppm;
    uint pl = rem / A.nslice, is = rem % A.nslice;
    ulong row = (ulong)(tg.x * RF_ROWS + r);
    ulong off = (ulong)tg.y * ((ulong)N * N) + row * N;
    ulong soff = ((ulong)m * A.nslice + is) * ((ulong)N * N) + row * N;
    device float2* p = d + off;
    float2 a[REGS];
    if (A.npol == 2u) {
        float2 pol = A.pol[pl];
        for (uint n1 = 0; n1 < REGS; n1++){
            ulong q = LANES*n1 + lane;
            a[n1] = p[q] + 2.0f * cmul(pol, dec_fixed(si, soff + q, sinv));
        }
    } else {
        for (uint n1 = 0; n1 < REGS; n1++){
            ulong q = LANES*n1 + lane;
            a[n1] = p[q] + 2.0f * dec_fixed(si, soff + q, sinv);
        }
    }
    fftN(a, sh + r*NG, W, lane, sgn);
    for (uint cc = 0; cc < CHUNK; cc++)
        for (uint j = 0; j < LANES; j++) p[OUTK(cc, j)] = a[cc*LANES + j];
}
// The filter's first pass over the source: it reads the accumulator, converts, and
// transforms, so the conversion costs no dispatch of its own.
kernel void fft_rows_fix (device float2* d [[buffer(0)]], const device uint* si [[buffer(1)]],
                          const device float2* W [[buffer(2)]], constant float& sgn [[buffer(3)]],
                          constant float& sinv [[buffer(4)]],
                          uint2 tg [[threadgroup_position_in_grid]], uint t [[thread_index_in_threadgroup]]){
    threadgroup float2 sh[RF_ROWS * NG];
    uint lane = t % LANES, r = t / LANES;
    ulong row = (ulong)(tg.x * RF_ROWS + r);
    ulong base = (ulong)tg.y * ((ulong)N * N) + row * N;
    device float2* p = d + base;
    float2 a[REGS];
    for (uint n1 = 0; n1 < REGS; n1++) a[n1] = dec_fixed(si, base + LANES*n1 + lane, sinv);
    fftN(a, sh + r*NG, W, lane, sgn);
    for (uint cc = 0; cc < CHUNK; cc++)
        for (uint j = 0; j < LANES; j++) p[OUTK(cc, j)] = a[cc*LANES + j];
}
kernel void fft_rows_mul (device float2* d [[buffer(0)]], const device float2* W [[buffer(1)]],
                          constant float& sgn [[buffer(2)]], const device float2* expK [[buffer(3)]],
                          constant uint& ppm [[buffer(4)]],
                          uint2 tg [[threadgroup_position_in_grid]], uint t [[thread_index_in_threadgroup]]){
    threadgroup float2 sh[RF_ROWS * NG];
    uint lane = t % LANES, r = t / LANES;
    ulong row = (ulong)(tg.x * RF_ROWS + r);
    device float2* p = d + (ulong)tg.y * ((ulong)N * N) + row * N;
    const device float2* k = expK + (ulong)(tg.y / ppm) * ((ulong)N * N) + row * N;
    float2 a[REGS];
    for (uint n1 = 0; n1 < REGS; n1++) a[n1] = cmul(p[LANES*n1 + lane], k[LANES*n1 + lane]);
    fftN(a, sh + r*NG, W, lane, sgn);
    for (uint cc = 0; cc < CHUNK; cc++)
        for (uint j = 0; j < LANES; j++) p[OUTK(cc, j)] = a[cc*LANES + j];
}
kernel void fft_cols (device float2* d [[buffer(0)]], const device float2* W [[buffer(1)]],
                      constant float& sgn [[buffer(2)]], constant float& scale [[buffer(3)]],
                      uint2 tg [[threadgroup_position_in_grid]], uint t [[thread_index_in_threadgroup]]){
    threadgroup float2 sh[CC_COLS * NG];
    uint c = t % CC_COLS, lane = t / CC_COLS;
    device float2* b = d + (ulong)tg.y * ((ulong)N * N) + (ulong)tg.x * CC_COLS + c;
    float2 a[REGS];
    for (uint n1 = 0; n1 < REGS; n1++) a[n1] = b[(ulong)(LANES*n1 + lane) * N];
    fftN(a, sh + c*NG, W, lane, sgn);
    for (uint cc = 0; cc < CHUNK; cc++)
        for (uint j = 0; j < LANES; j++) b[(ulong)OUTK(cc, j) * N] = a[cc*LANES + j] * scale;
}
kernel void fft_cols_add (device float2* d [[buffer(0)]], const device float2* W [[buffer(1)]],
                          constant float& sgn [[buffer(2)]], constant float& scale [[buffer(3)]],
                          const device float2* src [[buffer(4)]],
                          constant AddPar& A [[buffer(5)]],
                          uint2 tg [[threadgroup_position_in_grid]], uint t [[thread_index_in_threadgroup]]){
    threadgroup float2 sh[CC_COLS * NG];
    uint c = t % CC_COLS, lane = t / CC_COLS;
    uint m = tg.y / A.ppm, rem = tg.y % A.ppm;
    uint pl = rem / A.nslice, is = rem % A.nslice;
    ulong off = (ulong)tg.y * ((ulong)N * N) + (ulong)tg.x * CC_COLS + c;
    ulong soff = ((ulong)m * A.nslice + is) * ((ulong)N * N) + (ulong)tg.x * CC_COLS + c;
    device float2* b = d + off;
    const device float2* sb = src + soff;
    float2 a[REGS];
    for (uint n1 = 0; n1 < REGS; n1++) a[n1] = b[(ulong)(LANES*n1 + lane) * N];
    fftN(a, sh + c*NG, W, lane, sgn);
    if (A.npol == 2u) {
        float2 pol = A.pol[pl];
        for (uint cc = 0; cc < CHUNK; cc++)
            for (uint j = 0; j < LANES; j++){
                ulong q = (ulong)OUTK(cc, j) * N;
                b[q] = a[cc*LANES + j] * scale + 2.0f * cmul(pol, sb[q]);
            }
    } else {
        for (uint cc = 0; cc < CHUNK; cc++)
            for (uint j = 0; j < LANES; j++){
                ulong q = (ulong)OUTK(cc, j) * N;
                b[q] = a[cc*LANES + j] * scale + 2.0f * sb[q];
            }
    }
}

// ---------------- the exact-wrap probes ----------------
// Device arithmetic itself, not a host emulation of it: the wrap check shifts
// resident accumulators by whole buckets and extracts phases on the GPU.

// The unfiltered solve's last pass: with no filter nothing transforms the source, so
// this is where it is first read and therefore where it is converted.
kernel void fft_cols_add_fix (device float2* d [[buffer(0)]], const device float2* W [[buffer(1)]],
                              constant float& sgn [[buffer(2)]], constant float& scale [[buffer(3)]],
                              const device uint* si [[buffer(4)]],
                              constant AddPar& A [[buffer(5)]], constant float& sinv [[buffer(6)]],
                              uint2 tg [[threadgroup_position_in_grid]], uint t [[thread_index_in_threadgroup]]){
    threadgroup float2 sh[CC_COLS * NG];
    uint c = t % CC_COLS, lane = t / CC_COLS;
    uint m = tg.y / A.ppm, rem = tg.y % A.ppm;
    uint pl = rem / A.nslice, is = rem % A.nslice;
    ulong off = (ulong)tg.y * ((ulong)N * N) + (ulong)tg.x * CC_COLS + c;
    ulong soff = ((ulong)m * A.nslice + is) * ((ulong)N * N) + (ulong)tg.x * CC_COLS + c;
    device float2* b = d + off;
    float2 a[REGS];
    for (uint n1 = 0; n1 < REGS; n1++) a[n1] = b[(ulong)(LANES*n1 + lane) * N];
    fftN(a, sh + c*NG, W, lane, sgn);
    if (A.npol == 2u) {
        float2 pol = A.pol[pl];
        for (uint cc = 0; cc < CHUNK; cc++)
            for (uint j = 0; j < LANES; j++){
                ulong q = (ulong)OUTK(cc, j) * N;
                b[q] = a[cc*LANES + j] * scale + 2.0f * cmul(pol, dec_fixed(si, soff + q, sinv));
            }
    } else {
        for (uint cc = 0; cc < CHUNK; cc++)
            for (uint j = 0; j < LANES; j++){
                ulong q = (ulong)OUTK(cc, j) * N;
                b[q] = a[cc*LANES + j] * scale + 2.0f * dec_fixed(si, soff + q, sinv);
            }
    }
}
kernel void wrap_shift (device long* u [[buffer(0)]], constant long& s [[buffer(1)]],
                        uint gid [[thread_position_in_grid]]){
    u[gid] += s;
}
kernel void wrap_phase (const device long* u [[buffer(0)]], device float* ph [[buffer(1)]],
                        uint gid [[thread_position_in_grid]]){
    ph[gid] = phase_of(u[gid]);
}
)MSL";

// Host mirrors of the MSL structs above. Editing one side alone skews the
// buffer layout silently; the reference backend records the same hazard.

struct TrkPar {
    float delz, inv_e0, k1x, k1y, kz, cos_t, sin_t;
    uint32_t npart, helical, leading;
};

struct PushPar {
    float dz, ks, qres, gam0, p0_mc, kx, ky, ax, ay, cos_t, sin_t;
    float gridmax, dgrid, aw, qmut;
    uint32_t ngrid, npart, nslice, first, nf, npol;
};

struct FieldPar { float h, rtmp, scl, pad; };
struct SetPar {
    FieldPar f[LUC_DEV_MAX_FIELD];
    float polc[4];    // float2 polc[2] in MSL: conj(pol), re/im interleaved
    float pol[4];     // float2 pol[2]
};

struct DepPar {
    float gridmax, dgrid, kx, ky, ax, ay, cos_t, sin_t, gam0, dscale, gam_floor;
    uint32_t ngrid, npart, nslice, first, mutate, nf;
};

struct AddPar {
    uint32_t ppm, nslice, npol, pad;
    float pol[4];     // float2 pol[2]
};

struct UnavgPar {
    float h, dsub, ks, ku, aw, cos_t, sin_t;
    float gam0, beta0, g0inv2, me;
    float gridmax, dgrid, scl_u, dscale, u_bound;
    uint32_t ngrid, npart, nslice, first, helical, mutate;
};

// The spontaneous reduction's shape, mirroring SPONT_T in the shader. One
// threadgroup owns one accumulator slot, so the count is part of the layout.
static const int kSpontThreads = 256;
static const int kSpontGroups = 64;
static const int kMaxSub = 4096;

namespace {

struct Impl {
    id<MTLDevice> dev {nil};
    id<MTLCommandQueue> queue {nil};

    int nslice {0}, npart {0}, ngrid {0}, nfield {1}, npol {1};
    int lanes {0}, regs {0}, rowsPerTG {0}, colsPerTG {0};

    // Buffer arithmetic for the member-major set. A plane is one slice of one
    // member's one polarization, and the source is one slice of one member.
    size_t plane (int im, int ip, int is) const {
        return ((size_t) im * npol + ip) * nslice + is;
    }
    size_t srcPlane (int im, int is) const { return (size_t) im * nslice + is; }
    size_t nplanes () const { return (size_t) nfield * npol * nslice; }

    id<MTLBuffer> bX, bPX, bY, bPY, bG, bU, bW;
    id<MTLBuffer> bField, bSrc, bSrcI, bExpK, bSig, bTw, bBase, bBaseDep;
    id<MTLBuffer> bProbe, bPh, bFault;
    // The unaveraged mode's own buffers: the record step's baseline energies, the
    // per-substep stage factors and carrier rotators, and the ledger's partials.
    id<MTLBuffer> bG0 {nil}, bFQ {nil}, bCB {nil}, bSpont {nil};
    int nsub {0};

    id<MTLComputePipelineState> pTrk, pPush, pZero, pDep, pRow, pRowM, pCol, pColA;
    id<MTLComputePipelineState> pRowF, pColAF;   // the two that convert where they read
    id<MTLComputePipelineState> pRowAF;          // the unaveraged solve's first pass
    id<MTLComputePipelineState> pWShift, pWPhase;
    id<MTLComputePipelineState> pUPush, pUKick, pUSpont, pUGsave;

    int64_t bytes {0};
    double busy {0};
    double srcScale {1};   // the last step's deposit scale, for the source readback

    // Per-pass timing (lucifer_device.h). Off by default: it costs an encoder a pass,
    // which the production path does not pay. The sample buffer takes two timestamps
    // an encoder and the hardware caps it at 4096 samples, so kMaxPass encoders fit in
    // one command buffer and a step that would pass that commits early.
    static const int kMaxPass = 2048;
    id<MTLCounterSampleBuffer> bSample {nil};
    bool timing {false};
    int nPass {0};
    int passSlot[kMaxPass];
    double passSec[LUC_DEV_PASS_N] {0};
    int64_t passCount[LUC_DEV_PASS_N] {0};
    int earlyCommits {0};

    // One command buffer per step, MetalEngine.mm's discipline: dispatches
    // accumulate into an open encoder, and every host touch of a buffer drains
    // first. sync() is the only committer, so no GPU work is in flight while
    // the host reads or writes the shared-storage buffers.
    id<MTLCommandBuffer> cb {nil};
    id<MTLComputeCommandEncoder> enc {nil};

    id<MTLComputeCommandEncoder> encoder (){
        if (enc == nil) {
            cb = [queue commandBuffer];
            enc = [cb computeCommandEncoder];
        }
        return enc;
    }

    // One pass's encoder. With timing off this is the shared encoder every dispatch
    // has always accumulated into. With it on the open encoder is closed and a new one
    // begun carrying this pass's two counter samples, which is what stage-boundary
    // sampling can measure. Every dispatch site sets its own pipeline state and all of
    // its buffers, so no encoder here inherits state from the one before it.
    id<MTLComputeCommandEncoder> pass (int slot){
        if (!timing) return encoder();
        if (enc != nil) { [enc endEncoding];  enc = nil; }
        if (nPass >= kMaxPass) { sync();  earlyCommits++; }
        if (cb == nil) cb = [queue commandBuffer];
        MTLComputePassDescriptor *pd = [MTLComputePassDescriptor computePassDescriptor];
        pd.dispatchType = MTLDispatchTypeSerial;
        pd.sampleBufferAttachments[0].sampleBuffer = bSample;
        pd.sampleBufferAttachments[0].startOfEncoderSampleIndex = (NSUInteger) (2 * nPass);
        pd.sampleBufferAttachments[0].endOfEncoderSampleIndex = (NSUInteger) (2 * nPass + 1);
        passSlot[nPass] = slot;
        nPass++;
        enc = [cb computeCommandEncoderWithDescriptor:pd];
        return enc;
    }

    // The timestamps of the passes in the command buffer just drained. A sample the
    // hardware could not take comes back with MTLCounterErrorValue, which is dropped
    // rather than accumulated as a wild interval.
    void resolve (){
        if (!timing || nPass == 0) { nPass = 0;  return; }
        NSData *rd = [bSample resolveCounterRange:NSMakeRange(0, (NSUInteger) (2 * nPass))];
        if (rd != nil) {
            const MTLCounterResultTimestamp *r = (const MTLCounterResultTimestamp *) [rd bytes];
            for (int i = 0; i < nPass; i++) {
                MTLTimestamp t0 = r[2*i].timestamp, t1 = r[2*i + 1].timestamp;
                if (t0 == MTLCounterErrorValue || t1 == MTLCounterErrorValue || t1 < t0) continue;
                passSec[passSlot[i]] += (double) (t1 - t0) * 1.0e-9;
                passCount[passSlot[i]] += 1;
            }
        }
        nPass = 0;
    }
    // Not thread-safe in general: encoder() and a draining sync() both touch enc.
    // What the header's concurrency contract rests on is the narrow case: after one
    // serial drain with nothing encoded since, enc is nil, concurrent transfer calls
    // all take the early return below, and their memcpys of disjoint regions of the
    // shared-storage buffers race nothing.
    void sync (){
        if (enc == nil && cb == nil) return;
        if (enc != nil) [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
        busy += [cb GPUEndTime] - [cb GPUStartTime];
        resolve();
        enc = nil;
        cb = nil;
    }
};

Impl *gImpl = nullptr;

void put_str (char *dst, int len, const std::string &src)
{
    if (dst == nullptr || len < 1) return;
    snprintf (dst, (size_t) len, "%s", src.c_str());
}

// The transform decomposes ngrid into two register stages, so the supported
// sizes are the reference backend's: powers of two, 64 to 1024, each with its
// blocking. Ties in the nearest-size message go to the larger grid, because
// dropping resolution silently is the worse surprise.
bool pickFFTShape (int ng, int &lanes, int &regs, int &rows, int &cols)
{
    switch (ng) {
    case   64: lanes =  8; regs =  8; rows = 16; cols = 16; return true;
    case  128: lanes =  8; regs = 16; rows = 16; cols = 16; return true;
    case  256: lanes = 16; regs = 16; rows =  8; cols = 16; return true;
    case  512: lanes = 16; regs = 32; rows =  4; cols =  8; return true;
    case 1024: lanes = 32; regs = 32; rows =  2; cols =  4; return true;
    default:   return false;
    }
}

int nearestSupported (int ng)
{
    static const int sizes[] = {64, 128, 256, 512, 1024};
    int best = sizes[0];
    double bd = 1e300;
    for (int i = 0; i < 5; i++) {
        const double d = std::fabs(std::log((double) sizes[i]) -
                                   std::log((double) (ng > 0 ? ng : 1)));
        if (d <= bd) { bd = d; best = sizes[i]; }
    }
    return best;
}

}   // namespace

extern "C" {

int luc_dev_available (char *name, int name_len, char *reason, int reason_len)
{
    @autoreleasepool {
        id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
        if (dev == nil) {
            put_str(name, name_len, "none");
            put_str(reason, reason_len, "no Metal device on this machine");
            return 0;
        }
        if (![dev hasUnifiedMemory]) {
            put_str(name, name_len, [[dev name] UTF8String]);
            put_str(reason, reason_len,
                    "the Metal device has no unified memory (a discrete GPU "
                    "would need the staging copies this design avoids)");
            return 0;
        }
        put_str(name, name_len, [[dev name] UTF8String]);
        put_str(reason, reason_len, "");
        return 1;
    }
}

int luc_dev_init (int nslice, int npart, int ngrid, int nfield, int npol,
                  char *reason, int reason_len)
{
    @autoreleasepool {
        luc_dev_close();
        if (nfield < 1 || nfield > LUC_DEV_MAX_FIELD) {
            char msg[128];
            snprintf(msg, sizeof(msg), "the field set has %d members; the device carries 1 to %d",
                     nfield, LUC_DEV_MAX_FIELD);
            put_str(reason, reason_len, msg);
            return 1;
        }
        if (npol < 1 || npol > 2) {
            char msg[128];
            snprintf(msg, sizeof(msg), "%d polarization planes; the device carries 1 or 2", npol);
            put_str(reason, reason_len, msg);
            return 1;
        }
        Impl *p = new Impl;

        p->dev = MTLCreateSystemDefaultDevice();
        if (p->dev == nil || ![p->dev hasUnifiedMemory]) {
            put_str(reason, reason_len, "no unified-memory Metal device");
            delete p;
            return 1;
        }
        if (!pickFFTShape(ngrid, p->lanes, p->regs, p->rowsPerTG, p->colsPerTG)) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "grid_n_pts = %d is not supported by the Metal field solver, "
                     "which handles powers of two from 64 to 1024; the nearest "
                     "supported size is %d", ngrid, nearestSupported(ngrid));
            put_str(reason, reason_len, msg);
            delete p;
            return 1;
        }
        p->queue = [p->dev newCommandQueue];
        p->nslice = nslice;
        p->npart = npart;
        p->ngrid = ngrid;
        p->nfield = nfield;
        p->npol = npol;

        const size_t np = (size_t) nslice * npart;
        const size_t nn = (size_t) ngrid * ngrid;
        auto alloc = [&](size_t n) -> id<MTLBuffer> {
            id<MTLBuffer> b = [p->dev newBufferWithLength:n options:kShared];
            p->bytes += (int64_t) n;
            return b;
        };
        p->bX = alloc(np * 4);  p->bPX = alloc(np * 4);
        p->bY = alloc(np * 4);  p->bPY = alloc(np * 4);
        p->bG = alloc(np * 4);  p->bU = alloc(np * 8);  p->bW = alloc(np * 4);
        p->bField = alloc(p->nplanes() * nn * 8);
        p->bSrc = alloc((size_t) nfield * nslice * nn * 8);
        p->bSrcI = alloc((size_t) nfield * nslice * nn * 16);
        p->bExpK = alloc((size_t) nfield * nn * 8);
        p->bSig = alloc((size_t) nfield * nn * 8);
        p->bTw = alloc((size_t) ngrid * 8);
        p->bBase = alloc((size_t) nfield * nslice * 8);
        p->bBaseDep = alloc((size_t) nfield * nslice * 8);
        p->bProbe = alloc(4096 * 8);
        p->bFault = alloc(4);
        if (p->bFault != nil) *((uint32_t *) [p->bFault contents]) = 0u;
        p->bPh = alloc(4096 * 4);
        if (p->bX == nil || p->bU == nil || p->bField == nil || p->bSrc == nil ||
            p->bSrcI == nil) {
            char msg[256];
            snprintf(msg, sizeof(msg),
                     "resident buffers do not fit: %lld MB wanted on %s",
                     (long long) (p->bytes / (1 << 20)),
                     [[p->dev name] UTF8String]);
            put_str(reason, reason_len, msg);
            delete p;
            return 1;
        }

        // The twiddle table, W_N = e^{-2 pi i m / N}, FP64 host trig rounded once.
        float *tw = (float *) [p->bTw contents];
        for (int m = 0; m < ngrid; m++) {
            const double a = -2.0 * M_PI * m / ngrid;
            tw[2*m]     = (float) cos(a);
            tw[2*m + 1] = (float) sin(a);
        }

        // Compile specialised to this grid, precise math: fast-math sincos and
        // free FP contraction would move the kernels off the priced arithmetic
        // for nothing a lockstep instrument could excuse.
        NSError *err = nil;
        MTLCompileOptions *copt = [[MTLCompileOptions alloc] init];
        copt.preprocessorMacros = @{
            @"NG"      : @(ngrid),
            @"LANES"   : @(p->lanes),
            @"REGS"    : @(p->regs),
            @"CHUNK"   : @(p->regs / p->lanes),
            @"RF_ROWS" : @(p->rowsPerTG),
            @"CC_COLS" : @(p->colsPerTG),
        };
        if (@available(macOS 15.0, *)) {
            copt.mathMode = MTLMathModeSafe;
        } else {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
            copt.fastMathEnabled = NO;
#pragma clang diagnostic pop
        }
        id<MTLLibrary> lib = [p->dev newLibraryWithSource:kMSL options:copt error:&err];
        if (lib == nil) {
            put_str(reason, reason_len,
                    std::string("shader compilation failed: ") +
                    [[err localizedDescription] UTF8String]);
            delete p;
            return 1;
        }
        bool ok = true;
        auto pso = [&](NSString *fname) -> id<MTLComputePipelineState> {
            NSError *e = nil;
            id<MTLFunction> fn = [lib newFunctionWithName:fname];
            id<MTLComputePipelineState> s =
                (fn == nil) ? nil : [p->dev newComputePipelineStateWithFunction:fn error:&e];
            if (s == nil) ok = false;
            return s;
        };
        p->pTrk = pso(@"trk_half");
        p->pPush = pso(@"push");
        p->pZero = pso(@"zero_src");
        p->pDep = pso(@"deposit");
        p->pRow = pso(@"fft_rows");
        p->pRowM = pso(@"fft_rows_mul");
        p->pCol = pso(@"fft_cols");
        p->pColA = pso(@"fft_cols_add");
        p->pRowF = pso(@"fft_rows_fix");
        p->pColAF = pso(@"fft_cols_add_fix");
        p->pRowAF = pso(@"fft_rows_add_fix");
        p->pWShift = pso(@"wrap_shift");
        p->pWPhase = pso(@"wrap_phase");
        p->pUPush = pso(@"unavg_push");
        p->pUKick = pso(@"unavg_kick");
        p->pUSpont = pso(@"unavg_spont");
        p->pUGsave = pso(@"unavg_gsave");
        if (!ok) {
            put_str(reason, reason_len, "compute pipeline creation failed");
            delete p;
            return 1;
        }
        gImpl = p;
        return 0;
    }
}

void luc_dev_close (void)
{
    if (gImpl == nullptr) return;
    gImpl->sync();
    delete gImpl;
    gImpl = nullptr;
}

void luc_dev_resize_particles (int npart)
{
    Impl *p = gImpl;
    if (p == nullptr || npart <= p->npart) return;
    @autoreleasepool {
        p->sync();
        const size_t np_old = (size_t) p->nslice * p->npart;
        const size_t np = (size_t) p->nslice * npart;
        // Seven particle buffers, six of 4 bytes and the phase accumulator of 8.
        p->bytes += (int64_t) ((np - np_old) * (6 * 4 + 8));
        p->bX = [p->dev newBufferWithLength:np * 4 options:kShared];
        p->bPX = [p->dev newBufferWithLength:np * 4 options:kShared];
        p->bY = [p->dev newBufferWithLength:np * 4 options:kShared];
        p->bPY = [p->dev newBufferWithLength:np * 4 options:kShared];
        p->bG = [p->dev newBufferWithLength:np * 4 options:kShared];
        p->bU = [p->dev newBufferWithLength:np * 8 options:kShared];
        p->bW = [p->dev newBufferWithLength:np * 4 options:kShared];
        if (p->bG0 != nil) p->bG0 = [p->dev newBufferWithLength:np * 4 options:kShared];
        p->npart = npart;
    }
}

void luc_dev_upload_slice (int is, int n, const float *x, const float *px,
                           const float *y, const float *py, const float *goff,
                           const int64_t *uphase, const float *w)
{
    Impl *p = gImpl;
    p->sync();
    const size_t o = (size_t) is * p->npart;
    const size_t nb = (size_t) n * 4;
    memcpy((float *) [p->bX contents] + o, x, nb);
    memcpy((float *) [p->bPX contents] + o, px, nb);
    memcpy((float *) [p->bY contents] + o, y, nb);
    memcpy((float *) [p->bPY contents] + o, py, nb);
    memcpy((float *) [p->bG contents] + o, goff, nb);
    memcpy((int64_t *) [p->bU contents] + o, uphase, (size_t) n * 8);
    memcpy((float *) [p->bW contents] + o, w, nb);
}

void luc_dev_download_slice (int is, int n, float *x, float *px, float *y,
                             float *py, float *goff, int64_t *uphase)
{
    Impl *p = gImpl;
    p->sync();
    const size_t o = (size_t) is * p->npart;
    const size_t nb = (size_t) n * 4;
    memcpy(x, (const float *) [p->bX contents] + o, nb);
    memcpy(px, (const float *) [p->bPX contents] + o, nb);
    memcpy(y, (const float *) [p->bY contents] + o, nb);
    memcpy(py, (const float *) [p->bPY contents] + o, nb);
    memcpy(goff, (const float *) [p->bG contents] + o, nb);
    memcpy(uphase, (const int64_t *) [p->bU contents] + o, (size_t) n * 8);
}

void luc_dev_upload_field_slice (int im, int ip, int is, const float *e)
{
    Impl *p = gImpl;
    p->sync();
    const size_t nn = (size_t) p->ngrid * p->ngrid;
    memcpy((float *) [p->bField contents] + p->plane(im, ip, is) * nn * 2, e, nn * 8);
}

void luc_dev_download_field_slice (int im, int ip, int is, float *e)
{
    Impl *p = gImpl;
    p->sync();
    const size_t nn = (size_t) p->ngrid * p->ngrid;
    memcpy(e, (const float *) [p->bField contents] + p->plane(im, ip, is) * nn * 2, nn * 8);
}

void luc_dev_zero_field_slice (int im, int ip, int is)
{
    Impl *p = gImpl;
    p->sync();
    const size_t nn = (size_t) p->ngrid * p->ngrid;
    memset((float *) [p->bField contents] + p->plane(im, ip, is) * nn * 2, 0, nn * 8);
}

void luc_dev_download_source_slice (int im, int is, float *s)
{
    Impl *p = gImpl;
    p->sync();
    const size_t nn = (size_t) p->ngrid * p->ngrid;

    // Decoded from the fixed-point accumulator rather than copied from the float source.
    // With the filter off nothing transforms the source and the solve's last pass
    // converts as it reads, so the float buffer holds no plane to copy. The decode here
    // is the same arithmetic that pass does, in double, so the instrument's rows read
    // what the device deposited. With the filter on that is now the deposit rather than
    // the filtered source this used to copy, which is what the caller compares against:
    // its reference is fel_fp32_deposit64, and no filter runs there.

    const uint32_t *w = (const uint32_t *) [p->bSrcI contents] + p->srcPlane(im, is) * nn * 4;
    const double inv = 1.0 / p->srcScale;
    for (size_t i = 0; i < nn; i++) {
        uint64_t ur = ((uint64_t) w[4*i + 1] << 32) | (uint64_t) w[4*i];
        uint64_t ui = ((uint64_t) w[4*i + 3] << 32) | (uint64_t) w[4*i + 2];
        s[2*i]     = (float) ((double) (int64_t) ur * inv);
        s[2*i + 1] = (float) ((double) (int64_t) ui * inv);
    }
}

int luc_dev_dep_fault (void)
{
    // Sticky, and read after a drain: any deposit that met a particle under the gamma the
    // host bounded the scale against set it, dropped that particle, and left the rest of
    // the step alone. The host refuses the run on seeing it.
    if (gImpl == nil || gImpl->bFault == nil) return 0;
    gImpl->sync();
    return (int) (*((const uint32_t *) [gImpl->bFault contents]) != 0u);
}

void luc_dev_set_kernel (int im, const float *expk)
{
    Impl *p = gImpl;
    p->sync();
    const size_t nn = (size_t) p->ngrid * p->ngrid;
    memcpy((float *) [p->bExpK contents] + (size_t) im * nn * 2, expk, nn * 8);
}

void luc_dev_set_filter (int im, const float *sig)
{
    Impl *p = gImpl;
    p->sync();
    const size_t nn = (size_t) p->ngrid * p->ngrid;
    memcpy((float *) [p->bSig contents] + (size_t) im * nn * 2, sig, nn * 8);
}

void luc_dev_set_slice_phases (const float *base, const float *base_dep)
{
    Impl *p = gImpl;
    p->sync();
    memcpy([p->bBase contents], base, (size_t) p->nfield * p->nslice * 8);
    memcpy([p->bBaseDep contents], base_dep, (size_t) p->nfield * p->nslice * 8);
}

int luc_dev_step (const luc_dev_step_par *par, int64_t cret_ticks,
                  char *reason, int reason_len)
{
    Impl *p = gImpl;
    if (p == nullptr) {
        put_str(reason, reason_len, "device not initialized");
        return 1;
    }
    if (par->nfield != p->nfield || par->npol != p->npol) {
        char msg[160];
        snprintf(msg, sizeof(msg),
                 "the step describes %d members of %d planes, and the device was initialized "
                 "for %d of %d", par->nfield, par->npol, p->nfield, p->npol);
        put_str(reason, reason_len, msg);
        return 1;
    }
    @autoreleasepool {
        const size_t nthread = (size_t) p->nslice * p->npart;
        const size_t nn = (size_t) p->ngrid * p->ngrid;
        const size_t nplane = p->nplanes();
        const MTLSize grid = MTLSizeMake(nthread, 1, 1);
        const MTLSize tgp = MTLSizeMake(256, 1, 1);
        const float fwd = 1.0f, inv = -1.0f, one = 1.0f;
        const float nrm = 1.0f / (float) nn;
        const long long cret = (long long) cret_ticks;
        const uint32_t ppm = (uint32_t) (p->npol * p->nslice);

        TrkPar T;
        T.delz = (float) (0.5 * par->dz);
        T.inv_e0 = (float) (1.0 / par->e0);
        T.k1x = (float) par->k1x;
        T.k1y = (float) par->k1y;
        T.kz = (float) par->ku;
        T.cos_t = (float) par->cos_t;
        T.sin_t = (float) par->sin_t;
        T.npart = (uint32_t) p->npart;
        T.helical = (uint32_t) par->helical;
        T.leading = 1;

        PushPar P;
        P.dz = (float) par->dz;
        P.ks = (float) par->ks;
        P.qres = (float) par->qres;
        P.gam0 = (float) par->gam0;
        P.p0_mc = (float) par->p0_mc;
        P.kx = (float) par->kx;
        P.ky = (float) par->ky;
        P.ax = (float) par->ax;
        P.ay = (float) par->ay;
        P.cos_t = (float) par->cos_t;
        P.sin_t = (float) par->sin_t;
        P.gridmax = (float) par->gridmax;
        P.dgrid = (float) par->dgrid;
        P.aw = (float) par->aw;
        // The falsifiability hook: 2^-12 on the detuning resonance, a wrong
        // kernel constant by construction, which must move the theta level.
        P.qmut = (par->mutate != 0) ? 2.44140625e-4f : 0.0f;
        P.ngrid = (uint32_t) p->ngrid;
        P.npart = (uint32_t) p->npart;
        P.nslice = (uint32_t) p->nslice;
        P.first = (uint32_t) par->first;
        P.nf = (uint32_t) par->nfield;
        P.npol = (uint32_t) par->npol;

        // The set: per-member coupling, scale and harmonic factor, and the
        // polarization pair as the gather reads it (conjugated) and the add
        // writes it. Unused members are zero and never dispatched over.
        SetPar S;
        memset(&S, 0, sizeof(S));
        for (int m = 0; m < par->nfield; m++) {
            S.f[m].h = (float) par->harm[m];
            S.f[m].rtmp = (float) par->rtmp[m];
            S.f[m].scl = (float) par->scl_w[m];
        }
        for (int q = 0; q < 2; q++) {
            S.pol[2*q] = (float) par->pol_re[q];
            S.pol[2*q + 1] = (float) par->pol_im[q];
            S.polc[2*q] = (float) par->pol_re[q];
            S.polc[2*q + 1] = (float) -par->pol_im[q];
        }

        // A power of two, so the scale and its reciprocal are both exact and the
        // conversion costs one rounding where the float accumulation paid n.
        const float dscale = (float) par->dep_scale;
        const float sinv = 1.0f / dscale;
        p->srcScale = par->dep_scale;

        DepPar D;
        D.dscale = dscale;
        D.gam_floor = (float) par->dep_gam_floor;
        D.gridmax = (float) par->gridmax;
        D.dgrid = (float) par->dgrid;
        D.kx = (float) par->kx;
        D.ky = (float) par->ky;
        D.ax = (float) par->ax;
        D.ay = (float) par->ay;
        D.cos_t = (float) par->cos_t;
        D.sin_t = (float) par->sin_t;
        D.gam0 = (float) par->gam0;
        D.ngrid = (uint32_t) p->ngrid;
        D.npart = (uint32_t) p->npart;
        D.nslice = (uint32_t) p->nslice;
        D.first = (uint32_t) par->first;
        D.mutate = (uint32_t) par->mutate;
        D.nf = (uint32_t) par->nfield;

        AddPar A;
        A.ppm = ppm;
        A.nslice = (uint32_t) p->nslice;
        A.npol = (uint32_t) p->npol;
        A.pad = 0;
        for (int i = 0; i < 4; i++) A.pol[i] = S.pol[i];

        id<MTLComputeCommandEncoder> e = nil;

        auto encTrk = [&](uint32_t leading) {
            T.leading = leading;
            e = p->pass(LUC_DEV_PASS_TRK);
            [e setComputePipelineState:p->pTrk];
            [e setBuffer:p->bX offset:0 atIndex:0];
            [e setBuffer:p->bPX offset:0 atIndex:1];
            [e setBuffer:p->bY offset:0 atIndex:2];
            [e setBuffer:p->bPY offset:0 atIndex:3];
            [e setBuffer:p->bG offset:0 atIndex:4];
            [e setBytes:&T length:sizeof(T) atIndex:5];
            [e dispatchThreads:grid threadsPerThreadgroup:tgp];
        };

        encTrk(1);

        e = p->pass(LUC_DEV_PASS_PUSH);
        [e setComputePipelineState:p->pPush];
        [e setBuffer:p->bG offset:0 atIndex:0];
        [e setBuffer:p->bU offset:0 atIndex:1];
        [e setBuffer:p->bX offset:0 atIndex:2];
        [e setBuffer:p->bY offset:0 atIndex:3];
        [e setBuffer:p->bPX offset:0 atIndex:4];
        [e setBuffer:p->bPY offset:0 atIndex:5];
        [e setBuffer:p->bField offset:0 atIndex:6];
        [e setBuffer:p->bBase offset:0 atIndex:7];
        [e setBytes:&P length:sizeof(P) atIndex:8];
        [e setBytes:&cret length:sizeof(cret) atIndex:9];
        [e setBytes:&S length:sizeof(S) atIndex:10];
        [e dispatchThreads:grid threadsPerThreadgroup:tgp];

        encTrk(0);

        e = p->pass(LUC_DEV_PASS_ZERO);
        [e setComputePipelineState:p->pZero];
        [e setBuffer:p->bSrcI offset:0 atIndex:0];
        // nn is a square of a power of two and every accumulator is four words, so the
        // fill divides by four exactly.
        [e dispatchThreads:MTLSizeMake((size_t) p->nfield * p->nslice * nn, 1, 1)
             threadsPerThreadgroup:tgp];

        e = p->pass(LUC_DEV_PASS_DEP);
        [e setComputePipelineState:p->pDep];
        [e setBuffer:p->bSrcI offset:0 atIndex:0];
        [e setBuffer:p->bX offset:0 atIndex:1];
        [e setBuffer:p->bY offset:0 atIndex:2];
        [e setBuffer:p->bG offset:0 atIndex:3];
        [e setBuffer:p->bU offset:0 atIndex:4];
        [e setBuffer:p->bW offset:0 atIndex:5];
        [e setBuffer:p->bBaseDep offset:0 atIndex:6];
        [e setBytes:&D length:sizeof(D) atIndex:7];
        [e setBytes:&S length:sizeof(S) atIndex:8];
        [e setBuffer:p->bFault offset:0 atIndex:9];
        [e dispatchThreads:grid threadsPerThreadgroup:tgp];

        // The source filter (fel-physics.md sec-source-filter). The CPU adds the filtered
        // source in Fourier space, inside the field's own transform pair, because it has one
        // to spare. The fused solve here adds the source in real space in its last pass, so
        // the source is filtered in place first and the solve is left exactly as it was:
        // IFFT(FFT(src) * sigmoid)/N^2 is the same quantity, reached by four passes over the
        // source buffer with the kernels the field already uses. The source has no
        // polarization plane, so its planes are member-major over slices alone.
        if (par->source_filter) {
            const uint32_t srcPPM = (uint32_t) p->nslice;
            const size_t nsrcplane = (size_t) par->nfield * (size_t) p->nslice;
            const MTLSize sRowTG = MTLSizeMake((size_t) (p->ngrid / p->rowsPerTG), nsrcplane, 1);
            const MTLSize sRowT = MTLSizeMake((size_t) (p->rowsPerTG * p->lanes), 1, 1);
            const MTLSize sColTG = MTLSizeMake((size_t) (p->ngrid / p->colsPerTG), nsrcplane, 1);
            const MTLSize sColT = MTLSizeMake((size_t) (p->colsPerTG * p->lanes), 1, 1);

            e = p->pass(LUC_DEV_PASS_FILTER);
            [e setComputePipelineState:p->pRowF];
            [e setBuffer:p->bSrc offset:0 atIndex:0];
            [e setBuffer:p->bSrcI offset:0 atIndex:1];
            [e setBuffer:p->bTw offset:0 atIndex:2];
            [e setBytes:&fwd length:4 atIndex:3];
            [e setBytes:&sinv length:4 atIndex:4];
            [e dispatchThreadgroups:sRowTG threadsPerThreadgroup:sRowT];

            e = p->pass(LUC_DEV_PASS_FILTER);
            [e setComputePipelineState:p->pCol];
            [e setBuffer:p->bSrc offset:0 atIndex:0];
            [e setBuffer:p->bTw offset:0 atIndex:1];
            [e setBytes:&fwd length:4 atIndex:2];
            [e setBytes:&one length:4 atIndex:3];
            [e dispatchThreadgroups:sColTG threadsPerThreadgroup:sColT];

            e = p->pass(LUC_DEV_PASS_FILTER);
            [e setComputePipelineState:p->pRowM];
            [e setBuffer:p->bSrc offset:0 atIndex:0];
            [e setBuffer:p->bTw offset:0 atIndex:1];
            [e setBytes:&inv length:4 atIndex:2];
            [e setBuffer:p->bSig offset:0 atIndex:3];
            [e setBytes:&srcPPM length:4 atIndex:4];
            [e dispatchThreadgroups:sRowTG threadsPerThreadgroup:sRowT];

            e = p->pass(LUC_DEV_PASS_FILTER);
            [e setComputePipelineState:p->pCol];
            [e setBuffer:p->bSrc offset:0 atIndex:0];
            [e setBuffer:p->bTw offset:0 atIndex:1];
            [e setBytes:&inv length:4 atIndex:2];
            [e setBytes:&nrm length:4 atIndex:3];
            [e dispatchThreadgroups:sColTG threadsPerThreadgroup:sColT];
        }

        // The transform over every plane of the set: tg.y indexes the plane.
        const MTLSize rowTG = MTLSizeMake((size_t) (p->ngrid / p->rowsPerTG), nplane, 1);
        const MTLSize rowT = MTLSizeMake((size_t) (p->rowsPerTG * p->lanes), 1, 1);
        const MTLSize colTG = MTLSizeMake((size_t) (p->ngrid / p->colsPerTG), nplane, 1);
        const MTLSize colT = MTLSizeMake((size_t) (p->colsPerTG * p->lanes), 1, 1);

        e = p->pass(LUC_DEV_PASS_SOLVE);
        [e setComputePipelineState:p->pRow];
        [e setBuffer:p->bField offset:0 atIndex:0];
        [e setBuffer:p->bTw offset:0 atIndex:1];
        [e setBytes:&fwd length:4 atIndex:2];
        [e dispatchThreadgroups:rowTG threadsPerThreadgroup:rowT];

        e = p->pass(LUC_DEV_PASS_SOLVE);
        [e setComputePipelineState:p->pCol];
        [e setBuffer:p->bField offset:0 atIndex:0];
        [e setBuffer:p->bTw offset:0 atIndex:1];
        [e setBytes:&fwd length:4 atIndex:2];
        [e setBytes:&one length:4 atIndex:3];
        [e dispatchThreadgroups:colTG threadsPerThreadgroup:colT];

        e = p->pass(LUC_DEV_PASS_SOLVE);
        [e setComputePipelineState:p->pRowM];
        [e setBuffer:p->bField offset:0 atIndex:0];
        [e setBuffer:p->bTw offset:0 atIndex:1];
        [e setBytes:&inv length:4 atIndex:2];
        [e setBuffer:p->bExpK offset:0 atIndex:3];
        [e setBytes:&ppm length:4 atIndex:4];
        [e dispatchThreadgroups:rowTG threadsPerThreadgroup:rowT];

        e = p->pass(LUC_DEV_PASS_SOLVE);
        [e setComputePipelineState:(par->source_filter ? p->pColA : p->pColAF)];
        [e setBuffer:p->bField offset:0 atIndex:0];
        [e setBuffer:p->bTw offset:0 atIndex:1];
        [e setBytes:&inv length:4 atIndex:2];
        [e setBytes:&nrm length:4 atIndex:3];
        [e setBuffer:(par->source_filter ? p->bSrc : p->bSrcI) offset:0 atIndex:4];
        [e setBytes:&A length:sizeof(A) atIndex:5];
        if (!par->source_filter) [e setBytes:&sinv length:4 atIndex:6];
        [e dispatchThreadgroups:colTG threadsPerThreadgroup:colT];
    }
    return 0;
}

int luc_dev_unavg_begin (int nsub, char *reason, int reason_len)
{
    Impl *p = gImpl;
    if (p == nullptr) {
        put_str(reason, reason_len, "device not initialized");
        return 1;
    }
    if (nsub < 1 || nsub > kMaxSub) {
        char msg[160];
        snprintf(msg, sizeof(msg), "the record step holds %d substeps; the device carries 1 to %d",
                 nsub, kMaxSub);
        put_str(reason, reason_len, msg);
        return 1;
    }
    @autoreleasepool {
        p->sync();
        const size_t np = (size_t) p->nslice * p->npart;
        if (p->bG0 == nil) {
            p->bG0 = [p->dev newBufferWithLength:np * 4 options:kShared];
            p->bytes += (int64_t) (np * 4);
        }
        if (p->bSpont == nil) {
            const size_t n = (size_t) p->nslice * kSpontGroups * 8;
            p->bSpont = [p->dev newBufferWithLength:n options:kShared];
            p->bytes += (int64_t) n;
        }
        if (nsub > p->nsub) {
            const size_t nfq = (size_t) nsub * 8 * 16;
            const size_t ncb = (size_t) nsub * p->nslice * 8;
            p->bFQ = [p->dev newBufferWithLength:nfq options:kShared];
            p->bCB = [p->dev newBufferWithLength:ncb options:kShared];
            p->bytes += (int64_t) (nfq + ncb);
            p->nsub = nsub;
        }
        if (p->bG0 == nil || p->bSpont == nil || p->bFQ == nil || p->bCB == nil) {
            put_str(reason, reason_len, "the unaveraged work buffers would not allocate");
            return 1;
        }
        memset([p->bSpont contents], 0, (size_t) p->nslice * kSpontGroups * 8);
    }
    return 0;
}

int luc_dev_unavg_step (const luc_dev_unavg_par *par, const float *fq,
                        const float *cbase, char *reason, int reason_len)
{
    Impl *p = gImpl;
    if (p == nullptr) {
        put_str(reason, reason_len, "device not initialized");
        return 1;
    }
    if (par->nsub < 1 || par->nsub > p->nsub) {
        char msg[160];
        snprintf(msg, sizeof(msg), "the record step holds %d substeps and the buffers carry %d",
                 par->nsub, p->nsub);
        put_str(reason, reason_len, msg);
        return 1;
    }
    @autoreleasepool {
        // The stage factors and the carrier rotators are host writes, so the device
        // drains before they land, as every other transfer here does.
        p->sync();
        memcpy([p->bFQ contents], fq, (size_t) par->nsub * 8 * 16);
        memcpy([p->bCB contents], cbase, (size_t) par->nsub * p->nslice * 8);

        const size_t nthread = (size_t) p->nslice * p->npart;
        const size_t nn = (size_t) p->ngrid * p->ngrid;
        const MTLSize grid = MTLSizeMake(nthread, 1, 1);
        const MTLSize tgp = MTLSizeMake(256, 1, 1);
        const float inv = -1.0f, fwd = 1.0f, one = 1.0f;
        const float nrm = 1.0f / (float) nn;
        const float dscale = (float) par->dep_scale;
        const float sinv = 1.0f / dscale;
        const float sfac = (float) par->spont_fac;
        const uint32_t unn = (uint32_t) nn;
        p->srcScale = par->dep_scale;

        UnavgPar P;
        P.h = (float) (0.5 * par->dsub);
        P.dsub = (float) par->dsub;
        P.ks = (float) par->ks;
        P.ku = (float) par->ku;
        P.aw = (float) par->aw;
        P.cos_t = (float) par->cos_t;
        P.sin_t = (float) par->sin_t;
        P.gam0 = (float) par->gam0;
        P.beta0 = (float) par->beta0;
        P.g0inv2 = (float) par->g0inv2;
        P.me = (float) par->m_electron;
        P.gridmax = (float) par->gridmax;
        P.dgrid = (float) par->dgrid;
        P.scl_u = (float) par->scl_u;
        P.dscale = dscale;
        P.u_bound = (float) par->dep_u_bound;
        P.ngrid = (uint32_t) p->ngrid;
        P.npart = (uint32_t) p->npart;
        P.nslice = (uint32_t) p->nslice;
        P.first = (uint32_t) par->first;
        P.helical = (uint32_t) par->helical;
        P.mutate = (uint32_t) par->mutate;

        AddPar A;
        A.ppm = (uint32_t) p->nslice;
        A.nslice = (uint32_t) p->nslice;
        A.npol = 1u;
        A.pad = 0;
        for (int i = 0; i < 4; i++) A.pol[i] = 0.0f;
        A.pol[0] = 1.0f;

        const size_t nplane = p->nplanes();
        const MTLSize rowTG = MTLSizeMake((size_t) (p->ngrid / p->rowsPerTG), nplane, 1);
        const MTLSize rowT = MTLSizeMake((size_t) (p->rowsPerTG * p->lanes), 1, 1);
        const MTLSize colTG = MTLSizeMake((size_t) (p->ngrid / p->colsPerTG), nplane, 1);
        const MTLSize colT = MTLSizeMake((size_t) (p->colsPerTG * p->lanes), 1, 1);
        const MTLSize spTG = MTLSizeMake((size_t) kSpontGroups, (size_t) p->nslice, 1);
        const MTLSize spT = MTLSizeMake((size_t) kSpontThreads, 1, 1);

        id<MTLComputeCommandEncoder> e = nil;

        // The record step's energy baseline. gamma changes only in the kick, so the
        // step's beam-energy change is a difference of two gammas and needs no
        // reduction of its own.
        e = p->pass(LUC_DEV_PASS_UPUSH);
        [e setComputePipelineState:p->pUGsave];
        [e setBuffer:p->bG0 offset:0 atIndex:0];
        [e setBuffer:p->bG offset:0 atIndex:1];
        [e dispatchThreads:grid threadsPerThreadgroup:tgp];

        for (int j = 0; j < par->nsub; j++) {
            e = p->pass(LUC_DEV_PASS_ZERO);
            [e setComputePipelineState:p->pZero];
            [e setBuffer:p->bSrcI offset:0 atIndex:0];
            [e dispatchThreads:MTLSizeMake((size_t) p->nslice * nn, 1, 1)
                 threadsPerThreadgroup:tgp];

            for (int half = 0; half < 2; half++) {
                if (half == 1) {
                    e = p->pass(LUC_DEV_PASS_UKICK);
                    [e setComputePipelineState:p->pUKick];
                    [e setBuffer:p->bSrcI offset:0 atIndex:0];
                    [e setBuffer:p->bG offset:0 atIndex:1];
                    [e setBuffer:p->bX offset:0 atIndex:2];
                    [e setBuffer:p->bY offset:0 atIndex:3];
                    [e setBuffer:p->bPX offset:0 atIndex:4];
                    [e setBuffer:p->bPY offset:0 atIndex:5];
                    [e setBuffer:p->bW offset:0 atIndex:6];
                    [e setBuffer:p->bU offset:0 atIndex:7];
                    [e setBuffer:p->bField offset:0 atIndex:8];
                    [e setBuffer:p->bCB offset:(NSUInteger) ((size_t) j * p->nslice * 8) atIndex:9];
                    [e setBytes:&P length:sizeof(P) atIndex:10];
                    [e setBuffer:p->bFault offset:0 atIndex:11];
                    [e dispatchThreads:grid threadsPerThreadgroup:tgp];
                }
                e = p->pass(LUC_DEV_PASS_UPUSH);
                [e setComputePipelineState:p->pUPush];
                [e setBuffer:p->bX offset:0 atIndex:0];
                [e setBuffer:p->bY offset:0 atIndex:1];
                [e setBuffer:p->bPX offset:0 atIndex:2];
                [e setBuffer:p->bPY offset:0 atIndex:3];
                [e setBuffer:p->bU offset:0 atIndex:4];
                [e setBuffer:p->bG offset:0 atIndex:5];
                [e setBuffer:p->bFQ offset:(NSUInteger) (((size_t) j * 8 + (size_t) half * 4) * 16)
                       atIndex:6];
                [e setBytes:&P length:sizeof(P) atIndex:7];
                [e dispatchThreads:grid threadsPerThreadgroup:tgp];
            }

            e = p->pass(LUC_DEV_PASS_USPONT);
            [e setComputePipelineState:p->pUSpont];
            [e setBuffer:p->bSrcI offset:0 atIndex:0];
            [e setBuffer:p->bSpont offset:0 atIndex:1];
            [e setBytes:&unn length:4 atIndex:2];
            [e setBytes:&sinv length:4 atIndex:3];
            [e setBytes:&sfac length:4 atIndex:4];
            [e dispatchThreadgroups:spTG threadsPerThreadgroup:spT];

            // The slice's own source add and diffract, four passes, the first landing the
            // substep's source on the record the kick read and the rest the transform
            // kernels unchanged: field = IFFT(FFT(field + 2 src) exp(K2 dsub))/N^2, the
            // order the CPU's step takes (fel_unaveraged_mod, FINDINGS 7.86).
            e = p->pass(LUC_DEV_PASS_SOLVE);
            [e setComputePipelineState:p->pRowAF];
            [e setBuffer:p->bField offset:0 atIndex:0];
            [e setBuffer:p->bTw offset:0 atIndex:1];
            [e setBytes:&fwd length:4 atIndex:2];
            [e setBuffer:p->bSrcI offset:0 atIndex:3];
            [e setBytes:&A length:sizeof(A) atIndex:4];
            [e setBytes:&sinv length:4 atIndex:5];
            [e dispatchThreadgroups:rowTG threadsPerThreadgroup:rowT];

            e = p->pass(LUC_DEV_PASS_SOLVE);
            [e setComputePipelineState:p->pCol];
            [e setBuffer:p->bField offset:0 atIndex:0];
            [e setBuffer:p->bTw offset:0 atIndex:1];
            [e setBytes:&fwd length:4 atIndex:2];
            [e setBytes:&one length:4 atIndex:3];
            [e dispatchThreadgroups:colTG threadsPerThreadgroup:colT];

            e = p->pass(LUC_DEV_PASS_SOLVE);
            [e setComputePipelineState:p->pRowM];
            [e setBuffer:p->bField offset:0 atIndex:0];
            [e setBuffer:p->bTw offset:0 atIndex:1];
            [e setBytes:&inv length:4 atIndex:2];
            [e setBuffer:p->bExpK offset:0 atIndex:3];
            [e setBytes:&A.ppm length:4 atIndex:4];
            [e dispatchThreadgroups:rowTG threadsPerThreadgroup:rowT];

            e = p->pass(LUC_DEV_PASS_SOLVE);
            [e setComputePipelineState:p->pCol];
            [e setBuffer:p->bField offset:0 atIndex:0];
            [e setBuffer:p->bTw offset:0 atIndex:1];
            [e setBytes:&inv length:4 atIndex:2];
            [e setBytes:&nrm length:4 atIndex:3];
            [e dispatchThreadgroups:colTG threadsPerThreadgroup:colT];
        }
    }
    return 0;
}

void luc_dev_download_g0 (int is, int n, float *g0)
{
    Impl *p = gImpl;
    if (p == nullptr || p->bG0 == nil) return;
    p->sync();
    memcpy(g0, (const float *) [p->bG0 contents] + (size_t) is * p->npart, (size_t) n * 4);
}

double luc_dev_unavg_spont (void)
{
    Impl *p = gImpl;
    if (p == nullptr || p->bSpont == nil) return 0;
    p->sync();
    // FP64 over the partials in index order, so the sum is the same whatever the
    // device did: every partial was written by one threadgroup and by no other.
    const float *a = (const float *) [p->bSpont contents];
    double tot = 0;
    for (int i = 0; i < p->nslice * kSpontGroups; i++) tot += (double) a[2*i];
    return tot;
}

void luc_dev_sync (void)
{
    if (gImpl != nullptr) gImpl->sync();
}

int luc_dev_wrap_check (int64_t bucket_ticks)
{
    Impl *p = gImpl;
    if (p == nullptr) return 1;
    @autoreleasepool {
        const int n = 4096;
        p->sync();
        int64_t *u = (int64_t *) [p->bProbe contents];
        std::vector<int64_t> u0(n);
        // Probes spanning several buckets either side of zero, dense near the
        // wrap boundaries where an off-by-one would live.
        for (int i = 0; i < n; i++) {
            const int64_t k = (int64_t) (i % 7) - 3;
            u0[i] = k * bucket_ticks + ((int64_t) i * 2097169) - (int64_t) n;
            u[i] = u0[i];
        }

        const MTLSize grid = MTLSizeMake(n, 1, 1);
        const MTLSize tgp = MTLSizeMake(256, 1, 1);
        auto shift = [&](long long s) {
            id<MTLComputeCommandEncoder> e = p->encoder();
            [e setComputePipelineState:p->pWShift];
            [e setBuffer:p->bProbe offset:0 atIndex:0];
            [e setBytes:&s length:sizeof(s) atIndex:1];
            [e dispatchThreads:grid threadsPerThreadgroup:tgp];
        };
        auto phases = [&](std::vector<float> &out) {
            id<MTLComputeCommandEncoder> e = p->encoder();
            [e setComputePipelineState:p->pWPhase];
            [e setBuffer:p->bProbe offset:0 atIndex:0];
            [e setBuffer:p->bPh offset:0 atIndex:1];
            [e dispatchThreads:grid threadsPerThreadgroup:tgp];
            p->sync();
            out.assign((const float *) [p->bPh contents],
                       (const float *) [p->bPh contents] + n);
        };

        std::vector<float> ph0, ph1;
        phases(ph0);
        shift((long long) bucket_ticks);
        phases(ph1);
        shift(-(long long) bucket_ticks);
        p->sync();

        // Exact assertions, not tolerances: the bucket round trip returns the
        // accumulators bit for bit, and the extracted phase never saw the shift.
        for (int i = 0; i < n; i++) {
            if (u[i] != u0[i]) return 1;
            if (memcmp(&ph0[i], &ph1[i], 4) != 0) return 1;
        }
        return 0;
    }
}

double luc_dev_seconds (void)
{
    return (gImpl != nullptr) ? gImpl->busy : 0;
}

int luc_dev_timing (int on, char *reason, int reason_len)
{
    Impl *p = gImpl;
    if (p == nullptr) {
        put_str(reason, reason_len, "device not initialized");
        return 1;
    }
    if (on == 0) {
        p->timing = false;
        return 0;
    }

    // Stage-boundary sampling times an encoder, which is why a pass gets one of its
    // own. Dispatch-boundary sampling would time each dispatch inside a single encoder
    // and is absent here, so the refusal names what the device does carry.

    if (![p->dev supportsCounterSampling:MTLCounterSamplingPointAtStageBoundary]) {
        put_str(reason, reason_len, "this device samples no counters at an encoder boundary");
        return 1;
    }
    if (p->bSample == nil) {
        id<MTLCounterSet> ts = nil;
        for (id<MTLCounterSet> cs in [p->dev counterSets])
            if ([[cs name] isEqualToString:MTLCommonCounterSetTimestamp]) ts = cs;
        if (ts == nil) {
            put_str(reason, reason_len, "this device carries no timestamp counter set");
            return 1;
        }
        MTLCounterSampleBufferDescriptor *sd = [MTLCounterSampleBufferDescriptor new];
        sd.counterSet = ts;
        sd.sampleCount = (NSUInteger) (2 * Impl::kMaxPass);
        sd.storageMode = MTLStorageModeShared;
        NSError *err = nil;
        p->bSample = [p->dev newCounterSampleBufferWithDescriptor:sd error:&err];
        if (p->bSample == nil) {
            put_str(reason, reason_len, err != nil
                    ? std::string([[err localizedDescription] UTF8String])
                    : std::string("the counter sample buffer would not allocate"));
            return 1;
        }
    }
    p->timing = true;
    return 0;
}

int luc_dev_pass_seconds (double *sec, int64_t *count, int n)
{
    Impl *p = gImpl;
    if (p == nullptr) return 0;
    if (n > LUC_DEV_PASS_N) n = LUC_DEV_PASS_N;
    for (int i = 0; i < n; i++) {
        if (sec != nullptr) sec[i] = p->passSec[i];
        if (count != nullptr) count[i] = p->passCount[i];
    }
    return p->earlyCommits;
}

int64_t luc_dev_bytes (void)
{
    return (gImpl != nullptr) ? gImpl->bytes : 0;
}

}   // extern "C"
