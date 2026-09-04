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

// Ticks of phase: 2^32 per radiation period. The extraction takes the low 32
// bits as a SIGNED fraction, so the working angle sits in [-pi, pi) where FP32
// spacing is finest, and is exactly invariant under whole-period shifts.
constant float TICKS_PER_RAD = 683565275.576f;   // 2^32 / 2 pi
constant float RAD_PER_TICK  = 1.46291808e-9f;   // 2 pi / 2^32

inline float phase_of (long u){
    return float(as_type<int>(uint(u & 0xffffffffL))) * RAD_PER_TICK;
}

inline float2 cmul (float2 a, float2 b){ return float2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x); }

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
// per-slice base rotator, bilinear scatter. The accumulation is a device
// atomic add per corner; the order threads reach a cell is not fixed, so two
// runs of the same step differ in the last bit or two of the source, exactly
// as both reference backends do (manual/GPU.md records it). Every member of the
// set is written from the one particle at its own phase h*theta and scale,
// fel_field_step's harm*theta, into its own source plane.

struct DepPar {
    float gridmax, dgrid, kx, ky, ax, ay, cos_t, sin_t, gam0;
    uint  ngrid, npart, nslice, first, mutate, nf;
};

kernel void zero_src (device float* s [[buffer(0)]], uint gid [[thread_position_in_grid]]){
    s[gid] = 0.0f;
}

kernel void deposit (device atomic_float* S [[buffer(0)]],
                     const device float* X [[buffer(1)]], const device float* Y [[buffer(2)]],
                     const device float* G [[buffer(3)]], const device long* U [[buffer(4)]],
                     const device float* W [[buffer(5)]],
                     const device float2* BASE [[buffer(6)]],
                     constant DepPar& P [[buffer(7)]],
                     constant SetPar& SP [[buffer(8)]],
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
    float sq = sqrt(1.0f + P.kx * ddx * ddx + P.ky * ddy * ddy);

    float d = phase_of(U[gid]);
    // The check's mutation reaches the source row here: coarsen the residual
    // angle by eight mantissa bits, fel_fp32_mod's own hook transcribed.
    if (P.mutate != 0u) {
        float sp = fabs(d) * 3.05175781e-5f;    // 256 ulps of a [-pi,pi) angle
        if (sp > 0.0f) d = rint(d / sp) * sp;
    }

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

        uint idx = (m * P.nslice + fs) * nn + cell;
        float w;
        uint dcell;
        w = wx * wy;                     dcell = 2u * idx;
        atomic_fetch_add_explicit(&S[dcell],      w * cp.x, memory_order_relaxed);
        atomic_fetch_add_explicit(&S[dcell + 1u], w * cp.y, memory_order_relaxed);
        w = (1.0f - wx) * wy;            dcell = 2u * (idx + 1u);
        atomic_fetch_add_explicit(&S[dcell],      w * cp.x, memory_order_relaxed);
        atomic_fetch_add_explicit(&S[dcell + 1u], w * cp.y, memory_order_relaxed);
        w = wx * (1.0f - wy);            dcell = 2u * (idx + P.ngrid);
        atomic_fetch_add_explicit(&S[dcell],      w * cp.x, memory_order_relaxed);
        atomic_fetch_add_explicit(&S[dcell + 1u], w * cp.y, memory_order_relaxed);
        w = (1.0f - wx) * (1.0f - wy);   dcell = 2u * (idx + P.ngrid + 1u);
        atomic_fetch_add_explicit(&S[dcell],      w * cp.x, memory_order_relaxed);
        atomic_fetch_add_explicit(&S[dcell + 1u], w * cp.y, memory_order_relaxed);
    }
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
    float gridmax, dgrid, kx, ky, ax, ay, cos_t, sin_t, gam0;
    uint32_t ngrid, npart, nslice, first, mutate, nf;
};

struct AddPar {
    uint32_t ppm, nslice, npol, pad;
    float pol[4];     // float2 pol[2]
};

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
    id<MTLBuffer> bField, bSrc, bExpK, bTw, bBase, bBaseDep;
    id<MTLBuffer> bProbe, bPh;

    id<MTLComputePipelineState> pTrk, pPush, pZero, pDep, pRow, pRowM, pCol, pColA;
    id<MTLComputePipelineState> pWShift, pWPhase;

    int64_t bytes {0};
    double busy {0};

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
    // Not thread-safe in general: encoder() and a draining sync() both touch enc.
    // What the header's concurrency contract rests on is the narrow case: after one
    // serial drain with nothing encoded since, enc is nil, concurrent transfer calls
    // all take the early return below, and their memcpys of disjoint regions of the
    // shared-storage buffers race nothing.
    void sync (){
        if (enc == nil) return;
        [enc endEncoding];
        [cb commit];
        [cb waitUntilCompleted];
        busy += [cb GPUEndTime] - [cb GPUStartTime];
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
        p->bExpK = alloc((size_t) nfield * nn * 8);
        p->bTw = alloc((size_t) ngrid * 8);
        p->bBase = alloc((size_t) nfield * nslice * 8);
        p->bBaseDep = alloc((size_t) nfield * nslice * 8);
        p->bProbe = alloc(4096 * 8);
        p->bPh = alloc(4096 * 4);
        if (p->bX == nil || p->bU == nil || p->bField == nil || p->bSrc == nil) {
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
        p->pWShift = pso(@"wrap_shift");
        p->pWPhase = pso(@"wrap_phase");
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
    memcpy(s, (const float *) [p->bSrc contents] + p->srcPlane(im, is) * nn * 2, nn * 8);
}

void luc_dev_set_kernel (int im, const float *expk)
{
    Impl *p = gImpl;
    p->sync();
    const size_t nn = (size_t) p->ngrid * p->ngrid;
    memcpy((float *) [p->bExpK contents] + (size_t) im * nn * 2, expk, nn * 8);
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

        DepPar D;
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

        id<MTLComputeCommandEncoder> e = p->encoder();

        auto encTrk = [&](uint32_t leading) {
            T.leading = leading;
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

        [e setComputePipelineState:p->pZero];
        [e setBuffer:p->bSrc offset:0 atIndex:0];
        [e dispatchThreads:MTLSizeMake((size_t) p->nfield * p->nslice * nn * 2, 1, 1)
             threadsPerThreadgroup:tgp];

        [e setComputePipelineState:p->pDep];
        [e setBuffer:p->bSrc offset:0 atIndex:0];
        [e setBuffer:p->bX offset:0 atIndex:1];
        [e setBuffer:p->bY offset:0 atIndex:2];
        [e setBuffer:p->bG offset:0 atIndex:3];
        [e setBuffer:p->bU offset:0 atIndex:4];
        [e setBuffer:p->bW offset:0 atIndex:5];
        [e setBuffer:p->bBaseDep offset:0 atIndex:6];
        [e setBytes:&D length:sizeof(D) atIndex:7];
        [e setBytes:&S length:sizeof(S) atIndex:8];
        [e dispatchThreads:grid threadsPerThreadgroup:tgp];

        // The transform over every plane of the set: tg.y indexes the plane.
        const MTLSize rowTG = MTLSizeMake((size_t) (p->ngrid / p->rowsPerTG), nplane, 1);
        const MTLSize rowT = MTLSizeMake((size_t) (p->rowsPerTG * p->lanes), 1, 1);
        const MTLSize colTG = MTLSizeMake((size_t) (p->ngrid / p->colsPerTG), nplane, 1);
        const MTLSize colT = MTLSizeMake((size_t) (p->colsPerTG * p->lanes), 1, 1);

        [e setComputePipelineState:p->pRow];
        [e setBuffer:p->bField offset:0 atIndex:0];
        [e setBuffer:p->bTw offset:0 atIndex:1];
        [e setBytes:&fwd length:4 atIndex:2];
        [e dispatchThreadgroups:rowTG threadsPerThreadgroup:rowT];

        [e setComputePipelineState:p->pCol];
        [e setBuffer:p->bField offset:0 atIndex:0];
        [e setBuffer:p->bTw offset:0 atIndex:1];
        [e setBytes:&fwd length:4 atIndex:2];
        [e setBytes:&one length:4 atIndex:3];
        [e dispatchThreadgroups:colTG threadsPerThreadgroup:colT];

        [e setComputePipelineState:p->pRowM];
        [e setBuffer:p->bField offset:0 atIndex:0];
        [e setBuffer:p->bTw offset:0 atIndex:1];
        [e setBytes:&inv length:4 atIndex:2];
        [e setBuffer:p->bExpK offset:0 atIndex:3];
        [e setBytes:&ppm length:4 atIndex:4];
        [e dispatchThreadgroups:rowTG threadsPerThreadgroup:rowT];

        [e setComputePipelineState:p->pColA];
        [e setBuffer:p->bField offset:0 atIndex:0];
        [e setBuffer:p->bTw offset:0 atIndex:1];
        [e setBytes:&inv length:4 atIndex:2];
        [e setBytes:&nrm length:4 atIndex:3];
        [e setBuffer:p->bSrc offset:0 atIndex:4];
        [e setBytes:&A length:sizeof(A) atIndex:5];
        [e dispatchThreadgroups:colTG threadsPerThreadgroup:colT];
    }
    return 0;
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

int64_t luc_dev_bytes (void)
{
    return (gImpl != nullptr) ? gImpl->bytes : 0;
}

}   // extern "C"
