#!/usr/bin/env node
/*
 * midiprec.js — "raw precision (24h)" of a MIDI file played on a 128-key piano.
 *
 * NaNDL precision model (w_i = N_i/f, s_i = ½·w_i·L·Πλ, p_i = P(|X| <= s_i), X ~ N(0,1),
 * E[T_C] = E[T_A]/P(C), precision L* such that E[T_C] = 24 h), translated to MIDI:
 *   inputs  = hand gestures: every tick with note-ons is split between the two hands;
 *             overlapped copies (same key, same tick) are extra "swift" re-strikes
 *   windows = time to the neighbouring notes (ticks -> seconds via PPQ + tempo map),
 *             at most ±50 ms; divided by (velocity/127) so soft notes matter less
 *   λ       = piano physics: stretch / clusters / arm length, arm spread, Fitts jumps,
 *             hand speed (glissando aware), same-key repetition (capped at 1e4)
 * Reads .mid and .zip/.7z/.xz/.gz/.bz2/.zst/.lzma/.br/.tar(.*), .rar via 7-Zip.
 * Single file, no dependencies, Node.js >= 18 (22+ for .zst). See README.md.
 */
'use strict';
const fs = require('fs');
const path = require('path');
const os = require('os');

const VERSION = '1.0.0';

// =====================================================================================
// Math: erf / Phi(v) = -ln(erf(x/sqrt2)), v = ln x, and fast tables
// =====================================================================================
const SQRT2 = Math.SQRT2;
const SQRT_2_PI = Math.sqrt(2 / Math.PI);
const INV_SQRT_PI = 1 / Math.sqrt(Math.PI);
const LN2 = Math.LN2;
const HALF_LN_2_PI = 0.5 * Math.log(2 / Math.PI);

function erfSeries(z) { // erf(z), non-alternating series, accurate for 0 <= z < ~3
  const z2 = z * z;
  let term = z, sum = z;
  for (let n = 1; n < 1000; n++) {
    term *= 2 * z2 / (2 * n + 1);
    sum += term;
    if (term < sum * 1e-17) break;
  }
  return 2 * INV_SQRT_PI * Math.exp(-z2) * sum;
}
function erfcCF(z) { // erfc(z) for z >= 2 (continued fraction, modified Lentz)
  const tiny = 1e-300;
  let f = z, C = z, D = 0;
  for (let j = 1; j < 10000; j++) {
    const a = j * 0.5;
    D = z + a * D; if (D === 0) D = tiny; D = 1 / D;
    C = z + a / C; if (C === 0) C = tiny;
    const d = C * D;
    f *= d;
    if (Math.abs(d - 1) < 2e-16) break;
  }
  return Math.exp(-z * z) * INV_SQRT_PI / f;
}
function erfcHP(z) { return z < 2 ? 1 - erfSeries(z) : erfcCF(z); }

// Exact Phi and derivatives w.r.t. v = ln x.   Phi(v) = F(x) = -ln P(|X| <= x), X~N(0,1)
function phiExact(x) {
  if (!(x > 0)) return Infinity;
  if (x > 40) return 0;
  const z = x / SQRT2;
  if (z < 1) return -Math.log(erfSeries(z));
  return -Math.log1p(-erfcHP(z));
}
function dphiExact(x) { // dPhi/dv = -x E'(x)/E(x)
  if (x > 40) return 0;
  const z = x / SQRT2;
  const E = z < 1 ? erfSeries(z) : 1 - erfcHP(z);
  return -x * SQRT_2_PI * Math.exp(-0.5 * x * x) / E;
}

// Table of psi(v) = ln Phi(v) (+ derivatives) for cubic Hermite interpolation.
const TV0 = -12, TV1 = 2.5, TSTEP = 1 / 128;
const TN = Math.round((TV1 - TV0) / TSTEP) + 1;
const T_PSI = new Float64Array(TN), T_D1 = new Float64Array(TN), T_D2 = new Float64Array(TN);
(function buildTable() {
  for (let i = 0; i < TN; i++) {
    const v = TV0 + i * TSTEP, x = Math.exp(v);
    const p = phiExact(x), d1 = dphiExact(x);
    const d2 = d1 * (1 - x * x) + d1 * d1; // Phi''(v)
    T_PSI[i] = Math.log(p);
    T_D1[i] = d1 / p;
    T_D2[i] = d2 / p - (d1 / p) * (d1 / p);
  }
})();
const V_CUT = TV1; // Φ(v) < 1e-32 beyond x≈12: treated as 0
// Fast Phi(v); also writes derivatives into PH[1], PH[2]
const PH = new Float64Array(3);
function phiV(v) {
  if (v >= TV1) { PH[0] = 0; PH[1] = 0; PH[2] = 0; return 0; } // Φ < 1e-32 beyond x≈12
  if (v < TV0) {
    const x2 = Math.exp(2 * v);
    const p = -v - HALF_LN_2_PI + x2 / 6;
    PH[0] = p; PH[1] = -1 + x2 / 3; PH[2] = 2 * x2 / 3; return p;
  }
  const f = (v - TV0) / TSTEP;
  let i = f | 0; if (i >= TN - 1) i = TN - 2;
  const t = f - i, t2 = t * t, t3 = t2 * t;
  const p0 = T_PSI[i], p1 = T_PSI[i + 1], m0 = T_D1[i] * TSTEP, m1 = T_D1[i + 1] * TSTEP;
  const psi = (2 * t3 - 3 * t2 + 1) * p0 + (t3 - 2 * t2 + t) * m0 + (-2 * t3 + 3 * t2) * p1 + (t3 - t2) * m1;
  const dpsi = ((6 * t2 - 6 * t) * p0 + (3 * t2 - 4 * t + 1) * m0 + (-6 * t2 + 6 * t) * p1 + (3 * t2 - 2 * t) * m1) / TSTEP;
  const d2psi = T_D2[i] + (T_D2[i + 1] - T_D2[i]) * t;
  const p = Math.exp(psi);
  PH[0] = p; PH[1] = p * dpsi; PH[2] = p * (d2psi + dpsi * dpsi);
  return p;
}
// Exact version with the same output convention (used for the final answer)
function phiVexact(v) {
  const x = Math.exp(v);
  if (x > 40) { PH[0] = 0; PH[1] = 0; PH[2] = 0; return 0; }
  const p = phiExact(x), d1 = dphiExact(x);
  PH[0] = p; PH[1] = d1; PH[2] = d1 * (1 - x * x) + d1 * d1;
  return p;
}

// float64 -> bits (bin indices without log())
const F64 = new Float64Array(1), U32 = new Uint32Array(F64.buffer);
F64[0] = 1; if (U32[1] !== 0x3FF00000) throw new Error('big-endian platforms are not supported');

// =====================================================================================
// Model parameters ("average trained pianist")
// =====================================================================================
const DEFAULTS = {
  hours: 24,          // completion target (NaNDL uses 24h)
  restart: 1.0,       // s per attempt overhead (reset/respawn); attempt time counts from the first note
  window: 0.050,      // s, max half-window (isolated notes may be off by ±50 ms)
  fingers: 5,         // keys (finger units) one hand can press at once
  reach: 8,           // white keys, comfortable hand stretch (a 9th)
  betaFingers: 1.5, betaReach: 2,
  armReach: 30,       // white keys (~70 cm): longest cluster one arm can lie on
  spread: 100,        // cm, comfortable distance between both hands' outer keys
  betaSpread: 2,
  armSpan: 170,       // cm, fingertip-to-fingertip limit for both hands together
  wallPower: 8,       // steepness of the arm-length / arm-span walls
  engage: 0.25,       // s, the other hand counts as "in use" for this long
  keyWidth: 2.35,     // cm, white key width
  fittsA: 0.040,      // s, Fitts intercept   MT = a + b·log2(1 + D/W)
  fittsB: 0.050,      // s/bit, Fitts slope
  handRate: 12,       // gestures/s one hand can sustain comfortably
  glissRate: 40,      // gestures/s while gliding over adjacent same-colour keys
  rateTau: 0.5,       // s, averaging time for hand speed
  betaRate: 1,
  repeatTime: 0.080,  // s, fastest comfortable re-strike of the same key
  betaRepeat: 1,
  cross: Math.LN2,    // hand-crossing penalty (λ = 0.5)
  twoHands: 1.0,      // preference for one hand when both work (selection only)
  maxPenalty: 1e4,    // cap of the total physical multiplier (1/λ)
  velocityPower: 1,   // window ÷ (velocity/127)^velocityPower  (soft notes matter less)
};

// Key geometry (white-key units)
const KX = new Float64Array(128), KBLACK = new Uint8Array(128), KIDX = new Int32Array(128);
(function () {
  const PCX = [0, 0.5, 1, 1.5, 2, 3, 3.5, 4, 4.5, 5, 5.5, 6];
  const PCB = [0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0];
  let w = 0, b = 0;
  for (let k = 0; k < 128; k++) {
    KX[k] = 7 * Math.floor(k / 12) + PCX[k % 12];
    KBLACK[k] = PCB[k % 12];
    KIDX[k] = PCB[k % 12] ? b++ : w++;
  }
})();

// log tables indexed by half-white-key steps / finger units
const LOG_HALF = new Float64Array(512), LOG_HALF1 = new Float64Array(512), LOG_INT = new Float64Array(257);
for (let i = 0; i < 512; i++) { LOG_HALF[i] = Math.log(i / 2); LOG_HALF1[i] = Math.log(i / 2 + 1); }
for (let i = 0; i < 257; i++) LOG_INT[i] = Math.log(i);
function logHalf(v) { const q = v * 2; return (q < 512 && q === Math.floor(q)) ? LOG_HALF[q] : Math.log(v); }

// =====================================================================================
// Hand model: splits an onset between two hands, returns physical multipliers μ
// =====================================================================================
class HandModel {
  constructor(P) {
    this.P = P;
    this.anchor = new Float64Array([KX[52], KX[76]]); // E3 / E5
    this.lastT = new Float64Array([-Infinity, -Infinity]);
    this.load = new Float64Array(2);
    this.lastLo = new Float64Array([KX[52], KX[76]]);
    this.lastHi = new Float64Array([KX[52], KX[76]]);
    this.gKey = new Int32Array([-1, -1]);
    this.gDir = new Int32Array(2);
    this.lastPress = new Float64Array(128).fill(-Infinity);
    this.K = new Int32Array(128);
    this.X = new Float64Array(128);
    this.fuPre = new Int32Array(129); this.fuSuf = new Int32Array(129);
    this.dPre = new Float64Array(129); this.dSuf = new Float64Array(129); this.DL = new Float64Array(128);
    this.RS = new Int32Array(128); this.RL = new Int32Array(128);
    this.rpP = new Float64Array(129); this.rpS = new Float64Array(129); this.rk = new Float64Array(128);
    this.cFuS = new Float64Array(16); this.cRL = new Float64Array(16); this.cRR = new Float64Array(16);
    this.nc = 0; this.r1 = 0; this.pm = 0; this.pw0 = 0; this.pw1 = 0; this.pw2 = 0; this.pw3 = 0;
    this.rmD = NaN; this.rmC = 0;
    this.cand = new Int32Array(160); this.g1 = new Float64Array(4); this.gj = new Int32Array(4);
    this.mu = new Float64Array(2); this.nG = 0;
    this.splitJ = 0; this.muH = new Float64Array(2); // keys [0,splitJ) -> left hand; μ per hand
    this.dt = new Float64Array(2); this.ld = new Float64Array(2);
    this.rcN = new Float64Array(2); this.rcG = new Float64Array(2);
    this.adj = new Uint8Array(2); this.adjDir = new Int32Array(2);
    this.F = P.fingers; this.S = P.reach; this.half = P.reach / 2;
    this.lnF = Math.log(P.fingers); this.lnS = Math.log(P.reach); this.lnSp = Math.log(P.spread); this.lnKW = Math.log(P.keyWidth);
    this.invArm = 1 / P.armReach; this.invSpan = 1 / P.armSpan; this.wp8 = P.wallPower === 8;
    this.invFB = LN2 / P.fittsB; this.invTau = 1 / P.rateTau;
    this.capLn = Math.log(P.maxPenalty);
    this.fittsTab = new Float64Array(1024);
    for (let i = 1; i < 1024; i++) this.fittsTab[i] = P.fittsA + P.fittsB * Math.log2(1 + i / 4);
    this._na = 0;
    // statistics: committed raw cost per factor [reach, move, speed, repeat, spread, cross]
    this.cs = new Float64Array(6);
    this.stats = { gestures: 0, twoHand: 0, capped: 0, sumCost: 0 };
  }
  wall(x) { if (this.wp8) { const x2 = x * x, x4 = x2 * x2; return x4 * x4; } return Math.pow(x, this.P.wallPower); }
  // stretch / too many keys for one hand (cluster) / longer than an arm
  reachCost(u, span) {
    const P = this.P;
    let c = 0;
    if (u > this.F) c = P.betaFingers * ((u < 257 ? LOG_INT[u] : Math.log(u)) - this.lnF);
    if (span > this.S) c += P.betaReach * (logHalf(span) - this.lnS) + this.wall(span * this.invArm);
    return c;
  }
  // both arms in use, outer keys L..H (white-key units)
  spreadCost(L, H) {
    const P = this.P, w = H - L, dcm = (w + 1) * P.keyWidth;
    if (dcm <= P.spread) return 0;
    const q = w * 2;
    return P.betaSpread * ((q < 512 && q === Math.floor(q) ? LOG_HALF1[q] : Math.log(w + 1)) + this.lnKW - this.lnSp) + this.wall(dcm * this.invSpan);
  }
  // Fitts: time needed to move hand h so that [lo,hi] is under it vs. time available; sets this._na
  moveCost(h, lo, hi) {
    const anc = this.anchor[h];
    let na = anc;
    if (hi - lo <= this.S) { if (lo < anc - this.half) na = lo + this.half; else if (hi > anc + this.half) na = hi - this.half; }
    else na = 0.5 * (lo + hi);
    this._na = na;
    const d = na > anc ? na - anc : anc - na;
    if (d === 0) return 0;
    const q = d * 4;
    const P = this.P;
    const MT = (q < 1024 && q === Math.floor(q)) ? this.fittsTab[q] : P.fittsA + P.fittsB * Math.log2(1 + d);
    const dt = this.dt[h];
    return dt < MT ? (MT - dt) * this.invFB : 0;
  }
  repCost(delta) {
    const P = this.P;
    if (delta >= P.repeatTime) return 0;
    return Math.log(1 + P.betaRepeat * Math.log(P.repeatTime / (delta > 1e-300 ? delta : 1e-300)));
  }
  crossCost(h) { return (h === 0 ? this._na > this.anchor[1] : this._na < this.anchor[0]) ? this.P.cross : 0; }
  commit(h, t, lo, hi, single, key, cr, cm, cq, cp, cs, cx) {
    this.anchor[h] = this._na; this.lastT[h] = t; this.load[h] = this.ld[h];
    this.lastLo[h] = lo; this.lastHi[h] = hi;
    if (single) { this.gDir[h] = this.adj[h] ? this.adjDir[h] : 0; this.gKey[h] = key; }
    else { this.gKey[h] = -1; this.gDir[h] = 0; }
    const C = this.cs;
    C[0] += cr; C[1] += cm; C[2] += cq; C[3] += cp; C[4] += cs; C[5] += cx;
    let c = cr + cm + cq + cp + cs + cx;
    if (c > this.capLn) { c = this.capLn; this.stats.capped++; }   // physical penalty cap (maxPenalty)
    this.stats.sumCost += c;
    return Math.exp(-c);
  }
  // rep cost with a one-entry memo (many keys share the same re-strike interval)
  repM(delta) {
    if (delta === this.rmD) return this.rmC;
    this.rmD = delta; this.rmC = this.repCost(delta);
    return this.rmC;
  }
  // words: Int32Array(4) key bitmap; t: onset time (s). Fills this.mu / this.nG; returns #keys
  onset(words, t) {
    const P = this.P, K = this.K, X = this.X, lp = this.lastPress, DL = this.DL;
    const fuPre = this.fuPre, dPre = this.dPre, dSuf = this.dSuf, RS = this.RS, RL = this.RL;
    const w0 = words[0], w1 = words[1], w2 = words[2], w3 = words[3];
    let m;
    // ---- chord geometry (cached when the same key set repeats)
    if (w0 === this.pw0 && w1 === this.pw1 && w2 === this.pw2 && w3 === this.pw3) {
      m = this.pm;
      let dmin = Infinity;
      for (let j = 0; j < m; j++) { const k = K[j], dl = t - lp[k]; DL[j] = dl; if (dl < dmin) dmin = dl; dPre[j + 1] = dmin; lp[k] = t; }
    } else {
      m = 0;
      for (let w = 0; w < 4; w++) {
        let bits = words[w];
        while (bits !== 0) {
          const low = bits & -bits;
          K[m++] = (w << 5) | (31 - Math.clz32(low));
          bits ^= low;
        }
      }
      // forward pass: positions, runs of adjacent keys (finger units: 2 adjacent keys per finger)
      let done = 0, run = 0, rs = 0, dmin = Infinity;
      fuPre[0] = 0;
      for (let j = 0; j < m; j++) {
        const k = K[j], x = KX[k];
        X[j] = x;
        if (j > 0 && x - X[j - 1] <= 1) run++;
        else { if (run > 0) RL[rs] = run; done += (run + 1) >> 1; run = 1; rs = j; }
        RS[j] = rs;
        fuPre[j + 1] = done + ((run + 1) >> 1);
        const dl = t - lp[k]; DL[j] = dl; if (dl < dmin) dmin = dl; dPre[j + 1] = dmin; lp[k] = t;
      }
      RL[rs] = run;
      this.pm = m; this.pw0 = w0; this.pw1 = w1; this.pw2 = w2; this.pw3 = w3;
      // split candidates
      const cand = this.cand;
      let nc = 0;
      if (m >= 2) {
        if (m <= 8) { for (let j = 1; j < m; j++) cand[nc++] = j; }
        else {
          const lo = X[0], hi = X[m - 1], mid = 0.5 * (lo + hi);
          let jb = 1; while (jb < m - 1 && X[jb] <= mid) jb++;
          cand[nc++] = jb; if ((m >> 1) !== jb) cand[nc++] = m >> 1;
          const g1 = this.g1, gj = this.gj; g1[0] = g1[1] = g1[2] = g1[3] = 0; gj[0] = gj[1] = gj[2] = gj[3] = 1;
          for (let j = 1; j < m; j++) {
            const g = X[j] - X[j - 1];
            if (g > g1[3]) {
              let q = 3; while (q > 0 && g > g1[q - 1]) { g1[q] = g1[q - 1]; gj[q] = gj[q - 1]; q--; }
              g1[q] = g; gj[q] = j;
            }
          }
          for (let q = 0; q < 4; q++) if (g1[q] > 0 && gj[q] !== jb && gj[q] !== (m >> 1)) cand[nc++] = gj[q];
        }
      }
      this.nc = nc;
      // per-candidate geometry: suffix finger units and reach costs
      const fuTot = fuPre[m];
      for (let ci = 0; ci < nc; ci++) {
        const j = this.cand[ci];
        let fs = fuTot - fuPre[j];
        if (RS[j] === RS[j - 1]) { // the run crossing the split is cut in two
          const L = RL[RS[j]], L1 = j - RS[j];
          fs += ((L1 + 1) >> 1) - ((L + 1) >> 1) + ((L - L1 + 1) >> 1);
        }
        this.cFuS[ci] = fs;
        this.cRL[ci] = this.reachCost(fuPre[j], X[j - 1] - X[0]);
        this.cRR[ci] = this.reachCost(fs, X[m - 1] - X[j]);
      }
      this.r1 = this.reachCost(fuTot, X[m - 1] - X[0]);
    }
    const single = m === 1, key0 = K[0], nc = this.nc;
    // re-strike cost of each part = cost of its smallest interval (monotone) → prefix/suffix maxima
    const rpP = this.rpP, rpS = this.rpS, small = m <= 8;
    if (small) {
      const RT = P.repeatTime;
      rpP[0] = 0;
      for (let j = 0; j < m; j++) { const dl = DL[j]; const c = dl < RT ? this.repM(dl) : 0; this.rk[j] = c; rpP[j + 1] = c > rpP[j] ? c : rpP[j]; }
      if (nc > 0) { rpS[m] = 0; for (let j = m - 1; j >= 0; j--) { const c = this.rk[j]; rpS[j] = c > rpS[j + 1] ? c : rpS[j + 1]; } }
    } else if (nc > 0) { let dmin = Infinity; for (let j = m - 1; j >= 0; j--) { if (DL[j] < dmin) dmin = DL[j]; dSuf[j] = dmin; } }
    // per-hand speed (EWMA gestures/s), glissando check — same for every option
    for (let h = 0; h < 2; h++) {
      const dt = t - this.lastT[h];
      this.dt[h] = dt;
      const ld = this.load[h] * Math.exp(-dt * this.invTau) + 1;
      this.ld[h] = ld;
      const r = (ld - 0.5) * this.invTau;
      const rn = r > P.handRate ? Math.log(1 + P.betaRate * Math.log(r / P.handRate)) : 0;
      this.rcN[h] = rn; this.rcG[h] = rn; this.adj[h] = 0;
      if (single) {
        const g = this.gKey[h];
        if (g >= 0 && dt < 0.25 && KBLACK[g] === KBLACK[key0]) {
          const di = KIDX[key0] - KIDX[g];
          if (di === 1 || di === -1) {
            this.adj[h] = 1; this.adjDir[h] = di;
            if (this.gDir[h] === 0 || this.gDir[h] === di) this.rcG[h] = r > P.glissRate ? Math.log(1 + P.betaRate * Math.log(r / P.glissRate)) : 0;
          }
        }
      }
    }
    const lo = X[0], hi = X[m - 1];
    const eng0 = this.dt[0] < P.engage, eng1 = this.dt[1] < P.engage;
    // ---- option 0/1: one hand plays everything
    const rpAll = small ? rpP[m] : this.repM(dPre[m]);
    const base1 = this.r1 + rpAll;
    const mvL = this.moveCost(0, lo, hi);
    let spL = 0, cxL = 0;
    if (eng1) { spL = this.spreadCost(lo < this.lastLo[1] ? lo : this.lastLo[1], hi > this.lastHi[1] ? hi : this.lastHi[1]); cxL = this.crossCost(0); }
    const cL = base1 + mvL + this.rcG[0] + spL + cxL;
    const mvR = this.moveCost(1, lo, hi);
    let spR = 0, cxR = 0;
    if (eng0) { spR = this.spreadCost(lo < this.lastLo[0] ? lo : this.lastLo[0], hi > this.lastHi[0] ? hi : this.lastHi[0]); cxR = this.crossCost(1); }
    const cR = base1 + mvR + this.rcG[1] + spR + cxR;
    let best = cL <= cR ? 0 : 1, bestCost = best === 0 ? cL : cR, bestC = -1;
    // ---- option 2: left hand takes keys < j, right hand keys >= j
    let sp2 = 0;
    if (nc > 0) {
      sp2 = this.spreadCost(lo, hi);
      const base2 = P.twoHands + 2 * sp2 + this.rcN[0] + this.rcN[1];
      const cand = this.cand;
      for (let ci = 0; ci < nc; ci++) {
        const j = cand[ci];
        const c0 = base2 + this.cRL[ci] + (small ? rpP[j] : this.repM(dPre[j])) + this.moveCost(0, lo, X[j - 1]);
        if (c0 >= bestCost) continue;
        const c = c0 + this.cRR[ci] + (small ? rpS[j] : this.repM(dSuf[j])) + this.moveCost(1, X[j], hi);
        if (c < bestCost) { bestCost = c; best = 2; bestC = ci; }
      }
    }
    // ---- commit the chosen option
    const st = this.stats;
    if (best < 2) {
      const h = best;
      this.moveCost(h, lo, hi); // sets _na
      this.mu[0] = h === 0
        ? this.commit(0, t, lo, hi, single, key0, this.r1, mvL, this.rcG[0], rpAll, spL, cxL)
        : this.commit(1, t, lo, hi, single, key0, this.r1, mvR, this.rcG[1], rpAll, spR, cxR);
      this.nG = 1;
      this.splitJ = h === 0 ? m : 0; this.muH[h] = this.mu[0]; this.muH[1 - h] = 0;
    } else {
      const j = this.cand[bestC], lHi = X[j - 1], rLo = X[j];
      let cm = this.moveCost(0, lo, lHi);
      this.mu[0] = this.commit(0, t, lo, lHi, false, 0, this.cRL[bestC], cm, this.rcN[0], small ? rpP[j] : this.repM(dPre[j]), sp2, 0);
      cm = this.moveCost(1, rLo, hi);
      this.mu[1] = this.commit(1, t, rLo, hi, false, 0, this.cRR[bestC], cm, this.rcN[1], small ? rpS[j] : this.repM(dSuf[j]), sp2, 0);
      this.nG = 2; st.twoHand++;
      this.splitJ = j; this.muH[0] = this.mu[0]; this.muH[1] = this.mu[1];
    }
    st.gestures += this.nG;
    return m;
  }
}

// =====================================================================================
// Accumulator: histograms of a_i = ½·w_i·μ_i, time segments, online CSV, solver
// =====================================================================================
const FINE_SHIFT = 14, FINE_N = 1 << 17;   // 64 bins per octave
const COARSE_SHIFT = 18, COARSE_N = 1 << 13; // 4 bins per octave
const GRID_STEP = LN2 / 32;

class Accumulator {
  constructor(P, opts) {
    this.P = P;
    this.lnTarget = Math.log(P.hours * 3600);
    this.restart = P.restart;
    // fine histogram (global)
    this.fn = new Float64Array(FINE_N); this.fs = new Float64Array(FINE_N); this.fs2 = new Float64Array(FINE_N);
    this.fDirty = new Uint8Array(FINE_N); this.fu = new Float64Array(FINE_N); this.fv = new Float64Array(FINE_N);
    this.fList = new Int32Array(4096); this.fLen = 0;
    // current batch (one CSV frame)
    this.bn = new Float64Array(COARSE_N); this.bs = new Float64Array(COARSE_N); this.bt = new Float64Array(COARSE_N);
    this.bTouched = new Int32Array(COARSE_N); this.bLen = 0; this.bLastT = 0;
    // segment store
    this.segStart = new Int32Array(1024); this.segLast = new Float64Array(1024); this.nSeg = 0;
    this.eBin = new Uint16Array(1 << 16); this.eN = new Float64Array(1 << 16); this.eU = new Float64Array(1 << 16); this.eT = new Float64Array(1 << 16);
    this.nEnt = 0; this.entCap = opts.entCap || (1 << 21);
    this.segDur = 1e-6; this.segEndT = -Infinity; // fine time buckets, ×4 on each compaction
    // counters
    this.inputs = 0; this.zero = 0; this.tFirst = NaN; this.tLast = NaN;
    // frames
    this.fps = opts.fps || 60;
    this.csv = !!opts.csv;
    this.nextLine = 0; this.nextLineT = 0;
    this.lineVal = this.csv ? new Float64Array(1024) : null;
    this.lineNotes = this.csv ? new Float64Array(1024) : null;
    this.batchNotes = 0; this.totalNotes = 0;
    this.nLines = 0;
    // max NPS (sliding 1 s over frames)
    this.npsRing = new Float64Array(Math.max(1, Math.round(this.fps)));
    this.npsSum = 0; this.npsMax = 0; this.npsMaxT = 0;
    // online solver state
    this.lam = -Infinity; this.curL = 0; this.changed = false;
    this.grid = new Map();
    this.solveCount = 0; this.gridInits = 0;
  }
  // ---- input emission (hot) ----
  addInput(t, a) {
    if (t >= this.nextLineT) this.advance(t);
    if (t >= this.segEndT && this.bLen > 0) this.flushBatch();
    if (this.inputs === 0) this.tFirst = t;
    this.inputs++;
    this.tLast = t;
    if (!(a > 0)) { this.zero++; this.changed = true; return; }
    F64[0] = a;
    const hi = U32[1];
    const fb = hi >>> FINE_SHIFT;
    if (this.fn[fb] === 0) this.fInsert(fb);
    this.fn[fb]++; this.fs[fb] += a; this.fs2[fb] += a * a; this.fDirty[fb] = 1;
    const cb = hi >>> COARSE_SHIFT;
    if (this.bLen === 0) this.segEndT = t + this.segDur;
    if (this.bn[cb] === 0) this.bTouched[this.bLen++] = cb;
    this.bn[cb]++; this.bs[cb] += a; this.bt[cb] += t;
    this.bLastT = t;
    this.changed = true;
  }
  addInputN(t, a, count) { // `count` identical inputs at time t
    if (t >= this.nextLineT) this.advance(t);
    if (t >= this.segEndT && this.bLen > 0) this.flushBatch();
    if (this.inputs === 0) this.tFirst = t;
    this.inputs += count;
    this.tLast = t;
    if (!(a > 0)) { this.zero += count; this.changed = true; return; }
    F64[0] = a;
    const hi = U32[1];
    const fb = hi >>> FINE_SHIFT;
    if (this.fn[fb] === 0) this.fInsert(fb);
    this.fn[fb] += count; this.fs[fb] += a * count; this.fs2[fb] += a * a * count; this.fDirty[fb] = 1;
    const cb = hi >>> COARSE_SHIFT;
    if (this.bLen === 0) this.segEndT = t + this.segDur;
    if (this.bn[cb] === 0) this.bTouched[this.bLen++] = cb;
    this.bn[cb] += count; this.bs[cb] += a * count; this.bt[cb] += t * count;
    this.bLastT = t;
    this.changed = true;
  }
  addNotes(t, n) { // raw note-on count at an onset (for {n}/{nps})
    if (t >= this.nextLineT) this.advance(t);
    this.batchNotes += n;
  }
  fInsert(b) {
    if (this.fLen === this.fList.length) { const n = new Int32Array(this.fLen * 2); n.set(this.fList); this.fList = n; }
    const L = this.fList;
    let lo = 0, hi = this.fLen;
    while (lo < hi) { const mid = (lo + hi) >> 1; if (L[mid] < b) lo = mid + 1; else hi = mid; }
    L.copyWithin(lo + 1, lo, this.fLen);
    L[lo] = b; this.fLen++;
  }
  // close all frames whose time is <= t
  advance(t) {
    this.closeLine();
    this.nextLine++;
    this.nextLineT = this.nextLine / this.fps;
    if (this.nextLineT > t) return;
    // remaining frames up to t are empty: fast-forward
    let kmax = Math.floor(t * this.fps);
    while ((kmax + 1) / this.fps <= t) kmax++;
    while (kmax / this.fps > t) kmax--;
    const k0 = this.nextLine, cnt = kmax - k0 + 1;
    if (cnt > 0) {
      const ring = this.npsRing, rl = ring.length;
      if (cnt >= rl) { ring.fill(0); this.npsSum = 0; }
      else for (let k = k0; k <= kmax; k++) { const ri = k % rl; this.npsSum -= ring[ri]; ring[ri] = 0; }
      if (this.csv) {
        this.ensureLines(kmax + 1);
        this.lineVal.fill(this.curL, k0, kmax + 1);
        this.lineNotes.fill(this.totalNotes, k0, kmax + 1);
        this.nLines = kmax + 1;
      }
      this.nextLine = kmax + 1;
      this.nextLineT = this.nextLine / this.fps;
    }
  }
  ensureLines(n) {
    if (n <= this.lineVal.length) return;
    let c = this.lineVal.length; while (c < n) c *= 2;
    const a = new Float64Array(c); a.set(this.lineVal); this.lineVal = a;
    const b = new Float64Array(c); b.set(this.lineNotes); this.lineNotes = b;
  }
  closeLine() {
    // max NPS bookkeeping (notes in the last 1 s worth of frames)
    const ring = this.npsRing, ri = this.nextLine % ring.length;
    this.npsSum += this.batchNotes - ring[ri]; ring[ri] = this.batchNotes;
    if (this.npsSum > this.npsMax) { this.npsMax = this.npsSum; this.npsMaxT = this.nextLine / this.fps; }
    this.totalNotes += this.batchNotes; this.batchNotes = 0;
    if (this.bLen > 0) this.flushBatch();
    if (this.csv) {
      const k = this.nextLine;
      this.ensureLines(k + 1);
      if (this.changed) { this.solveOnline(); this.changed = false; }
      this.lineVal[k] = this.curL;
      this.lineNotes[k] = this.totalNotes;
      this.nLines = k + 1;
    }
  }
  flushBatch() {
    const n = this.bLen;
    if (this.nEnt + n > this.eBin.length) this.growEnt(this.nEnt + n);
    if (this.nSeg + 1 >= this.segStart.length) this.growSeg();
    const s = this.nSeg;
    this.segStart[s] = this.nEnt;
    let e = this.nEnt;
    for (let i = 0; i < n; i++) {
      const b = this.bTouched[i];
      const cnt = this.bn[b];
      this.eBin[e] = b; this.eN[e] = cnt; this.eU[e] = Math.log(this.bs[b] / cnt); this.eT[e] = this.bt[b];
      this.bn[b] = 0; this.bs[b] = 0; this.bt[b] = 0;
      e++;
    }
    this.nEnt = e; this.segLast[s] = this.bLastT; this.nSeg = s + 1;
    this.segStart[this.nSeg] = e;
    // keep online grid points current
    if (this.grid.size) for (const gp of this.grid.values()) this.segIntegrate(gp, s);
    this.bLen = 0;
    if (this.nEnt > this.entCap) {
      do { this.compact(); } while (this.nEnt > this.entCap * 0.5 && this.nSeg > 1024);
      if (this.nEnt > this.entCap * 0.5) this.entCap = this.nEnt * 2;
    }
  }
  growEnt(min) {
    let c = this.eBin.length; while (c < min) c *= 2;
    const nb = new Uint16Array(c); nb.set(this.eBin); this.eBin = nb;
    let a = new Float64Array(c); a.set(this.eN); this.eN = a;
    a = new Float64Array(c); a.set(this.eU); this.eU = a;
    a = new Float64Array(c); a.set(this.eT); this.eT = a;
  }
  growSeg() {
    const c = this.segStart.length * 2;
    const a = new Int32Array(c); a.set(this.segStart); this.segStart = a;
    const b = new Float64Array(c); b.set(this.segLast); this.segLast = b;
  }
  // merge groups of 4 neighbouring segments (bounded memory for very long / dense MIDIs)
  compact() {
    const sn = new Float64Array(COARSE_N), ss = new Float64Array(COARSE_N), st = new Float64Array(COARSE_N);
    const touched = new Int32Array(COARSE_N);
    let w = 0, ws = 0;
    for (let s = 0; s < this.nSeg; s += 4) {
      const s2 = Math.min(s + 4, this.nSeg);
      let nt = 0;
      for (let e = this.segStart[s]; e < this.segStart[s2]; e++) {
        const b = this.eBin[e], n = this.eN[e];
        if (sn[b] === 0) touched[nt++] = b;
        sn[b] += n; ss[b] += n * Math.exp(this.eU[e]); st[b] += this.eT[e];
      }
      const start = w;
      for (let i = 0; i < nt; i++) {
        const b = touched[i];
        this.eBin[w] = b; this.eN[w] = sn[b]; this.eU[w] = Math.log(ss[b] / sn[b]); this.eT[w] = st[b];
        sn[b] = 0; ss[b] = 0; st[b] = 0; w++;
      }
      this.segStart[ws] = start; this.segLast[ws] = this.segLast[s2 - 1]; ws++;
    }
    this.nSeg = ws; this.nEnt = w; this.segStart[ws] = w;
    this.segDur *= 4; this.compactions = (this.compactions || 0) + 1;
  }
  // integrate survival over one stored segment for grid point gp (table Phi), tracking d/dλ too
  segIntegrate(gp, s) {
    const lam = gp.lam;
    let h = 0, hd = 0, ht = 0, htd = 0;
    const e1 = this.segStart[s + 1];
    for (let e = this.segStart[s]; e < e1; e++) {
      const p = phiV(this.eU[e] + lam);
      if (p > 0) {
        const n = this.eN[e], T = this.eT[e], d = PH[1];
        h += n * p; hd += n * d; ht += p * T; htd += d * T;
      }
    }
    const tl = this.segLast[s];
    const S0 = gp.S;
    if (h > 0) {
      const te = ht / h, ted = (htd - te * hd) / h, eh = Math.exp(-h);
      const I = S0 * ((te - gp.tEnd) + (tl - te) * eh);
      gp.Ad += -gp.Hd * I + S0 * (ted * (1 - eh) - (tl - te) * eh * hd);
      gp.A += I;
      gp.S = S0 * eh; gp.H += h; gp.Hd += hd;
    } else {
      const I = S0 * (tl - gp.tEnd);
      gp.Ad += -gp.Hd * I;
      gp.A += I;
    }
    gp.tEnd = tl;
  }
  gridPoint(g) {
    let gp = this.grid.get(g);
    if (gp) return gp;
    gp = { lam: g * GRID_STEP, A: 0, Ad: 0, S: 1, H: 0, Hd: 0, tEnd: this.tFirst };
    for (let s = 0; s < this.nSeg; s++) this.segIntegrate(gp, s);
    this.grid.set(g, gp); this.gridInits++;
    return gp;
  }
  lnAttempt(lam) { // ln(restart + E[T_A]): cubic Hermite in λ between grid points (values + exact slopes)
    const f = lam / GRID_STEP;
    const g = Math.floor(f);
    const a = this.gridPoint(g), b = this.gridPoint(g + 1);
    const ra = this.restart + a.A, rb = this.restart + b.A;
    const fa = Math.log(ra), fb = Math.log(rb);
    const ma = a.Ad / ra * GRID_STEP, mb = b.Ad / rb * GRID_STEP;
    const t = f - g, t2 = t * t, t3 = t2 * t;
    return (2 * t3 - 3 * t2 + 1) * fa + (t3 - 2 * t2 + t) * ma + (-2 * t3 + 3 * t2) * fb + (t3 - t2) * mb;
  }
  binStats(b) { // cached mean log and relative variance of a fine bin
    if (this.fDirty[b]) {
      const n = this.fn[b], m = this.fs[b] / n;
      this.fu[b] = Math.log(m);
      let v = this.fs2[b] / n / (m * m) - 1; if (!(v > 0)) v = 0;
      this.fv[b] = v; this.fDirty[b] = 0;
    }
  }
  // H(λ) = Σ n_b Φ(u_b + λ) (+ 2nd-order bin-width correction); also dH/dλ in this.dH
  hazard(lam, exact) {
    let h = 0, dh = 0;
    const L = this.fList, n = this.fLen;
    for (let i = 0; i < n; i++) {
      const b = L[i];
      this.binStats(b);
      const v = this.fu[b] + lam;
      if (v >= V_CUT) break;
      if (exact) phiVexact(v); else phiV(v);
      const cnt = this.fn[b];
      h += cnt * (PH[0] + 0.5 * this.fv[b] * (PH[2] - PH[1]));
      dh += cnt * PH[1];
    }
    this.dH = dh;
    return h;
  }
  // ---- online (per CSV line) solve: prefix L* ----
  solveOnline() {
    this.solveCount++;
    if (this.zero > 0) { this.curL = Infinity; return; }
    if (this.inputs === 0) { this.curL = 0; return; }
    const span = this.tLast - this.tFirst;
    if (this.restart + span >= Math.exp(this.lnTarget)) { this.curL = Infinity; return; }
    const G = (lam) => this.lnAttempt(lam) + this.hazard(lam, false) - this.lnTarget;
    let lo = this.lam;
    if (!isFinite(lo)) lo = this.initialGuess() - 8;
    let glo = G(lo);
    if (glo <= 0) { // bracket from below
      if (isFinite(this.lam)) { this.setLam(this.lam); return; }
      let step = 4; let hi = lo; lo = hi - step;
      while ((glo = G(lo)) <= 0) { hi = lo; step *= 2; lo -= step; }
    }
    let step = 0.02, hi = lo + step, ghi;
    while ((ghi = G(hi)) > 0) { lo = hi; glo = ghi; step *= 2; hi = lo + step; }
    const r = brent(G, lo, hi, glo, ghi, 1e-10);
    this.setLam(r);
  }
  setLam(lam) {
    this.lam = lam; this.curL = Math.exp(lam);
    const gmin = Math.floor(lam / GRID_STEP);
    for (const g of this.grid.keys()) if (g < gmin || g > gmin + 6) this.grid.delete(g);
  }
  initialGuess() { // λ where typical inputs have x ≈ 1
    let tot = 0;
    for (let i = 0; i < this.fLen; i++) tot += this.fn[this.fList[i]];
    let acc = 0;
    for (let i = 0; i < this.fLen; i++) {
      const b = this.fList[i]; acc += this.fn[b];
      if (acc >= tot * 0.05) { this.binStats(b); return -this.fu[b]; }
    }
    return 0;
  }
  // ---- final exact solve over the full MIDI ----
  finalAttempt(lam) { // E[T_A] over all stored segments
    const gp = { lam, A: 0, Ad: 0, S: 1, H: 0, Hd: 0, tEnd: this.tFirst };
    for (let s = 0; s < this.nSeg; s++) { this.segIntegrate(gp, s); if (gp.S < 1e-300) break; }
    return gp.A;
  }
  finish() {
    // flush the trailing (partial) frame
    const ring = this.npsRing, ri = this.nextLine % ring.length;
    this.npsSum += this.batchNotes - ring[ri]; ring[ri] = this.batchNotes;
    if (this.npsSum > this.npsMax) { this.npsMax = this.npsSum; this.npsMaxT = this.nextLine / this.fps; }
    this.totalNotes += this.batchNotes; this.batchNotes = 0;
    if (this.bLen > 0) this.flushBatch();
    const res = { L: 0, lam: -Infinity, inputs: this.inputs, attempt: 0, pc: 1, reason: '' };
    if (this.inputs === 0) { res.reason = 'no notes'; return res; }
    if (this.zero > 0) { res.L = Infinity; res.reason = `${this.zero} notes are 0 seconds apart (tempo 0 / infinite BPM)`; return res; }
    const span = this.tLast - this.tFirst;
    if (this.restart + span >= Math.exp(this.lnTarget)) { res.L = Infinity; res.reason = 'one perfect run is longer than the target time'; return res; }
    // Fixed-point iteration: λ = H⁻¹(lnT − ln(restart + E[T_A](λ))). H is cheap (fine histogram),
    // E[T_A] needs a pass over the time segments. E[T_A] depends weakly on λ → converges in a few passes.
    let nEval = 0;
    const Hinv = (target, guess) => { // solve hazard(λ, exact) = target
      const f = (lam) => this.hazard(lam, true) - target;
      let lo = guess - 0.25, hi = guess + 0.25, flo = f(lo), fhi = f(hi), st = 0.5;
      while (flo <= 0) { hi = lo; fhi = flo; lo -= st; st *= 2; flo = f(lo); }
      st = 0.5;
      while (fhi > 0) { lo = hi; flo = fhi; hi += st; st *= 2; fhi = f(hi); }
      return brent(f, lo, hi, flo, fhi, 1e-14);
    };
    const attempt = (lam) => { nEval++; return this.finalAttempt(lam); };
    let lam = Hinv(this.lnTarget - Math.log(this.restart + 0.5 * span), isFinite(this.lam) ? this.lam : this.initialGuess());
    let prev = NaN, prev2 = NaN;
    for (let it = 0; it < 40; it++) {
      const nl = Hinv(this.lnTarget - Math.log(this.restart + attempt(lam)), lam);
      if (Math.abs(nl - lam) < 1e-13 * Math.max(1, Math.abs(lam))) { lam = nl; break; }
      prev2 = prev; prev = lam; lam = nl;
      // Aitken Δ² acceleration
      if (it >= 2 && isFinite(prev2)) {
        const d1 = prev - prev2, d2 = lam - prev, den = d2 - d1;
        if (Math.abs(den) > 1e-300) { const acc = lam - d2 * d2 / den; if (isFinite(acc) && Math.abs(acc - lam) < Math.abs(d2) * 10) { lam = acc; prev = NaN; prev2 = NaN; } }
      }
    }
    res.lam = lam; res.L = Math.exp(lam); res.nEval = nEval;
    res.attempt = this.finalAttempt(lam);
    res.pc = Math.exp(-this.hazard(lam, true));
    return res;
  }
}

// Brent's root finder (f(a) > 0 > f(b) or vice versa)
function brent(f, a, b, fa, fb, tol) {
  let c = a, fc = fa, d = b - a, e = d;
  for (let it = 0; it < 200; it++) {
    if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) { c = a; fc = fa; d = b - a; e = d; }
    if (Math.abs(fc) < Math.abs(fb)) { a = b; b = c; c = a; fa = fb; fb = fc; fc = fa; }
    const tol1 = 2 * 2.2e-16 * Math.abs(b) + 0.5 * tol, xm = 0.5 * (c - b);
    if (Math.abs(xm) <= tol1 || fb === 0) return b;
    if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
      const s = fb / fa;
      let p, q;
      if (a === c) { p = 2 * xm * s; q = 1 - s; }
      else { const q0 = fa / fc, r = fb / fc; p = s * (2 * xm * q0 * (q0 - r) - (b - a) * (r - 1)); q = (q0 - 1) * (r - 1) * (s - 1); }
      if (p > 0) q = -q; p = Math.abs(p);
      const m1 = 3 * xm * q - Math.abs(tol1 * q), m2 = Math.abs(e * q);
      if (2 * p < (m1 < m2 ? m1 : m2)) { e = d; d = p / q; } else { d = xm; e = d; }
    } else { d = xm; e = d; }
    a = b; fa = fb;
    b += Math.abs(d) > tol1 ? d : (xm > 0 ? tol1 : -tol1);
    fb = f(b);
  }
  return b;
}

// =====================================================================================
// MIDI reading (plain files: one buffered reader per track, positional reads)
// =====================================================================================
const EV_BUF_MIN = 4096;

class Ctx { // shared state the track readers write into
  constructor() {
    this.words = new Int32Array(4);
    this.cnt = new Uint32Array(128); // note-ons per key at the current tick (overlaps)
    this.vm = new Uint8Array(128);    // max velocity per key at the current tick
    this.vs = new Float64Array(128);  // velocity sum per key at the current tick
    this.notes = 0;          // raw note-on count
    this.spt = 0;            // seconds per tick
    this.baseTick = 0; this.baseTime = 0;     // absolute time anchor
    this.gapAcc = Infinity; this.gapTick = 0;  // time since previous onset (exact)
    this.ppq = 0; this.smpte = false;
    this.tempoMin = Infinity; this.tempoMax = -Infinity; this.tempoCount = 0;
    this.chanMask = 0xFFFF;
  }
  tempo(tick, us) {
    if (this.smpte) return;
    this.tempoCount++;
    if (us < this.tempoMin) this.tempoMin = us;
    if (us > this.tempoMax) this.tempoMax = us;
    this.baseTime += (tick - this.baseTick) * this.spt; this.baseTick = tick;
    this.gapAcc += (tick - this.gapTick) * this.spt; this.gapTick = tick;
    this.spt = us * 1e-6 / this.ppq;
  }
}

const PRIME_CTX = { words: new Int32Array(4), cnt: new Uint32Array(128), vm: new Uint8Array(128), vs: new Float64Array(128), notes: 0, tempo() {} };

class FileTrackReader {
  constructor(fd, start, length, bufSize, chanMask) {
    this.fd = fd; this.fpos = start; this.left = length;
    this.buf = Buffer.allocUnsafe(Math.max(EV_BUF_MIN, Math.min(bufSize, length + 16)));
    this.pos = 0; this.end = 0;
    this.tick = 0; this.rs = 0; this.endTick = 0;
    this.mask = chanMask === undefined ? 0xFFFF : chanMask;
    this.pend = 0; // pending relevant event: 1 = note(key), 2 = tempo, 0 = none
    this.pv = 0; this.pvel = 0;
    this.next = -1;            // prime: nothing happens at tick -1, stops at the first relevant event
    this.step(PRIME_CTX);
  }
  refill() {
    const keep = this.end - this.pos;
    if (keep > 0) this.buf.copyWithin(0, this.pos, this.end);
    this.pos = 0; this.end = keep;
    const want = Math.min(this.buf.length - keep, this.left);
    if (want > 0) {
      let got = 0;
      while (got < want) {
        const r = fs.readSync(this.fd, this.buf, keep + got, want - got, this.fpos + got);
        if (r <= 0) break;
        got += r;
      }
      this.fpos += got; this.left -= got; this.end += got;
      if (got < want) this.left = 0;
    }
  }
  skip(n) { // skip n bytes of payload
    const avail = this.end - this.pos;
    if (n <= avail) { this.pos += n; return; }
    n -= avail; this.pos = this.end;
    if (n >= this.left) { this.fpos += this.left; this.left = 0; return; }
    this.fpos += n; this.left -= n;
  }
  // Consume every event at tick this.next (notes -> ctx.words, tempo -> ctx), then parse ahead to
  // the first relevant event (note-on vel>0 / tempo) at a later tick, keep it pending, return its tick.
  step(ctx) {
    const T = this.next;
    const words = ctx.words, cnt = ctx.cnt, vm = ctx.vm, vs = ctx.vs;
    if (this.pend === 1) { const k = this.pv, v = this.pvel; words[k >> 5] |= 1 << (k & 31); cnt[k]++; vs[k] += v; if (v > vm[k]) vm[k] = v; ctx.notes++; }
    else if (this.pend === 2) ctx.tempo(T, this.pv);
    this.pend = 0;
    const buf = this.buf, mask = this.mask;
    let pos = this.pos, end = this.end, tick = this.tick, rs = this.rs, notes = 0;
    for (;;) {
      if (end - pos < 16 && this.left > 0) { this.pos = pos; this.refill(); pos = this.pos; end = this.end; }
      if (pos >= end) break;
      let b = buf[pos++], d = b & 0x7F;
      while (b & 0x80) { if (pos >= end) { pos = end; break; } b = buf[pos++]; d = d * 128 + (b & 0x7F); }
      if (pos >= end) break;
      let st = buf[pos];
      if (st & 0x80) pos++; else st = rs;
      const hi = st >> 4;
      if (hi >= 8 && hi <= 0xE) { // channel message
        if (hi === 9) {
          if (pos + 2 > end) { pos = end; break; }
          rs = st; tick += d;
          const v = buf[pos + 1];
          if (v !== 0 && ((mask >> (st & 15)) & 1)) {
            const k = buf[pos] & 127; pos += 2;
            if (tick === T) { words[k >> 5] |= 1 << (k & 31); cnt[k]++; vs[k] += v; if (v > vm[k]) vm[k] = v; notes++; continue; }
            this.pend = 1; this.pv = k; this.pvel = v;
            this.pos = pos; this.tick = tick; this.rs = rs; ctx.notes += notes; this.next = tick;
            return tick;
          }
          pos += 2; continue;
        }
        const len = (hi === 0xC || hi === 0xD) ? 1 : 2;
        if (pos + len > end) { pos = end; break; }
        rs = st; tick += d; pos += len; continue;
      }
      tick += d;
      if (st === 0xFF || st === 0xF0 || st === 0xF7) {
        let type = 0;
        if (st === 0xFF) { if (pos >= end) { pos = end; break; } type = buf[pos++]; }
        let len = 0, c = 0x80;
        while (c & 0x80) { if (pos >= end) break; c = buf[pos++]; len = len * 128 + (c & 0x7F); }
        if (c & 0x80) { pos = end; break; }
        if (st === 0xFF) {
          if (type === 0x51 && len >= 3) {
            if (pos + 3 > end) { pos = end; break; }
            const us = (buf[pos] << 16) | (buf[pos + 1] << 8) | buf[pos + 2];
            this.pos = pos; this.skip(len); pos = this.pos; end = this.end;
            if (tick === T) { ctx.tempo(T, us); continue; }
            this.pend = 2; this.pv = us;
            this.tick = tick; this.rs = rs; ctx.notes += notes; this.next = tick;
            return tick;
          }
          if (type === 0x2F) { this.left = 0; pos = end; break; } // end of track
        }
        this.pos = pos; this.skip(len); pos = this.pos; end = this.end;
        continue;
      }
      if (st < 0x80) { pos++; continue; }  // data byte without running status
      pos += (st === 0xF2 ? 2 : (st === 0xF1 || st === 0xF3) ? 1 : 0); // stray system common
    }
    // end of track
    this.pos = end; this.left = 0; this.tick = tick; this.endTick = tick; this.rs = rs;
    ctx.notes += notes; this.next = Infinity;
    return Infinity;
  }
}

function findMThd(buf) {
  for (let i = 0; i + 14 <= buf.length; i++) {
    if (buf[i] === 0x4D && buf[i + 1] === 0x54 && buf[i + 2] === 0x68 && buf[i + 3] === 0x64) return i;
  }
  return -1;
}

function openMidiFile(file, base, length) {
  const fd = fs.openSync(file, 'r');
  const fsize = fs.fstatSync(fd).size;
  base = base || 0;
  const size = base + (length === undefined ? fsize - base : Math.min(length, fsize - base)); // end offset
  const head = Buffer.alloc(Math.max(0, Math.min(size - base, 1 << 16)));
  fs.readSync(fd, head, 0, head.length, base);
  const at = findMThd(head);
  if (at < 0) { fs.closeSync(fd); throw new Error('not a MIDI file (no MThd header)'); }
  const hlen = head.readUInt32BE(at + 4);
  const fmt = head.readUInt16BE(at + 8), ntrk = head.readUInt16BE(at + 10), div = head.readUInt16BE(at + 12);
  const tracks = [];
  let p = base + at + 8 + hlen;
  const h8 = Buffer.alloc(8);
  while (p + 8 <= size) {
    fs.readSync(fd, h8, 0, 8, p);
    const id = h8.toString('latin1', 0, 4);
    let len = h8.readUInt32BE(4);
    if (id === 'MTrk') {
      if (p + 8 + len > size) len = size - p - 8;
      tracks.push({ start: p + 8, length: len });
    } else if (!/^[\x20-\x7E]{4}$/.test(id)) break; // garbage
    p += 8 + len;
  }
  return { fd, size: size - base, format: fmt, declaredTracks: ntrk, division: div, tracks };
}

// =====================================================================================
// Analysis driver: merge tracks by tick, build onsets, run hand model + accumulator
// =====================================================================================
function analyzeReaders(readers, division, P, opts) {
  const ctx = new Ctx();
  ctx.chanMask = opts.chanMask;
  if (division & 0x8000) {
    const fps = 256 - (division >> 8), tpf = division & 0xFF;
    ctx.smpte = true;
    ctx.spt = 1 / ((fps === 29 ? 29.97 : fps) * (tpf || 1));
  } else {
    ctx.ppq = division || 96;
    ctx.spt = 0.5 / ctx.ppq; // 120 BPM default
  }
  const hands = new HandModel(P);
  const acc = new Accumulator(P, opts);
  const timingOnly = !!opts.timingOnly;
  const H = P.window;
  // ---- merge all tracks by tick: calendar queue (O(1) per track-tick) + heap for far-away ticks
  const n = readers.length;
  const R = 4096, RM = R - 1;
  const head = new Int32Array(R).fill(-1), nxt = new Int32Array(Math.max(1, n)), bits = new Int32Array(R >> 5);
  const tk = new Float64Array(Math.max(1, n)), far = new Int32Array(Math.max(1, n));
  let farN = 0, cur = Infinity;
  for (let i = 0; i < n; i++) { tk[i] = readers[i].next; if (tk[i] < cur) cur = tk[i]; }
  const farPush = (i) => { // min-heap on tk
    let c = farN++;
    while (c > 0) { const p = (c - 1) >> 1; if (tk[far[p]] <= tk[i]) break; far[c] = far[p]; c = p; }
    far[c] = i;
  };
  const farPop = () => {
    const top = far[0], x = far[--farN];
    if (farN > 0) {
      const kx = tk[x]; let i = 0;
      for (;;) { let c = 2 * i + 1; if (c >= farN) break; if (c + 1 < farN && tk[far[c + 1]] < tk[far[c]]) c++; if (tk[far[c]] >= kx) break; far[i] = far[c]; i = c; }
      far[i] = x;
    }
    return top;
  };
  const toBucket = (i) => { const b = tk[i] & RM; nxt[i] = head[b]; head[b] = i; bits[b >> 5] |= 1 << (b & 31); };
  if (cur !== Infinity) for (let i = 0; i < n; i++) if (tk[i] !== Infinity) { if (tk[i] - cur < R) toBucket(i); else farPush(i); }
  let curTick = 0, onsets = 0, presses = 0, extraInputs = 0;
  let notesAtLast = 0;
  // pending onset (window needs the next onset's time)
  let pT = 0, pGap = Infinity, pN = 0, pMu0 = 1, pMu1 = 1, pNotes = 0, havePend = false;
  const words = ctx.words, cnt = ctx.cnt, vm = ctx.vm, vs = ctx.vs;
  // overlapped notes (same key, same tick) as swift re-strikes: multiplicity c, hand μ, velocity weight
  const overlaps = opts.overlaps !== false;
  const useVel = opts.velocity !== false;
  const OMEGA = new Float64Array(128);   // 1/ω(v): window multiplier for velocity v
  for (let v = 0; v < 128; v++) OMEGA[v] = useVel ? Math.pow(127 / Math.max(1, v), P.velocityPower) : 1;
  const omega = (v) => (v === (v | 0) ? OMEGA[v] : (useVel ? Math.pow(127 / Math.max(1, v), P.velocityPower) : 1));
  // pending re-strike groups: identical (multiplicity, weight) are merged; pXk = number of keys
  const pXc = new Int32Array(128), pXmu = new Float64Array(128), pXk = new Float64Array(128);
  let pXn = 0;
  const emitPending = (gapNext) => {
    let h = pGap < gapNext ? pGap : gapNext; h *= 0.5; if (h > H) h = H;
    acc.addNotes(pT, pNotes);
    acc.addInput(pT, h * pMu0);
    if (pN === 2) acc.addInput(pT, h * pMu1);
    let lastC = -1, rf = 1;
    for (let x = 0; x < pXn; x++) {
      // c strikes of one key share the window: each gets 1/c of it, re-striking every 2h/c seconds
      const c = pXc[x], hc = h / c;
      if (c !== lastC) { rf = timingOnly ? 1 : Math.exp(-hands.repCost(2 * hc)); lastC = c; }
      const n = (c - 1) * pXk[x];
      acc.addInputN(pT, hc * pXmu[x] * rf, n);
      extraInputs += n;
    }
  };

  const finalizeOnset = () => {
    const t = ctx.baseTime + (curTick - ctx.baseTick) * ctx.spt;
    const gap = ctx.gapAcc + (curTick - ctx.gapTick) * ctx.spt;
    ctx.gapAcc = 0; ctx.gapTick = curTick;
    if (havePend) emitPending(gap);
    const m = hands.onset(words, t);
    words[0] = 0; words[1] = 0; words[2] = 0; words[3] = 0;
    onsets++; presses += m;
    const rawNotes = ctx.notes - notesAtLast; notesAtLast = ctx.notes;
    // per key: loudest velocity per hand, overlapped copies
    const K = hands.K, sj = hands.splitJ;
    let vL = 0, vR = 0, e0 = 0, e1 = 0;
    pXn = 0;
    const dupes = overlaps && rawNotes > m;
    for (let j = 0; j < m; j++) {
      const k = K[j], v = vm[k];
      if (j < sj) { if (v > vL) vL = v; } else if (v > vR) vR = v;
      if (dupes) {
        const c = cnt[k];
        if (c > 1) {
          const hh = j < sj ? 0 : 1;
          // the loudest copy is the chord's own strike; the others are re-strikes (mean velocity)
          const vr = (vs[k] - v) / (c - 1);
          const w = (timingOnly ? 1 : hands.muH[hh]) * omega(vr);
          let x = 0;
          while (x < pXn && !(pXc[x] === c && pXmu[x] === w)) x++;
          if (x === pXn) { pXc[x] = c; pXmu[x] = w; pXk[x] = 1; pXn++; } else pXk[x]++;
          if (hh === 0) e0 += c - 1; else e1 += c - 1;
        }
      }
      cnt[k] = 0; vm[k] = 0; vs[k] = 0;
    }
    if (dupes) { hands.load[0] += e0; hands.load[1] += e1; } // re-strikes also tire the hand
    pT = t; pGap = gap; pN = hands.nG; pNotes = rawNotes; havePend = true;
    if (pN === 2) {
      pMu0 = (timingOnly ? 1 : hands.mu[0]) * OMEGA[vL];
      pMu1 = (timingOnly ? 1 : hands.mu[1]) * OMEGA[vR];
    } else {
      pMu0 = (timingOnly ? 1 : hands.mu[0]) * OMEGA[vL > vR ? vL : vR]; pMu1 = 1;
    }
  };

  while (cur !== Infinity) {
    // every track with an event at tick `cur`
    const b = cur & RM;
    let i = head[b]; head[b] = -1; bits[b >> 5] &= ~(1 << (b & 31));
    while (i !== -1) {
      const ni = nxt[i];
      const nt = readers[i].step(ctx);
      if (nt !== Infinity) { tk[i] = nt; if (nt - cur < R) toBucket(i); else farPush(i); }
      i = ni;
    }
    curTick = cur;
    if ((words[0] | words[1] | words[2] | words[3]) !== 0) finalizeOnset();
    // next non-empty tick
    while (farN > 0 && tk[far[0]] - cur < R) toBucket(farPop());
    let idx = (b + 1) & RM, found = -1;
    for (let scanned = 0; scanned < R;) {
      const w = idx >> 5, off = idx & 31;
      const word = bits[w] & (-1 << off);
      if (word !== 0) { found = (w << 5) + (31 - Math.clz32(word & -word)); break; }
      scanned += 32 - off; idx = ((w + 1) << 5) & RM;
    }
    if (found >= 0) cur += (found - b) & RM;
    else if (farN > 0) { cur = tk[far[0]]; while (farN > 0 && tk[far[0]] - cur < R) toBucket(farPop()); }
    else cur = Infinity;
  }
  if (havePend) emitPending(Infinity);
  let endTick = 0;
  for (const rd of readers) if (rd.endTick > endTick) endTick = rd.endTick;
  if (curTick > endTick) endTick = curTick;
  if (!(endTick >= 0)) endTick = 0;
  const endTime = ctx.baseTime + (endTick - ctx.baseTick) * ctx.spt;
  return { ctx, hands, acc, onsets, presses, notes: ctx.notes, endTick, endTime, extraInputs, overlaps };
}

// =====================================================================================
// CSV export
// =====================================================================================
function writeCsv(file, acc, endTime, finalL, opts) {
  const fps = acc.fps;
  // lines k = 0..K: state of all inputs with t < k/fps
  let K = Math.ceil(endTime * fps - 1e-9);
  const needK = Math.floor(acc.tLast * fps + 1e-12) + 1;
  if (acc.inputs > 0 && K < needK) K = needK;
  if (K < 0) K = 0;
  // make sure all lines exist (advance with the final state)
  if (acc.nLines < K + 1) {
    // lines after the last input: everything is included
    const old = acc.lineVal, oldN = acc.lineNotes;
    const nv = new Float64Array(K + 1), nn = new Float64Array(K + 1);
    nv.set(old.subarray(0, Math.min(acc.nLines, K + 1))); nn.set(oldN.subarray(0, Math.min(acc.nLines, K + 1)));
    for (let k = acc.nLines; k <= K; k++) { nv[k] = finalL; nn[k] = acc.totalNotes; }
    acc.lineVal = nv; acc.lineNotes = nn; acc.nLines = K + 1;
  }
  const V = acc.lineVal, N = acc.lineNotes;
  // lines whose prefix already contains every input use the exact final value
  const kAll = Math.floor(acc.tLast * fps + 1e-12) + 1;
  for (let k = kAll; k <= K; k++) V[k] = finalL;
  // the true prefix precision never decreases and never exceeds the final value:
  // project the (≈1e-8 accurate) online values onto that set
  let run = 0;
  for (let k = 0; k <= K; k++) { let v = V[k]; if (v < run) v = run; if (v > finalL) v = finalL; V[k] = v; run = v; }
  const dec = opts.decimals;
  const tpl = opts.format;
  const parts = compileTemplate(tpl);
  const per = Math.max(1, Math.round(fps));
  const out = [];
  let mps = 0, mnps = 0;
  const fmt = (x) => (x === Infinity ? 'inf' : x.toFixed(dec));
  let chunk = '';
  const fd = fs.openSync(file, 'w');
  for (let k = 0; k <= K; k++) {
    const p = V[k];
    const prev = k >= per ? V[k - per] : 0;
    let ps = p - prev; if (!(ps >= 0)) ps = p === Infinity ? Infinity : 0;
    if (ps > mps) mps = ps;
    const nps = N[k] - (k >= per ? N[k - per] : 0);
    if (nps > mnps) mnps = nps;
    let line = '';
    for (const pt of parts) {
      switch (pt) {
        case 0: line += fmt(p); break;
        case 1: line += fmt(ps); break;
        case 2: line += fmt(mps); break;
        case 3: line += (k / fps).toFixed(3); break;
        case 4: line += k; break;
        case 5: line += N[k]; break;
        case 6: line += nps; break;
        case 7: line += mnps; break;
        default: line += pt;
      }
    }
    chunk += k < K ? line + '\n' : line;
    if (chunk.length > 1 << 20) { fs.writeSync(fd, chunk); chunk = ''; }
  }
  if (chunk.length) fs.writeSync(fd, chunk);
  fs.closeSync(fd);
  return { lines: K + 1, maxPS: mps };
}
function compileTemplate(t) {
  const map = { p: 0, ps: 1, mps: 2, t: 3, f: 4, n: 5, nps: 6, mnps: 7 };
  const parts = [];
  const re = /\{(p|ps|mps|t|f|n|nps|mnps)\}/g;
  let last = 0, m;
  while ((m = re.exec(t))) {
    if (m.index > last) parts.push(t.slice(last, m.index));
    parts.push(map[m[1]]);
    last = re.lastIndex;
  }
  if (last < t.length) parts.push(t.slice(last));
  return parts.map((x) => (typeof x === 'string' ? x.replace(/\\n/g, '\n').replace(/\\t/g, '\t') : x));
}

// =====================================================================================
// Pretty output
// =====================================================================================
function fmtNum(x) {
  if (x === Infinity) return '∞';
  if (!isFinite(x)) return String(x);
  if (Math.abs(x) >= 1e15) return x.toExponential(4);
  return x.toLocaleString('en-US', { maximumFractionDigits: 0 });
}
function fmtPrec(x, dec = 3) {
  if (x === Infinity) return '∞';
  if (x >= 1e15) return x.toExponential(6);
  return x.toLocaleString('en-US', { minimumFractionDigits: dec, maximumFractionDigits: dec });
}
function fmtTime(s) {
  if (!isFinite(s)) return String(s);
  if (s < 1e-3) return (s * 1e6).toPrecision(4) + ' µs';
  if (s < 60) return s.toFixed(3) + ' s';
  const h = Math.floor(s / 3600), m = Math.floor((s % 3600) / 60), sec = s % 60;
  return (h ? h + ':' + String(m).padStart(2, '0') : String(m)) + ':' + sec.toFixed(3).padStart(6, '0');
}
function fmtSigma(L) {
  if (!(L > 0)) return '—';
  if (L === Infinity) return '0 (impossible)';
  const s = 1 / L;
  if (s >= 1) return s.toPrecision(3) + ' s';
  if (s >= 1e-3) return (s * 1e3).toPrecision(3) + ' ms';
  if (s >= 1e-6) return (s * 1e6).toPrecision(3) + ' µs';
  if (s >= 1e-9) return (s * 1e9).toPrecision(3) + ' ns';
  if (s >= 1e-12) return (s * 1e12).toPrecision(3) + ' ps';
  return s.toExponential(2) + ' s';
}
function bpm(us) { return us > 0 ? 60e6 / us : Infinity; }
function fmtBpm(us) { const b = bpm(us); return b === Infinity ? '∞' : String(+b.toPrecision(b >= 1e6 ? 8 : 5)); }

// =====================================================================================
// Decompression building blocks (pure JS): buffered input, LZMA/LZMA2, XZ, BZip2, filters
// =====================================================================================
// Buffered sequential input over a file region (or an in-memory buffer). Keeps >= `pad`
// readable bytes past `end` (zero filled at EOF) so decoders can read without bounds checks.
class InStream {
  constructor(fd, start, length, bufSize) {
    this.fd = fd; this.fpos = start; this.left = length === undefined ? Infinity : length;
    this.buf = new Uint8Array((bufSize || (1 << 20)) + 128);
    this.pos = 0; this.end = 0; this.consumedBefore = 0; this.eof = false;
  }
  static fromBuffer(u8) {
    const s = new InStream(-1, 0, 0, 16);
    s.buf = new Uint8Array(u8.length + 128); s.buf.set(u8); s.end = u8.length; s.left = 0; s.eof = true;
    return s;
  }
  tell() { return this.consumedBefore + this.pos; }
  // make at least n bytes available after pos (unless EOF); returns available count
  fill(n) {
    if (this.end - this.pos >= n || this.eof) return this.end - this.pos;
    const keep = this.end - this.pos;
    if (this.pos > 0) { this.buf.copyWithin(0, this.pos, this.end); this.consumedBefore += this.pos; this.pos = 0; this.end = keep; }
    while (this.end - this.pos < n && !this.eof) {
      const room = this.buf.length - 128 - this.end;
      const want = Math.min(room, this.left);
      if (want <= 0) { if (this.left <= 0) { this.eof = true; break; } else break; }
      const got = fs.readSync(this.fd, this.buf, this.end, want, this.fpos);
      if (got <= 0) { this.eof = true; break; }
      this.fpos += got; this.left -= got; this.end += got;
      if (this.left <= 0) this.eof = true;
    }
    if (this.eof) this.buf.fill(0, this.end, Math.min(this.buf.length, this.end + 128));
    return this.end - this.pos;
  }
  byte() { if (this.pos >= this.end && this.fill(1) < 1) throw new Error('unexpected end of compressed data'); return this.buf[this.pos++]; }
  bytes(n) {
    const out = new Uint8Array(n);
    let o = 0;
    while (o < n) {
      if (this.pos >= this.end && this.fill(1) < 1) throw new Error('unexpected end of compressed data');
      const k = Math.min(n - o, this.end - this.pos);
      out.set(this.buf.subarray(this.pos, this.pos + k), o); this.pos += k; o += k;
    }
    return out;
  }
  skip(n) {
    while (n > 0) {
      if (this.pos >= this.end) {
        // jump in the file directly
        if (this.fd >= 0 && n > this.buf.length) { const k = Math.min(n, this.left); this.consumedBefore += this.pos + k; this.pos = 0; this.end = 0; this.fpos += k; this.left -= k; n -= k; if (this.left <= 0) this.eof = true; if (k === 0) throw new Error('unexpected end of compressed data'); continue; }
        if (this.fill(1) < 1) throw new Error('unexpected end of compressed data');
      }
      const k = Math.min(n, this.end - this.pos); this.pos += k; n -= k;
    }
  }
  // remaining payload (for "copy" coders)
  pipeTo(sink, n) {
    while (n > 0) {
      if (this.pos >= this.end && this.fill(1) < 1) throw new Error('unexpected end of compressed data');
      const k = Math.min(n, this.end - this.pos);
      sink(this.buf, this.pos, this.pos + k); this.pos += k; n -= k;
    }
  }
}

// ---- LZMA / LZMA2 --------------------------------------------------------------------
const LZ_IS_MATCH = 0, LZ_IS_REP = 192, LZ_IS_REP_G0 = 204, LZ_IS_REP_G1 = 216, LZ_IS_REP_G2 = 228,
  LZ_IS_REP0_LONG = 240, LZ_POS_SLOT = 432, LZ_SPEC_POS = 688, LZ_ALIGN = 803, LZ_LEN = 819, LZ_REP_LEN = 1333,
  LZ_LITERAL = 1847;

class LzmaDecoder {
  // sink(buf, start, end) receives the decoded bytes in order
  constructor(dictSize, sink, maxOut) {
    let ws = Math.max(1 << 16, dictSize >>> 0);
    if (maxOut !== undefined && maxOut >= 0 && maxOut + 1 < ws) ws = Math.max(1 << 16, Math.ceil(maxOut) + 1);
    this.ws = ws; this.win = new Uint8Array(ws);
    this.wp = 0; this.fp = 0; this.filled = 0; // write pos, flushed pos, bytes since dict reset (capped)
    this.total = 0;                             // bytes decoded (for lp/pb position bits)
    this.sink = sink;
    this.pr = null; this.lc = 3; this.lp = 0; this.pb = 2;
    this.state = 0; this.rep0 = 0; this.rep1 = 0; this.rep2 = 0; this.rep3 = 0;
    this.range = 0; this.code = 0;
  }
  setProps(b) {
    if (b >= 225) throw new Error('bad LZMA properties');
    const lc = b % 9; b = (b / 9) | 0; const lp = b % 5; const pb = (b / 5) | 0;
    this.lc = lc; this.lp = lp; this.pb = pb;
    const n = LZ_LITERAL + (0x300 << (lc + lp));
    if (!this.pr || this.pr.length < n) this.pr = new Uint16Array(n);
  }
  resetState() {
    this.pr.fill(1024, 0, LZ_LITERAL + (0x300 << (this.lc + this.lp)));
    this.state = 0; this.rep0 = this.rep1 = this.rep2 = this.rep3 = 0;
  }
  resetDict() { this.filled = 0; }
  initRc(src) {
    src.fill(5);
    if (src.byte() !== 0) throw new Error('corrupt LZMA stream');
    let c = 0; for (let i = 0; i < 4; i++) c = c * 256 + src.byte();
    this.code = c; this.range = 0xFFFFFFFF;
  }
  flush() {
    if (this.wp > this.fp) { this.sink(this.win, this.fp, this.wp); this.fp = this.wp; }
    if (this.wp === this.ws) { this.wp = 0; this.fp = 0; }
  }
  putUncompressed(src, n) { // LZMA2 stored chunk
    while (n > 0) {
      if (src.pos >= src.end && src.fill(1) < 1) throw new Error('unexpected end of compressed data');
      let k = Math.min(n, src.end - src.pos, this.ws - this.wp);
      this.win.set(src.buf.subarray(src.pos, src.pos + k), this.wp);
      src.pos += k; this.wp += k; n -= k; this.total += k;
      this.filled = Math.min(this.ws, this.filled + k);
      if (this.wp === this.ws || this.wp - this.fp >= (1 << 20)) this.flush();
    }
  }
  // Decode until `outLimit` bytes are produced (Infinity: until end marker). Returns true on end marker.
  decode(src, outLimit, allowEnd) {
    const pr = this.pr, win = this.win, ws = this.ws;
    const lc = this.lc, lpMask = (1 << this.lp) - 1, pbMask = (1 << this.pb) - 1, lcShift = 8 - lc;
    let range = this.range, code = this.code, state = this.state;
    let rep0 = this.rep0, rep1 = this.rep1, rep2 = this.rep2, rep3 = this.rep3;
    let wp = this.wp, filled = this.filled, total = this.total;
    let ib = src.buf, ip = src.pos, ie = src.end;
    let left = outLimit, ended = false;
    let p = 0, bound = 0, bit = 0;
    while (left > 0) {
      if (ie - ip < 64) { src.pos = ip; src.fill(64); ib = src.buf; ip = src.pos; ie = src.end; if (ip > ie + 32) throw new Error('unexpected end of compressed data'); }
      if (wp === ws || wp - this.fp >= (1 << 20)) { this.wp = wp; this.flush(); wp = this.wp; }
      const posState = total & pbMask;
      // isMatch
      let i = LZ_IS_MATCH + (state << 4) + posState;
      p = pr[i]; bound = (range >>> 11) * p;
      if (code < bound) {
        range = bound; pr[i] = p + ((2048 - p) >> 5);
        if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
        // ---- literal
        const prev = filled > 0 ? win[wp > 0 ? wp - 1 : ws - 1] : 0;
        const base = LZ_LITERAL + 0x300 * (((total & lpMask) << lc) + (prev >> lcShift));
        let sym = 1;
        if (state >= 7) {
          let mi = wp - rep0 - 1; if (mi < 0) mi += ws;
          let mb = win[mi];
          do {
            const mbit = (mb >> 7) & 1; mb <<= 1;
            i = base + ((1 + mbit) << 8) + sym;
            p = pr[i]; bound = (range >>> 11) * p;
            if (code < bound) { range = bound; pr[i] = p + ((2048 - p) >> 5); bit = 0; } else { range -= bound; code -= bound; pr[i] = p - (p >> 5); bit = 1; }
            if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
            sym = (sym << 1) | bit;
            if (mbit !== bit) break;
          } while (sym < 0x100);
        }
        while (sym < 0x100) {
          i = base + sym;
          p = pr[i]; bound = (range >>> 11) * p;
          if (code < bound) { range = bound; pr[i] = p + ((2048 - p) >> 5); sym <<= 1; } else { range -= bound; code -= bound; pr[i] = p - (p >> 5); sym = (sym << 1) | 1; }
          if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
        }
        win[wp++] = sym & 0xFF; total++; left--; if (filled < ws) filled++;
        state = state < 4 ? 0 : (state < 10 ? state - 3 : state - 6);
        continue;
      }
      range -= bound; code -= bound; pr[i] = p - (p >> 5);
      if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
      let lenBase;
      // isRep
      i = LZ_IS_REP + state;
      p = pr[i]; bound = (range >>> 11) * p;
      if (code < bound) {
        range = bound; pr[i] = p + ((2048 - p) >> 5);
        if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
        rep3 = rep2; rep2 = rep1; rep1 = rep0;
        lenBase = LZ_LEN;
        state = state < 7 ? 7 : 10;
        rep0 = -1; // decoded after the length
      } else {
        range -= bound; code -= bound; pr[i] = p - (p >> 5);
        if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
        if (filled === 0) throw new Error('corrupt LZMA stream');
        i = LZ_IS_REP_G0 + state;
        p = pr[i]; bound = (range >>> 11) * p;
        if (code < bound) {
          range = bound; pr[i] = p + ((2048 - p) >> 5);
          if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
          i = LZ_IS_REP0_LONG + (state << 4) + posState;
          p = pr[i]; bound = (range >>> 11) * p;
          if (code < bound) { // short rep: one byte at rep0
            range = bound; pr[i] = p + ((2048 - p) >> 5);
            if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
            state = state < 7 ? 9 : 11;
            let si = wp - rep0 - 1; if (si < 0) si += ws;
            win[wp++] = win[si]; total++; left--; if (filled < ws) filled++;
            continue;
          }
          range -= bound; code -= bound; pr[i] = p - (p >> 5);
          if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
        } else {
          range -= bound; code -= bound; pr[i] = p - (p >> 5);
          if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
          let dist;
          i = LZ_IS_REP_G1 + state;
          p = pr[i]; bound = (range >>> 11) * p;
          if (code < bound) {
            range = bound; pr[i] = p + ((2048 - p) >> 5);
            if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
            dist = rep1;
          } else {
            range -= bound; code -= bound; pr[i] = p - (p >> 5);
            if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
            i = LZ_IS_REP_G2 + state;
            p = pr[i]; bound = (range >>> 11) * p;
            if (code < bound) {
              range = bound; pr[i] = p + ((2048 - p) >> 5);
              if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
              dist = rep2;
            } else {
              range -= bound; code -= bound; pr[i] = p - (p >> 5);
              if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
              dist = rep3; rep3 = rep2;
            }
            rep2 = rep1;
          }
          rep1 = rep0; rep0 = dist;
        }
        lenBase = LZ_REP_LEN;
        state = state < 7 ? 8 : 11;
      }
      // ---- length (lenBase: choice, choice2, low[16][8], mid[16][8], high[256])
      let len, tb, nb;
      i = lenBase;
      p = pr[i]; bound = (range >>> 11) * p;
      if (code < bound) {
        range = bound; pr[i] = p + ((2048 - p) >> 5);
        if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
        tb = lenBase + 2 + (posState << 3); nb = 3; len = 0;
      } else {
        range -= bound; code -= bound; pr[i] = p - (p >> 5);
        if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
        i = lenBase + 1;
        p = pr[i]; bound = (range >>> 11) * p;
        if (code < bound) {
          range = bound; pr[i] = p + ((2048 - p) >> 5);
          if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
          tb = lenBase + 130 + (posState << 3); nb = 3; len = 8;
        } else {
          range -= bound; code -= bound; pr[i] = p - (p >> 5);
          if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
          tb = lenBase + 258; nb = 8; len = 16;
        }
      }
      let m = 1;
      for (let k = 0; k < nb; k++) {
        i = tb + m;
        p = pr[i]; bound = (range >>> 11) * p;
        if (code < bound) { range = bound; pr[i] = p + ((2048 - p) >> 5); m <<= 1; } else { range -= bound; code -= bound; pr[i] = p - (p >> 5); m = (m << 1) | 1; }
        if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
      }
      len += m - (1 << nb);
      // ---- distance for plain matches
      if (rep0 === -1) {
        const lenState = len < 4 ? len : 3;
        tb = LZ_POS_SLOT + (lenState << 6); m = 1;
        for (let k = 0; k < 6; k++) {
          i = tb + m;
          p = pr[i]; bound = (range >>> 11) * p;
          if (code < bound) { range = bound; pr[i] = p + ((2048 - p) >> 5); m <<= 1; } else { range -= bound; code -= bound; pr[i] = p - (p >> 5); m = (m << 1) | 1; }
          if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
        }
        const slot = m - 64;
        let dist;
        if (slot < 4) dist = slot;
        else {
          const ndb = (slot >> 1) - 1;
          dist = (2 | (slot & 1)) * Math.pow(2, ndb);
          if (slot < 14) {
            tb = LZ_SPEC_POS + dist - slot - 1; m = 1; let r = 0;
            for (let k = 0; k < ndb; k++) {
              i = tb + m;
              p = pr[i]; bound = (range >>> 11) * p;
              if (code < bound) { range = bound; pr[i] = p + ((2048 - p) >> 5); m <<= 1; } else { range -= bound; code -= bound; pr[i] = p - (p >> 5); m = (m << 1) | 1; r |= 1 << k; }
              if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
            }
            dist += r;
          } else {
            let r = 0;
            for (let k = ndb - 4; k > 0; k--) {
              range = range >>> 1;
              if (code >= range) { code -= range; r = r * 2 + 1; } else r *= 2;
              if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
            }
            dist += r * 16;
            m = 1; r = 0;
            for (let k = 0; k < 4; k++) {
              i = LZ_ALIGN + m;
              p = pr[i]; bound = (range >>> 11) * p;
              if (code < bound) { range = bound; pr[i] = p + ((2048 - p) >> 5); m <<= 1; } else { range -= bound; code -= bound; pr[i] = p - (p >> 5); m = (m << 1) | 1; r |= 1 << k; }
              if (range < 16777216) { range *= 256; code = code * 256 + ib[ip++]; }
            }
            dist += r;
          }
        }
        if (dist === 0xFFFFFFFF) { ended = true; break; }   // end marker
        rep0 = dist;
        if (rep0 >= filled) throw new Error('corrupt LZMA stream (distance)');
      }
      // ---- copy len+2 bytes from distance rep0+1
      len += 2;
      if (len > left) { if (outLimit === Infinity) throw new Error('corrupt LZMA stream'); len = left; }
      let si = wp - rep0 - 1; if (si < 0) si += ws;
      left -= len; total += len;
      filled = filled + len < ws ? filled + len : ws;
      while (len > 0) {
        if (wp === ws) { this.wp = wp; this.flush(); wp = this.wp; }
        let k = len;
        if (ws - wp < k) k = ws - wp;
        if (ws - si < k) k = ws - si;
        if (rep0 + 1 >= k && k > 32) win.copyWithin(wp, si, si + k);
        else for (let q = 0; q < k; q++) win[wp + q] = win[si + q];
        wp += k; si += k; if (si === ws) si = 0; len -= k;
      }
    }
    src.pos = ip;
    this.range = range; this.code = code; this.state = state;
    this.rep0 = rep0; this.rep1 = rep1; this.rep2 = rep2; this.rep3 = rep3;
    this.wp = wp; this.filled = filled; this.total = total;
    this.flush();
    if (ended && !allowEnd) throw new Error('unexpected LZMA end marker');
    return ended;
  }
}

function lzma2DictSize(b) {
  if (b > 40) throw new Error('bad LZMA2 dictionary size');
  if (b === 40) return 0xFFFFFFFF;
  return (2 | (b & 1)) * Math.pow(2, (b >> 1) + 11);
}

// LZMA2 stream (xz / 7z): decode until the end-of-stream chunk
function decodeLzma2(src, dictSize, sink, maxOut) {
  const d = new LzmaDecoder(dictSize, sink, maxOut);
  let needProps = true;
  for (;;) {
    const c = src.byte();
    if (c === 0) break;
    if (c === 1 || c === 2) {
      const n = ((src.byte() << 8) | src.byte()) + 1;
      if (c === 1) d.resetDict();
      d.putUncompressed(src, n);
      d.flush();
      continue;
    }
    if (c < 0x80) throw new Error('corrupt LZMA2 stream');
    const un = ((c & 0x1F) * 65536) + (src.byte() << 8) + src.byte() + 1;
    const pk = ((src.byte() << 8) | src.byte()) + 1;
    const reset = (c >> 5) & 3;
    if (reset >= 2) { d.setProps(src.byte()); needProps = false; }
    else if (needProps) throw new Error('corrupt LZMA2 stream (no properties)');
    if (reset === 3) d.resetDict();
    if (reset >= 1) d.resetState();
    const start = src.tell();
    d.initRc(src);
    d.decode(src, un, false);
    const used = src.tell() - start;
    if (used > pk) throw new Error('corrupt LZMA2 stream (chunk size)');
    if (used < pk) src.skip(pk - used);
  }
  d.flush();
}

// LZMA (7z coder / .lzma / zip method 14): props byte + dict size given
function decodeLzma1(src, propsByte, dictSize, outSize, sink) {
  const d = new LzmaDecoder(dictSize, sink, outSize >= 0 ? outSize : undefined);
  d.setProps(propsByte); d.resetState(); d.initRc(src);
  d.decode(src, outSize >= 0 ? outSize : Infinity, true);
  d.flush();
}

// ---- simple filters (applied to decoded output) ------------------------------------
function deltaFilter(dist, sink) {
  const hist = new Uint8Array(256); let hp = 0;
  const out = new Uint8Array(1 << 16);
  return (buf, a, b) => {
    while (a < b) {
      const n = Math.min(b - a, out.length);
      for (let i = 0; i < n; i++) { const v = (buf[a + i] + hist[(hp - dist) & 255]) & 255; hist[hp & 255] = v; hp++; out[i] = v; }
      sink(out, 0, n); a += n;
    }
  };
}

// ---- XZ container ----------------------------------------------------------------------
function readVli(src) {
  let v = 0, mul = 1;
  for (let i = 0; i < 9; i++) { const b = src.byte(); v += (b & 0x7F) * mul; if (!(b & 0x80)) return v; mul *= 128; }
  throw new Error('corrupt xz (VLI)');
}
function decodeXz(src, sink) {
  let streams = 0;
  for (;;) {
    // skip stream padding (zeros) between/after streams
    if (streams > 0) {
      let z = 0;
      while (src.fill(1) >= 1 && src.buf[src.pos] === 0) { src.pos++; z++; }
      if (src.fill(1) < 1) break;
    }
    const h = src.bytes(12);
    if (!(h[0] === 0xFD && h[1] === 0x37 && h[2] === 0x7A && h[3] === 0x58 && h[4] === 0x5A && h[5] === 0)) {
      if (streams > 0) break;
      throw new Error('not an xz stream');
    }
    const ct = h[7] & 0x0F;
    const checkSize = ct === 0 ? 0 : (4 << Math.floor((ct - 1) / 3));
    streams++;
    for (;;) {
      const blockStart = src.tell();
      const b0 = src.byte();
      if (b0 === 0) { // index
        const n = readVli(src);
        for (let i = 0; i < n; i++) { readVli(src); readVli(src); }
        const used = src.tell() - blockStart;
        src.skip(((4 - (used % 4)) % 4) + 4); // padding + CRC32
        src.skip(12);                          // stream footer
        break;
      }
      const hsize = (b0 + 1) * 4;
      const hdr = src.bytes(hsize - 1);
      const flags = hdr[0];
      const nf = (flags & 3) + 1;
      let q = 1;
      const vli = () => { let v = 0, mul = 1; for (let i = 0; i < 9; i++) { const b = hdr[q++]; v += (b & 0x7F) * mul; if (!(b & 0x80)) return v; mul *= 128; } throw new Error('corrupt xz block header'); };
      if (flags & 0x40) vli();
      if (flags & 0x80) vli();
      const filters = [];
      for (let i = 0; i < nf; i++) { const id = vli(); const ps = vli(); filters.push({ id, props: hdr.subarray(q, q + ps) }); q += ps; }
      const last = filters[filters.length - 1];
      if (last.id !== 0x21) throw new Error('unsupported xz filter chain');
      let out = sink;
      for (let i = filters.length - 2; i >= 0; i--) {
        const f = filters[i];
        if (f.id === 0x03) out = deltaFilter(f.props[0] + 1, out);
        else throw new Error(`unsupported xz filter 0x${f.id.toString(16)} (try 7-Zip)`);
      }
      const dataStart = src.tell();
      decodeLzma2(src, lzma2DictSize(last.props[0]), out);
      const csize = src.tell() - dataStart;
      src.skip(((4 - (csize % 4)) % 4) + checkSize);
    }
  }
}

// ---- BZip2 -------------------------------------------------------------------------------
function decodeBzip2(src, sink) {
  let bitBuf = 0, bitCnt = 0;
  const bits = (n) => {
    while (bitCnt < n) {
      if (src.pos >= src.end && src.fill(1) < 1) throw new Error('unexpected end of bzip2 data');
      bitBuf = ((bitBuf << 8) | src.buf[src.pos++]) >>> 0; bitCnt += 8;
      if (bitCnt > 24) { // keep at most 32 bits meaningful
        const v = (bitBuf >>> (bitCnt - n)) & ((1 << n) - 1);
        bitCnt -= n; bitBuf &= (1 << bitCnt) - 1;
        return v;
      }
    }
    const v = (bitBuf >>> (bitCnt - n)) & ((1 << n) - 1);
    bitCnt -= n; bitBuf &= bitCnt === 32 ? 0xFFFFFFFF : (1 << bitCnt) - 1;
    return v;
  };
  const out = new Uint8Array(1 << 16); let op = 0;
  const flushOut = () => { if (op) { sink(out, 0, op); op = 0; } };
  let tt = null;
  const MAXLEN = 20;
  let streamCount = 0;
  for (;;) {
    // stream header "BZh1".."BZh9" (multi-stream files allowed)
    if (streamCount > 0) {
      bitCnt = 0; bitBuf = 0;
      if (src.fill(4) < 4) break;
      if (!(src.buf[src.pos] === 0x42 && src.buf[src.pos + 1] === 0x5A && src.buf[src.pos + 2] === 0x68)) break;
    }
    if (bits(8) !== 0x42 || bits(8) !== 0x5A || bits(8) !== 0x68) throw new Error('not a bzip2 stream');
    const level = bits(8) - 0x30;
    if (level < 1 || level > 9) throw new Error('bad bzip2 block size');
    const bsMax = level * 100000;
    if (!tt || tt.length < bsMax) tt = new Uint32Array(bsMax);
    streamCount++;
    for (;;) {
      const m1 = bits(24), m2 = bits(24);
      if (m1 === 0x177245 && m2 === 0x385090) { bits(16); bits(16); break; } // end of stream + CRC
      if (m1 !== 0x314159 || m2 !== 0x265359) throw new Error('corrupt bzip2 block');
      bits(16); bits(16); // block CRC
      if (bits(1)) throw new Error('randomised bzip2 blocks are not supported');
      const origPtr = bits(24);
      // symbol map
      const used = bits(16);
      const seqToUnseq = new Uint8Array(256); let nInUse = 0;
      for (let i = 0; i < 16; i++) if (used & (0x8000 >> i)) { const u = bits(16); for (let j = 0; j < 16; j++) if (u & (0x8000 >> j)) seqToUnseq[nInUse++] = i * 16 + j; }
      if (nInUse === 0) throw new Error('corrupt bzip2 block');
      const alpha = nInUse + 2;
      const nGroups = bits(3);
      if (nGroups < 2 || nGroups > 6) throw new Error('corrupt bzip2 block');
      const nSel = bits(15);
      if (nSel < 1) throw new Error('corrupt bzip2 block');
      const mtfG = [0, 1, 2, 3, 4, 5];
      const selectors = new Uint8Array(nSel);
      for (let i = 0; i < nSel; i++) {
        let j = 0; while (bits(1)) { j++; if (j >= nGroups) throw new Error('corrupt bzip2 block'); }
        const v = mtfG[j]; for (let k = j; k > 0; k--) mtfG[k] = mtfG[k - 1]; mtfG[0] = v; selectors[i] = v;
      }
      // Huffman tables
      const limit = [], base = [], perm = [], minLens = [];
      for (let t = 0; t < nGroups; t++) {
        const len = new Uint8Array(alpha);
        let cur = bits(5);
        for (let i = 0; i < alpha; i++) {
          for (;;) { if (cur < 1 || cur > MAXLEN) throw new Error('corrupt bzip2 block'); if (!bits(1)) break; cur += bits(1) ? -1 : 1; }
          len[i] = cur;
        }
        let mn = 32, mx = 0;
        for (let i = 0; i < alpha; i++) { if (len[i] > mx) mx = len[i]; if (len[i] < mn) mn = len[i]; }
        const pm = new Int32Array(alpha), bs = new Int32Array(MAXLEN + 2), lim = new Int32Array(MAXLEN + 2).fill(-1);
        let pp = 0;
        for (let i = mn; i <= mx; i++) for (let j = 0; j < alpha; j++) if (len[j] === i) pm[pp++] = j;
        for (let i = 0; i < alpha; i++) bs[len[i] + 1]++;
        for (let i = 1; i < MAXLEN + 2; i++) bs[i] += bs[i - 1];
        let vec = 0;
        for (let i = mn; i <= mx; i++) { vec += bs[i + 1] - bs[i]; lim[i] = vec - 1; vec <<= 1; }
        for (let i = mn + 1; i <= mx; i++) bs[i] = ((lim[i - 1] + 1) << 1) - bs[i];
        limit.push(lim); base.push(bs); perm.push(pm); minLens.push(mn);
      }
      // MTF / RLE2 decode
      const EOB = nInUse + 1;
      const mtf = new Uint8Array(256); for (let i = 0; i < 256; i++) mtf[i] = i;
      const unzftab = new Int32Array(256);
      let nblock = 0, groupPos = 0, groupNo = -1, lim = null, bs = null, pm = null, mn = 0;
      let es = 0, N = 1;
      for (;;) {
        if (groupPos === 0) {
          groupNo++; if (groupNo >= nSel) throw new Error('corrupt bzip2 block');
          groupPos = 50; const g = selectors[groupNo]; lim = limit[g]; bs = base[g]; pm = perm[g]; mn = minLens[g];
        }
        groupPos--;
        let zn = mn, zvec = bits(zn);
        while (zn <= MAXLEN && zvec > lim[zn]) { zn++; zvec = (zvec << 1) | bits(1); }
        if (zn > MAXLEN) throw new Error('corrupt bzip2 block');
        const sym = pm[zvec - bs[zn]];
        if (sym <= 1) { es += (sym + 1) * N; N <<= 1; if (es > bsMax) throw new Error('corrupt bzip2 block'); continue; }
        if (es > 0) {
          const uc = seqToUnseq[mtf[0]];
          if (nblock + es > bsMax) throw new Error('corrupt bzip2 block');
          unzftab[uc] += es;
          for (let k = 0; k < es; k++) tt[nblock++] = uc;
          es = 0; N = 1;
        }
        if (sym === EOB) break;
        const idx = sym - 1;
        const v = mtf[idx]; mtf.copyWithin(1, 0, idx); mtf[0] = v;
        const uc = seqToUnseq[v];
        if (nblock >= bsMax) throw new Error('corrupt bzip2 block');
        unzftab[uc]++; tt[nblock++] = uc;
      }
      if (origPtr >= nblock) throw new Error('corrupt bzip2 block');
      // inverse BWT
      const cftab = new Int32Array(257);
      for (let i = 0; i < 256; i++) cftab[i + 1] = cftab[i] + unzftab[i];
      for (let i = 0; i < nblock; i++) { const uc = tt[i] & 0xFF; tt[cftab[uc]] |= (i << 8); cftab[uc]++; }
      let tPos = tt[origPtr] >>> 8;
      // output with RLE1 (4 equal bytes + count)
      let last = -1, run = 0;
      for (let i = 0; i < nblock; i++) {
        tPos = tt[tPos]; const ch = tPos & 0xFF; tPos >>>= 8;
        if (run === 4) {
          for (let k = 0; k < ch; k++) { out[op++] = last; if (op === out.length) flushOut(); }
          run = 0; last = -1; continue;
        }
        if (ch === last) run++; else { run = 1; last = ch; }
        out[op++] = ch; if (op === out.length) flushOut();
      }
    }
  }
  flushOut();
}

// =====================================================================================
// Containers: zip, 7z, tar, single-stream compressors; worker thread + shared ring buffer
// =====================================================================================
const WT = require('worker_threads');
const zlib = require('zlib');
const cp = require('child_process');

const MIDI_EXT = /\.(mid|midi|kar|rmi|smf)$/i;
const RING_HDR = 4096 + 16;

// ---- shared ring (worker writes, main thread reads synchronously) ----------------------
function makeRing(size) {
  const sab = new SharedArrayBuffer(RING_HDR + size);
  return sab;
}
class RingWriter { // worker side
  constructor(sab) {
    this.h = new Int32Array(sab, 0, 4); this.err = new Uint8Array(sab, 16, 4096);
    this.d = new Uint8Array(sab, RING_HDR); this.cap = this.d.length; this.head = 0;
    this.stopped = false;
  }
  write(buf, a, b) {
    const h = this.h, d = this.d, cap = this.cap;
    while (a < b) {
      if (Atomics.load(h, 3) === 1) { this.stopped = true; throw STOP; } // reader closed
      let used = Atomics.load(h, 0);
      while (used === cap) { Atomics.wait(h, 0, cap, 200); if (Atomics.load(h, 3) === 1) { this.stopped = true; throw STOP; } used = Atomics.load(h, 0); }
      let n = Math.min(cap - used, b - a, cap - this.head);
      d.set(buf.subarray(a, a + n), this.head);
      this.head = (this.head + n) % cap; a += n;
      Atomics.add(h, 0, n); Atomics.notify(h, 0);
    }
  }
  finish() { Atomics.store(this.h, 1, 1); Atomics.notify(this.h, 0); }
  fail(msg) {
    const b = Buffer.from(String(msg)).subarray(0, 4000);
    this.err.set(b); Atomics.store(this.h, 2, b.length); Atomics.store(this.h, 1, 2); Atomics.notify(this.h, 0);
  }
}
const STOP = { stop: true };

class RingSource { // main side
  constructor(sab, worker) {
    this.h = new Int32Array(sab, 0, 4); this.err = new Uint8Array(sab, 16, 4096);
    this.d = new Uint8Array(sab, RING_HDR); this.cap = this.d.length; this.tail = 0; this.worker = worker;
  }
  read(dst, off, len) {
    const h = this.h;
    for (;;) {
      const avail = Atomics.load(h, 0);
      if (avail > 0) {
        const n = Math.min(avail, len, this.cap - this.tail);
        dst.set(this.d.subarray(this.tail, this.tail + n), off);
        this.tail = (this.tail + n) % this.cap;
        Atomics.sub(h, 0, n); Atomics.notify(h, 0);
        return n;
      }
      const st = Atomics.load(h, 1);
      if (st === 1) { if (Atomics.load(h, 0) > 0) continue; return 0; }
      if (st === 2) {
        const n = Atomics.load(h, 2);
        throw new Error(Buffer.from(this.err.subarray(0, n)).toString());
      }
      Atomics.wait(h, 0, 0, 500);
    }
  }
  close() { Atomics.store(this.h, 3, 1); Atomics.notify(this.h, 0); if (this.worker) { this.worker.terminate(); this.worker = null; } }
}
class FileSource { // plain file region as a sequential source
  constructor(fd, start, length) { this.fd = fd; this.pos = start; this.left = length; }
  read(dst, off, len) {
    const n = Math.min(len, this.left);
    if (n <= 0) return 0;
    const r = fs.readSync(this.fd, dst, off, n, this.pos);
    this.pos += r; this.left -= r;
    return r;
  }
  close() {}
}

function startJob(job, ringMB) {
  const sab = makeRing((ringMB || 32) << 20);
  const worker = new WT.Worker(__filename, { workerData: { midiprecJob: job, sab } });
  worker.unref();
  worker.on('error', () => {}); // errors are reported through the ring
  return new RingSource(sab, worker);
}

// ---- buffered sequential reader (main thread) ----------------------------------------
class SeqReader {
  constructor(src, bufSize) {
    this.src = src; this.buf = Buffer.allocUnsafe(bufSize || (4 << 20));
    this.pos = 0; this.end = 0; this.base = 0; this.eof = false;
  }
  tell() { return this.base + this.pos; }
  fill(n) {
    if (this.end - this.pos >= n) return true;
    if (this.eof) return false;
    if (this.pos > 0) { this.buf.copyWithin(0, this.pos, this.end); this.base += this.pos; this.end -= this.pos; this.pos = 0; }
    if (n > this.buf.length) { const nb = Buffer.allocUnsafe(n * 2); this.buf.copy(nb, 0, 0, this.end); this.buf = nb; }
    while (this.end - this.pos < n) {
      const r = this.src.read(this.buf, this.end, this.buf.length - this.end);
      if (r <= 0) { this.eof = true; return false; }
      this.end += r;
    }
    return true;
  }
  skip(n) {
    while (n > 0) {
      if (this.pos >= this.end) { if (!this.fill(1)) return false; }
      const k = Math.min(n, this.end - this.pos); this.pos += k; n -= k;
    }
    return true;
  }
}

// ---- compact note store (for MIDIs that arrive as a stream) -----------------------------
// per track: 0x00-0x7F note (same tick), 0x81-0xF9 note with delta 1..0x79 + key,
// 0xFA vel, 0xFB LEB delta + key, 0xFC LEB delta + tempo(3), 0xFD LEB delta end-of-track
const ARENA_BLOCK = 8 << 20;
class Arena {
  constructor(ramBudget, tmpDir) {
    this.blocks = []; this.starts = []; this.used = []; this.size = 0; this.ram = 0; this.budget = ramBudget; this.tmpDir = tmpDir;
    this.fd = -1; this.tmpPath = null; this.spilled = 0; this.lastBi = 0;
    this.cur = null; this.cp = 0; this.ce = 0; this.newBlock();
  }
  // seal the current block at its used length; data stays contiguous in the global address space
  newBlock() {
    let start = 0;
    if (this.cur) { const i = this.blocks.length - 1; this.used[i] = this.cp; start = this.starts[i] + this.cp; }
    this.cur = new Uint8Array(ARENA_BLOCK); this.cp = 0; this.ce = ARENA_BLOCK - 64;
    this.blocks.push(this.cur); this.starts.push(start); this.used.push(ARENA_BLOCK); this.ram += ARENA_BLOCK;
    if (this.ram > this.budget) this.spill();
  }
  spill() { // move the oldest in-memory sealed blocks to a temp file
    if (this.fd < 0) {
      this.tmpPath = path.join(this.tmpDir || os.tmpdir(), `midiprec-${process.pid}-${Date.now()}.tmp`);
      this.fd = fs.openSync(this.tmpPath, 'w+');
    }
    for (let i = 0; i < this.blocks.length - 1 && this.ram > this.budget; i++) {
      const b = this.blocks[i];
      if (!b) continue;
      fs.writeSync(this.fd, b, 0, this.used[i], i * ARENA_BLOCK);
      this.blocks[i] = null; this.ram -= ARENA_BLOCK; this.spilled += this.used[i];
    }
  }
  tell() { return this.starts[this.blocks.length - 1] + this.cp; }
  blockOf(g) {
    let i = this.lastBi;
    if (!(g >= this.starts[i] && g < this.starts[i] + this.used[i])) {
      let lo = 0, hi = this.blocks.length - 1;
      while (lo < hi) { const m = (lo + hi + 1) >> 1; if (this.starts[m] <= g) lo = m; else hi = m - 1; }
      i = lo;
    }
    this.lastBi = i;
    return i;
  }
  // copy n bytes at global position g into dst[off..]
  read(dst, off, g, n) {
    while (n > 0) {
      const bi = this.blockOf(g), bo = g - this.starts[bi];
      const k = Math.min(n, this.used[bi] - bo);
      if (k <= 0) throw new Error('internal: arena read out of range');
      const b = this.blocks[bi];
      if (b) dst.set(b.subarray(bo, bo + k), off);
      else fs.readSync(this.fd, dst, off, k, bi * ARENA_BLOCK + bo);
      off += k; g += k; n -= k;
    }
  }
  close() {
    this.blocks = []; this.cur = null;
    if (this.fd >= 0) { try { fs.closeSync(this.fd); fs.unlinkSync(this.tmpPath); } catch (e) { /* ignore */ } this.fd = -1; }
  }
}

// Parse one MIDI from a sequential reader into the arena. `limit` = byte length (Infinity if unknown).
function parseMidiStream(seq, limit, opts) {
  const startPos = seq.tell();
  const lim = limit === undefined ? Infinity : limit;
  const endPos = startPos + lim;
  // find MThd within the first 64 KB (RIFF/RMID wrappers, junk)
  let found = false;
  for (let scanned = 0; scanned < 65536; scanned++) {
    if (!seq.fill(4)) break;
    const b = seq.buf, p = seq.pos;
    if (b[p] === 0x4D && b[p + 1] === 0x54 && b[p + 2] === 0x68 && b[p + 3] === 0x64) { found = true; break; }
    seq.pos++;
  }
  if (!found || !seq.fill(14)) throw new Error('not a MIDI file (no MThd header)');
  const hb = seq.buf, hp = seq.pos;
  const hlen = hb.readUInt32BE(hp + 4);
  const format = hb.readUInt16BE(hp + 8), declared = hb.readUInt16BE(hp + 10), division = hb.readUInt16BE(hp + 12);
  seq.skip(8 + hlen);
  const arena = new Arena((opts.ramMB || 1024) * 1048576, opts.tmpDir);
  const tracks = [];
  const mask = opts.chanMask === undefined ? 0xFFFF : opts.chanMask;
  try {
    while (seq.tell() + 8 <= endPos && seq.fill(8)) {
      const b = seq.buf, p = seq.pos;
      const id = b.toString('latin1', p, p + 4);
      let len = b.readUInt32BE(p + 4);
      seq.pos += 8;
      if (seq.tell() + len > endPos) len = endPos - seq.tell();
      if (id === 'MTrk') tracks.push(encodeTrack(seq, len, arena, mask));
      else if (/^[\x20-\x7E]{4}$/.test(id)) { if (!seq.skip(len)) break; }
      else break;
    }
  } catch (e) { arena.close(); throw e; }
  if (lim !== Infinity) { const rest = endPos - seq.tell(); if (rest > 0) seq.skip(rest); }
  arena.size = arena.tell();
  return { format, declaredTracks: declared, division, tracks, arena };
}

function encodeTrack(seq, len, A, mask) {
  const start = A.tell();
  const chunkEnd = seq.tell() + len;
  let ob = A.cur, op = A.cp, oe = A.ce;
  let tick = 0, last = 0, rs = 0, vel = -1;
  const lebOut = (v) => { while (v >= 128) { ob[op++] = (v % 128) | 0x80; v = Math.floor(v / 128); } ob[op++] = v; };
  for (;;) {
    if (op > oe) { A.cp = op; A.newBlock(); ob = A.cur; op = A.cp; oe = A.ce; }
    const remain = chunkEnd - seq.tell();
    if (remain <= 0) break;
    if (seq.end - seq.pos < 16) seq.fill(Math.min(16, remain));
    const buf = seq.buf;
    let pos = seq.pos;
    const end = Math.min(seq.end, pos + remain);
    if (pos >= end) break;
    let b = buf[pos++], d = b & 0x7F;
    while (b & 0x80) { if (pos >= end) { pos = end; break; } b = buf[pos++]; d = d * 128 + (b & 0x7F); }
    if (pos >= end) { seq.pos = pos; break; }
    let st = buf[pos];
    if (st & 0x80) pos++; else st = rs;
    const hi = st >> 4;
    if (hi >= 8 && hi <= 0xE) {
      const n = (hi === 0xC || hi === 0xD) ? 1 : 2;
      if (pos + n > end) { seq.pos = end; break; }
      rs = st; tick += d;
      if (hi === 9) {
        const v = buf[pos + 1];
        if (v !== 0 && ((mask >> (st & 15)) & 1)) {
          const k = buf[pos] & 127;
          if (v !== vel) { ob[op++] = 0xFA; ob[op++] = v; vel = v; }
          const dd = tick - last; last = tick;
          if (dd === 0) ob[op++] = k;
          else if (dd < 0x7A) { ob[op++] = 0x80 + dd; ob[op++] = k; }
          else { ob[op++] = 0xFB; lebOut(dd); ob[op++] = k; }
        }
      }
      seq.pos = pos + n; continue;
    }
    tick += d;
    if (st === 0xFF || st === 0xF0 || st === 0xF7) {
      let type = 0;
      if (st === 0xFF) { if (pos >= end) { seq.pos = end; break; } type = buf[pos++]; }
      let l = 0, c = 0x80;
      while (c & 0x80) { if (pos >= end) break; c = buf[pos++]; l = l * 128 + (c & 0x7F); }
      if (c & 0x80) { seq.pos = end; break; }
      seq.pos = pos;
      if (st === 0xFF && type === 0x51 && l >= 3) {
        if (!seq.fill(3)) break;
        const q = seq.pos, bb = seq.buf;
        const us = (bb[q] << 16) | (bb[q + 1] << 8) | bb[q + 2];
        ob[op++] = 0xFC; lebOut(tick - last); last = tick;
        ob[op++] = us >> 16; ob[op++] = (us >> 8) & 255; ob[op++] = us & 255;
      }
      if (st === 0xFF && type === 0x2F) { seq.skip(Math.max(0, chunkEnd - seq.tell())); break; }
      const rem = chunkEnd - seq.tell();
      seq.skip(Math.min(l, rem));
      continue;
    }
    if (st < 0x80) { seq.pos = pos + 1; continue; }
    seq.pos = pos + (st === 0xF2 ? 2 : (st === 0xF1 || st === 0xF3) ? 1 : 0);
  }
  if (op > oe) { A.cp = op; A.newBlock(); ob = A.cur; op = A.cp; }
  ob[op++] = 0xFD; lebOut(tick - last);
  A.cp = op;
  return { start, length: A.tell() - start };
}

class CompactTrackReader {
  constructor(arena, start, length, bufSize) {
    this.A = arena; this.g = start; this.left = length;
    this.buf = new Uint8Array(Math.max(4096, Math.min(bufSize, length + 32)));
    this.pos = 0; this.end = 0;
    this.tick = 0; this.vel = 127; this.endTick = 0;
    this.pend = 0; this.pv = 0; this.pvel = 0;
    this.next = -1;
    this.step(PRIME_CTX);
  }
  refill() {
    const keep = this.end - this.pos;
    if (keep > 0) this.buf.copyWithin(0, this.pos, this.end);
    this.pos = 0; this.end = keep;
    const n = Math.min(this.buf.length - keep, this.left);
    if (n > 0) { this.A.read(this.buf, keep, this.g, n); this.g += n; this.left -= n; this.end += n; }
  }
  step(ctx) {
    const T = this.next;
    const words = ctx.words, cnt = ctx.cnt, vm = ctx.vm, vs = ctx.vs;
    if (this.pend === 1) { const k = this.pv, v = this.pvel; words[k >> 5] |= 1 << (k & 31); cnt[k]++; vs[k] += v; if (v > vm[k]) vm[k] = v; ctx.notes++; }
    else if (this.pend === 2) ctx.tempo(T, this.pv);
    this.pend = 0;
    const buf = this.buf;
    let pos = this.pos, end = this.end, tick = this.tick, vel = this.vel, notes = 0;
    for (;;) {
      if (end - pos < 16 && this.left > 0) { this.pos = pos; this.refill(); pos = this.pos; end = this.end; }
      if (pos >= end) { tick = this.tick; break; }
      const b = buf[pos++];
      if (b < 0x80) { // note, same tick as the previous event
        if (tick === T) { words[b >> 5] |= 1 << (b & 31); cnt[b]++; vs[b] += vel; if (vel > vm[b]) vm[b] = vel; notes++; continue; }
        this.pend = 1; this.pv = b; this.pvel = vel;
        this.pos = pos; this.tick = tick; this.vel = vel; ctx.notes += notes; this.next = tick; return tick;
      }
      if (b <= 0xF9) {
        tick += b - 0x80; const k = buf[pos++];
        this.pend = 1; this.pv = k; this.pvel = vel;
        this.pos = pos; this.tick = tick; this.vel = vel; ctx.notes += notes; this.next = tick; return tick;
      }
      if (b === 0xFA) { vel = buf[pos++]; continue; }
      let dd = 0, mul = 1, c;
      do { c = buf[pos++]; dd += (c & 0x7F) * mul; mul *= 128; } while (c & 0x80);
      tick += dd;
      if (b === 0xFB) {
        const k = buf[pos++];
        this.pend = 1; this.pv = k; this.pvel = vel;
        this.pos = pos; this.tick = tick; this.vel = vel; ctx.notes += notes; this.next = tick; return tick;
      }
      if (b === 0xFC) {
        const us = (buf[pos] << 16) | (buf[pos + 1] << 8) | buf[pos + 2]; pos += 3;
        if (dd === 0 && tick === T) { ctx.tempo(T, us); continue; }
        this.pend = 2; this.pv = us;
        this.pos = pos; this.tick = tick; this.vel = vel; ctx.notes += notes; this.next = tick; return tick;
      }
      // 0xFD end of track
      this.tick = tick; break;
    }
    this.pos = end; this.left = 0; this.endTick = tick; this.vel = vel;
    ctx.notes += notes; this.next = Infinity;
    return Infinity;
  }
}

// ---- 7z ------------------------------------------------------------------------------------
class BP { // byte parser for 7z headers
  constructor(b) { this.b = b; this.p = 0; }
  byte() { if (this.p >= this.b.length) throw new Error('corrupt 7z header'); return this.b[this.p++]; }
  bytes(n) { if (this.p + n > this.b.length) throw new Error('corrupt 7z header'); const r = this.b.subarray(this.p, this.p + n); this.p += n; return r; }
  num() {
    const first = this.byte(); let mask = 0x80, v = 0;
    for (let i = 0; i < 8; i++) {
      if ((first & mask) === 0) return v + (first & (mask - 1)) * Math.pow(2, 8 * i);
      v += this.byte() * Math.pow(2, 8 * i); mask >>= 1;
    }
    return v;
  }
  bits(n) { const r = []; let b = 0, m = 0; for (let i = 0; i < n; i++) { if (m === 0) { b = this.byte(); m = 0x80; } r.push((b & m) !== 0); m >>= 1; } return r; }
  defined(n) { const all = this.byte(); return all ? new Array(n).fill(true) : this.bits(n); }
}
function sz7Folder(p) {
  const nc = p.num();
  const coders = []; let tin = 0, tout = 0;
  for (let i = 0; i < nc; i++) {
    const fl = p.byte();
    const id = Buffer.from(p.bytes(fl & 0x0F)).toString('hex');
    let ni = 1, no = 1;
    if (fl & 0x10) { ni = p.num(); no = p.num(); }
    let props = null;
    if (fl & 0x20) props = Buffer.from(p.bytes(p.num()));
    if (fl & 0x80) throw new Error('unsupported 7z coder flags');
    coders.push({ id, ni, no, props, inBase: tin, outBase: tout });
    tin += ni; tout += no;
  }
  const bind = [];
  for (let i = 0; i < tout - 1; i++) bind.push({ in: p.num(), out: p.num() });
  const npk = tin - bind.length;
  const packed = [];
  if (npk === 1) { for (let i = 0; i < tin; i++) if (!bind.some((x) => x.in === i)) { packed.push(i); break; } }
  else for (let i = 0; i < npk; i++) packed.push(p.num());
  return { coders, bind, packed, tin, tout, unpack: [], crc: false };
}
function sz7StreamsInfo(p) {
  const si = { packPos: 0, packSizes: [], folders: [], nus: null, sub: null };
  for (;;) {
    const t = p.byte();
    if (t === 0) break;
    if (t === 0x06) {
      si.packPos = p.num(); const n = p.num();
      for (;;) { const u = p.byte(); if (u === 0) break; if (u === 0x09) for (let i = 0; i < n; i++) si.packSizes.push(p.num()); else if (u === 0x0A) { const d = p.defined(n); p.bytes(4 * d.filter(Boolean).length); } else throw new Error('corrupt 7z header'); }
    } else if (t === 0x07) {
      if (p.byte() !== 0x0B) throw new Error('corrupt 7z header');
      const nf = p.num();
      if (p.byte() !== 0) throw new Error('unsupported 7z header (external)');
      for (let i = 0; i < nf; i++) si.folders.push(sz7Folder(p));
      if (p.byte() !== 0x0C) throw new Error('corrupt 7z header');
      for (const f of si.folders) for (let j = 0; j < f.tout; j++) f.unpack.push(p.num());
      for (;;) { const u = p.byte(); if (u === 0) break; if (u === 0x0A) { const d = p.defined(nf); d.forEach((x, i) => { si.folders[i].crc = x; }); p.bytes(4 * d.filter(Boolean).length); } else throw new Error('corrupt 7z header'); }
    } else if (t === 0x08) {
      const nf = si.folders.length;
      const nus = new Array(nf).fill(1);
      let u = p.byte();
      if (u === 0x0D) { for (let i = 0; i < nf; i++) nus[i] = p.num(); u = p.byte(); }
      const sub = [];
      for (let i = 0; i < nf; i++) {
        const n = nus[i], arr = [];
        if (n > 0) {
          let sum = 0;
          if (u === 0x09) for (let j = 0; j < n - 1; j++) { const s = p.num(); arr.push(s); sum += s; }
          arr.push(folderUnpackSize(si.folders[i]) - sum);
        }
        sub.push(arr);
      }
      if (u === 0x09) u = p.byte();
      while (u !== 0) {
        if (u === 0x0A) {
          let nd = 0; for (let i = 0; i < nf; i++) if (!(nus[i] === 1 && si.folders[i].crc)) nd += nus[i];
          const d = p.defined(nd); p.bytes(4 * d.filter(Boolean).length);
        } else throw new Error('corrupt 7z header');
        u = p.byte();
      }
      si.nus = nus; si.sub = sub;
    } else throw new Error('corrupt 7z header');
  }
  if (!si.nus) { si.nus = si.folders.map(() => 1); si.sub = si.folders.map((f) => [folderUnpackSize(f)]); }
  return si;
}
function folderMainOut(f) { for (let i = 0; i < f.tout; i++) if (!f.bind.some((b) => b.out === i)) return i; return 0; }
function folderUnpackSize(f) { return f.unpack[folderMainOut(f)]; }
// Linear coder chain from the packed stream to the folder output: [main decoder, filter, filter, ...]
function folderChain(f) {
  if (f.packed.length !== 1) throw new Error('unsupported 7z coder graph (BCJ2 / multi-stream)');
  const chain = [];
  let outIdx = folderMainOut(f);
  for (let guard = 0; guard < 8; guard++) {
    const ci = f.coders.findIndex((c) => outIdx >= c.outBase && outIdx < c.outBase + c.no);
    const c = f.coders[ci];
    if (c.ni !== 1 || c.no !== 1) throw new Error('unsupported 7z coder graph');
    chain.unshift({ id: c.id, props: c.props ? Array.from(c.props) : null, size: f.unpack[c.outBase] });
    const inIdx = c.inBase;
    const bp = f.bind.find((b) => b.in === inIdx);
    if (!bp) { if (inIdx !== f.packed[0]) throw new Error('unsupported 7z coder graph'); return chain; }
    outIdx = bp.out;
  }
  throw new Error('unsupported 7z coder graph');
}
function read7z(file) {
  const fd = fs.openSync(file, 'r');
  try {
    const sig = Buffer.alloc(32); fs.readSync(fd, sig, 0, 32, 0);
    const nhOff = Number(sig.readBigUInt64LE(12)), nhSize = Number(sig.readBigUInt64LE(20));
    const fsize = fs.fstatSync(fd).size;
    if (!(nhSize > 0 && nhSize < (1 << 30) && 32 + nhOff + nhSize <= fsize)) throw new Error('corrupt or truncated 7z archive');
    let hdr = Buffer.alloc(nhSize); fs.readSync(fd, hdr, 0, nhSize, 32 + nhOff);
    let p = new BP(hdr), id = p.byte();
    while (id === 0x17) { // encoded header
      const si = sz7StreamsInfo(p);
      const f = si.folders[0];
      const chunks = [];
      decodeFolderSync(fd, 32 + si.packPos, si.packSizes[0], folderChain(f), (b, s, e) => chunks.push(Buffer.from(b.subarray(s, e))));
      hdr = Buffer.concat(chunks); p = new BP(hdr); id = p.byte();
    }
    if (id !== 0x01) throw new Error('corrupt 7z header');
    let si = null, files = [];
    for (;;) {
      const t = p.byte();
      if (t === 0) break;
      if (t === 0x02) { for (;;) { const u = p.byte(); if (u === 0) break; p.bytes(p.num()); } }
      else if (t === 0x03) sz7StreamsInfo(p);
      else if (t === 0x04) si = sz7StreamsInfo(p);
      else if (t === 0x05) files = sz7Files(p);
      else throw new Error('corrupt 7z header');
    }
    const entries = [];
    if (si) {
      let fi = 0, k = 0, off = 0;
      for (const f of files) {
        if (f.empty) continue;
        while (fi < si.folders.length && si.nus[fi] === 0) fi++;
        if (fi >= si.folders.length) break;
        const size = si.sub[fi][k];
        entries.push({ name: f.name, size, folder: fi, offset: off });
        off += size; k++;
        if (k >= si.nus[fi]) { fi++; k = 0; off = 0; }
      }
      // pack stream positions
      let pos = 32 + si.packPos, pi = 0;
      for (const f of si.folders) {
        f.packStart = pos; f.packSize = 0;
        for (let j = 0; j < f.packed.length; j++) { f.packSize += si.packSizes[pi]; pos += si.packSizes[pi]; pi++; }
      }
    }
    return { entries, folders: si ? si.folders : [] };
  } finally { fs.closeSync(fd); }
}
function sz7Files(p) {
  const n = p.num();
  const files = []; for (let i = 0; i < n; i++) files.push({ name: '', empty: false });
  let emptyStreams = [];
  for (;;) {
    const t = p.num();
    if (t === 0) break;
    const size = p.num(), end = p.p + size;
    if (t === 0x0E) { emptyStreams = p.bits(n); emptyStreams.forEach((e, i) => { files[i].empty = e; }); }
    else if (t === 0x11) {
      if (p.byte() !== 0) throw new Error('unsupported 7z header (external names)');
      for (let i = 0; i < n; i++) {
        let s = '';
        for (;;) { const lo = p.byte(), hi = p.byte(), c = lo | (hi << 8); if (c === 0) break; s += String.fromCharCode(c); }
        files[i].name = s;
      }
    }
    p.p = end;
  }
  return files;
}

// decode a 7z folder synchronously: packed bytes at [start, start+size) through `chain`, output to sink
function decodeFolderSync(fd, start, size, chain, sink) {
  const src = new InStream(fd, start, size);
  let out = sink;
  for (let i = chain.length - 1; i >= 1; i--) {
    const c = chain[i];
    if (c.id === '03') out = deltaFilter((c.props ? c.props[0] : 0) + 1, out);
    else if (c.id === '00') { /* copy */ }
    else throw new Error(`unsupported 7z filter ${c.id} (install 7-Zip to read this archive)`);
  }
  const m = chain[0];
  const total = chain[chain.length - 1].size;
  if (m.id === '21') decodeLzma2(src, lzma2DictSize(m.props[0]), out, m.size);
  else if (m.id === '030101') {
    const pr = m.props; const dict = pr[1] + pr[2] * 256 + pr[3] * 65536 + pr[4] * 16777216;
    decodeLzma1(src, pr[0], dict, m.size, out);
  } else if (m.id === '00') src.pipeTo(out, m.size);
  else if (m.id === '040202') decodeBzip2(src, out);
  else if (m.id === '040108') {
    const buf = Buffer.alloc(size); fs.readSync(fd, buf, 0, size, start);
    const r = zlib.inflateRawSync(buf); out(r, 0, r.length);
  } else throw new Error(`unsupported 7z method ${m.id} (install 7-Zip to read this archive)`);
  return total;
}

// ---- zip ----------------------------------------------------------------------------------
function readZip(file) {
  const fd = fs.openSync(file, 'r');
  try {
    const size = fs.fstatSync(fd).size;
    const tl = Math.min(size, 65557 + 22);
    const tail = Buffer.alloc(tl); fs.readSync(fd, tail, 0, tl, size - tl);
    let e = -1;
    for (let i = tl - 22; i >= 0; i--) if (tail.readUInt32LE(i) === 0x06054b50) { e = i; break; }
    if (e < 0) throw new Error('corrupt zip (no end of central directory)');
    let count = tail.readUInt16LE(e + 10), cdSize = tail.readUInt32LE(e + 12), cdOff = tail.readUInt32LE(e + 16);
    if (e >= 20 && tail.readUInt32LE(e - 20) === 0x07064b50) {
      const z64 = Number(tail.readBigUInt64LE(e - 20 + 8));
      const zb = Buffer.alloc(56); fs.readSync(fd, zb, 0, 56, z64);
      if (zb.readUInt32LE(0) === 0x06064b50) { count = Number(zb.readBigUInt64LE(32)); cdSize = Number(zb.readBigUInt64LE(40)); cdOff = Number(zb.readBigUInt64LE(48)); }
    }
    const cd = Buffer.alloc(cdSize); fs.readSync(fd, cd, 0, cdSize, cdOff);
    const entries = [];
    let p = 0;
    for (let i = 0; i < count && p + 46 <= cd.length; i++) {
      if (cd.readUInt32LE(p) !== 0x02014b50) break;
      const flags = cd.readUInt16LE(p + 8), method = cd.readUInt16LE(p + 10);
      let csize = cd.readUInt32LE(p + 20), usize = cd.readUInt32LE(p + 24);
      const nl = cd.readUInt16LE(p + 28), xl = cd.readUInt16LE(p + 30), cl = cd.readUInt16LE(p + 32);
      let loff = cd.readUInt32LE(p + 42);
      const name = cd.toString(flags & 0x800 ? 'utf8' : 'latin1', p + 46, p + 46 + nl);
      let x = p + 46 + nl; const xe = x + xl;
      while (x + 4 <= xe) {
        const hid = cd.readUInt16LE(x), hs = cd.readUInt16LE(x + 2);
        if (hid === 0x0001) {
          let q = x + 4;
          if (usize === 0xFFFFFFFF) { usize = Number(cd.readBigUInt64LE(q)); q += 8; }
          if (csize === 0xFFFFFFFF) { csize = Number(cd.readBigUInt64LE(q)); q += 8; }
          if (loff === 0xFFFFFFFF) { loff = Number(cd.readBigUInt64LE(q)); q += 8; }
        }
        x += 4 + hs;
      }
      entries.push({ name, method, flags, csize, size: usize, loff });
      p += 46 + nl + xl + cl;
    }
    for (const en of entries) {
      const lh = Buffer.alloc(30); fs.readSync(fd, lh, 0, 30, en.loff);
      if (lh.readUInt32LE(0) !== 0x04034b50) { en.dataStart = -1; continue; }
      en.dataStart = en.loff + 30 + lh.readUInt16LE(26) + lh.readUInt16LE(28);
    }
    return entries;
  } finally { fs.closeSync(fd); }
}

// ---- tar ----------------------------------------------------------------------------------
function tarNum(b, off, len) {
  if (b[off] & 0x80) { let v = 0; for (let i = 1; i < len; i++) v = v * 256 + b[off + i]; return v; }
  const s = b.toString('latin1', off, off + len).replace(/\0.*$/, '').trim();
  return s ? parseInt(s, 8) : 0;
}
function tarStr(b, off, len) { const s = b.toString('utf8', off, off + len); const z = s.indexOf('\0'); return z >= 0 ? s.slice(0, z) : s; }
function isTarHeader(b, off) { return b.length >= off + 512 && b.toString('latin1', off + 257, off + 262) === 'ustar'; }
// iterate tar entries from a SeqReader; cb(name, size) must consume exactly `size` bytes or return false
function tarWalk(seq, cb) {
  let longName = null, paxName = null, paxSize = null;
  for (;;) {
    if (!seq.fill(512)) return;
    const h = seq.buf.subarray(seq.pos, seq.pos + 512);
    if (h.every((x) => x === 0)) return;
    seq.pos += 512;
    let name = tarStr(h, 0, 100);
    const prefix = tarStr(h, 345, 155);
    if (prefix && h.toString('latin1', 257, 262) === 'ustar') name = prefix + '/' + name;
    let size = tarNum(h, 124, 12);
    const type = String.fromCharCode(h[156] || 48);
    const padded = Math.ceil(size / 512) * 512;
    if (type === 'L' || type === 'x' || type === 'g') {
      if (!seq.fill(Math.min(padded, 1 << 20))) return;
      const data = Buffer.from(seq.buf.subarray(seq.pos, seq.pos + size));
      seq.skip(padded);
      if (type === 'L') longName = data.toString('utf8').replace(/\0.*$/s, '');
      else if (type === 'x') {
        const txt = data.toString('utf8');
        for (const m of txt.matchAll(/\d+ (path|size)=([^\n]*)\n/g)) { if (m[1] === 'path') paxName = m[2]; else paxSize = +m[2]; }
      }
      continue;
    }
    if (longName) { name = longName; longName = null; }
    if (paxName) { name = paxName; paxName = null; }
    if (paxSize !== null) { size = paxSize; paxSize = null; }
    const pad = Math.ceil(size / 512) * 512 - size;
    if (type === '0' || type === '\0' || type === '7') {
      const before = seq.tell();
      const keep = cb(name, size);
      const used = seq.tell() - before;
      if (used < size) seq.skip(size - used);
      seq.skip(pad);
      if (keep === false) return;
    } else seq.skip(size + pad);
  }
}

// ---- external 7-Zip (fallback / speed) -------------------------------------------------
let SEVENZIP = undefined;
function find7z() {
  if (SEVENZIP !== undefined) return SEVENZIP;
  const cands = process.platform === 'win32'
    ? ['7z.exe', '7za.exe', '7zz.exe', 'C:\\Program Files\\7-Zip\\7z.exe', 'C:\\Program Files (x86)\\7-Zip\\7z.exe']
    : ['7zz', '7z', '7za'];
  for (const c of cands) {
    try { const r = cp.spawnSync(c, ['i'], { stdio: ['ignore', 'pipe', 'ignore'], timeout: 5000, windowsHide: true }); if (r.status === 0 || (r.stdout && r.stdout.length)) { SEVENZIP = c; return c; } } catch (e) { /* not found */ }
  }
  SEVENZIP = null; return null;
}
function hasTool(name) {
  try { const r = cp.spawnSync(name, ['--version'], { stdio: ['ignore', 'ignore', 'ignore'], timeout: 5000, windowsHide: true }); return r.status === 0; } catch (e) { return false; }
}
function list7zExternal(exe, file) {
  const r = cp.spawnSync(exe, ['l', '-slt', '-ba', '--', file], { maxBuffer: 1 << 28, windowsHide: true });
  if (r.status !== 0) throw new Error('7-Zip could not list ' + file + ': ' + String(r.stderr || '').trim());
  const out = r.stdout.toString('utf8');
  const entries = [];
  for (const block of out.split(/\r?\n\r?\n/)) {
    const m = /^Path = (.*)$/m.exec(block);
    if (!m) continue;
    const isDir = /^Folder = \+/m.test(block) || /^Attributes = D/m.test(block);
    const sz = /^Size = (\d+)/m.exec(block);
    if (!isDir) entries.push({ name: m[1].trim(), size: sz ? +sz[1] : -1 });
  }
  return entries;
}

// ---- worker side ------------------------------------------------------------------------
async function workerMain() {
  const { midiprecJob: job, sab } = WT.workerData;
  const w = new RingWriter(sab);
  const sink = (b, s, e) => w.write(b, s, e);
  try {
    if (job.kind === 'zlib') {
      await new Promise((resolve, reject) => {
        const rs = fs.createReadStream(job.file, { start: job.start, end: job.start + job.length - 1, highWaterMark: 1 << 20 });
        const z = job.method === 'gunzip' ? zlib.createGunzip()
          : job.method === 'brotli' ? zlib.createBrotliDecompress()
            : job.method === 'zstd' ? zlib.createZstdDecompress()
              : zlib.createInflateRaw();
        z.on('data', (c) => { try { w.write(c, 0, c.length); } catch (e) { rs.destroy(); z.destroy(); if (e === STOP) resolve(); else reject(e); } });
        z.on('end', resolve); z.on('error', reject); rs.on('error', reject);
        rs.pipe(z);
      });
    } else if (job.kind === 'external') {
      await new Promise((resolve, reject) => {
        const ch = cp.spawn(job.cmd, job.args, { stdio: ['ignore', 'pipe', 'pipe'], windowsHide: true });
        let err = '';
        ch.stderr.on('data', (d) => { if (err.length < 4000) err += d; });
        ch.stdout.on('data', (c) => { try { w.write(c, 0, c.length); } catch (e) { ch.kill(); if (e === STOP) resolve(); else reject(e); } });
        ch.on('error', reject);
        ch.on('close', (code) => { if (code === 0 || w.stopped) resolve(); else reject(new Error(`${path.basename(job.cmd)} failed (${code}): ${err.trim().split('\n').pop()}`)); });
      });
    } else {
      const fd = fs.openSync(job.file, 'r');
      try {
        const src = new InStream(fd, job.start, job.length);
        if (job.kind === 'xz') decodeXz(src, sink);
        else if (job.kind === 'bzip2') decodeBzip2(src, sink);
        else if (job.kind === 'lzma') {
          const h = src.bytes(13);
          const dict = h[1] + h[2] * 256 + h[3] * 65536 + h[4] * 16777216;
          let size = 0, unk = true; for (let i = 12; i >= 5; i--) { size = size * 256 + h[i]; if (h[i] !== 255) unk = false; }
          decodeLzma1(src, h[0], dict, unk ? -1 : size, sink);
        } else if (job.kind === 'zip-lzma') {
          const h = src.bytes(4); const ps = h[2] | (h[3] << 8); const pr = src.bytes(ps);
          const dict = pr[1] + pr[2] * 256 + pr[3] * 65536 + pr[4] * 16777216;
          decodeLzma1(src, pr[0], dict, job.eos ? -1 : job.size, sink);
        } else if (job.kind === '7z') {
          // pass through only the requested [offset, offset+size) ranges of the folder output
          const ranges = job.ranges; let ri = 0, pos = 0;
          const filt = (b, s, e) => {
            while (s < e && ri < ranges.length) {
              const [ro, rl] = ranges[ri];
              if (pos + (e - s) <= ro) { pos += e - s; return; }
              if (pos < ro) { s += ro - pos; pos = ro; }
              const take = Math.min(e - s, ro + rl - pos);
              if (take > 0) { w.write(b, s, s + take); s += take; pos += take; }
              if (pos >= ro + rl) ri++;
            }
            if (ri >= ranges.length) throw STOP;
          };
          decodeFolderSync(fd, job.start, job.length, job.chain, filt);
        } else throw new Error('unknown job ' + job.kind);
      } finally { fs.closeSync(fd); }
    }
    w.finish();
  } catch (e) {
    if (e === STOP) w.finish();
    else w.fail(e && e.message ? e.message : String(e));
  }
}

// =====================================================================================
// Input planning: plain MIDI / zip / 7z / tar / xz / gz / bz2 / zst / lzma / br / rar
// =====================================================================================
function sniff(file) {
  const fd = fs.openSync(file, 'r');
  const b = Buffer.alloc(65536);
  const n = fs.readSync(fd, b, 0, b.length, 0);
  fs.closeSync(fd);
  const h = b.subarray(0, n);
  const s4 = n >= 4 ? h.toString('latin1', 0, 4) : '';
  if (s4 === 'MThd') return 'midi';
  if (s4 === 'RIFF' && n >= 12 && h.toString('latin1', 8, 12) === 'RMID') return 'midi';
  if (n >= 4 && h[0] === 0x50 && h[1] === 0x4B && (h[2] === 3 || h[2] === 5 || h[2] === 7)) return 'zip';
  if (n >= 6 && h[0] === 0x37 && h[1] === 0x7A && h[2] === 0xBC && h[3] === 0xAF && h[4] === 0x27 && h[5] === 0x1C) return '7z';
  if (n >= 6 && h[0] === 0xFD && h[1] === 0x37 && h[2] === 0x7A && h[3] === 0x58 && h[4] === 0x5A && h[5] === 0) return 'xz';
  if (n >= 2 && h[0] === 0x1F && h[1] === 0x8B) return 'gzip';
  if (n >= 4 && h[0] === 0x42 && h[1] === 0x5A && h[2] === 0x68 && h[3] >= 0x31 && h[3] <= 0x39) return 'bzip2';
  if (n >= 4 && h[0] === 0x28 && h[1] === 0xB5 && h[2] === 0x2F && h[3] === 0xFD) return 'zstd';
  if (s4 === 'Rar!') return 'rar';
  if (isTarHeader(h, 0)) return 'tar';
  const ext = path.extname(file).toLowerCase();
  if (ext === '.lzma') return 'lzma';
  if (ext === '.br') return 'brotli';
  if (findMThd(h) >= 0) return 'midi';
  return 'unknown';
}

function wantEntry(name, opts) {
  if (opts.entry) return name.toLowerCase().includes(opts.entry.toLowerCase());
  return MIDI_EXT.test(name);
}

// time + result wrapper around analyzeReaders
function finishAnalysis(display, info, readers, division, P, opts, t0, extra) {
  const res = analyzeReaders(readers, division, P, opts);
  const fin = res.acc.finish();
  const t1 = process.hrtime.bigint();
  let csvInfo = null, csvPath = null;
  if (opts.csv) { csvPath = csvPathFor(opts, display); csvInfo = writeCsv(csvPath, res.acc, res.endTime, fin.L, opts); }
  const t2 = process.hrtime.bigint();
  if (process.env.MIDIPREC_DEBUG) {
    const a = res.acc;
    console.error(`[debug] segs=${a.nSeg} ents=${a.nEnt} compactions=${a.compactions || 0} fineBins=${a.fLen} finalEvals=${fin.nEval} onlineSolves=${a.solveCount} gridInits=${a.gridInits}`);
    console.error(`[debug] hands: gestures=${res.hands.stats.gestures} twoHand=${res.hands.stats.twoHand} capped=${res.hands.stats.capped} cost[reach,move,speed,repeat,spread,cross]=${Array.from(res.hands.cs).map((x) => (x / Math.max(1, res.hands.stats.gestures)).toFixed(3))}`);
  }
  return Object.assign({ file: display, mf: info, res, fin, csvInfo, csvPath, secs: Number(t1 - t0) / 1e9, csvSecs: Number(t2 - t1) / 1e9 }, extra || {});
}
let csvUsed = new Map();
function csvPathFor(opts, display) {
  let p = opts.csv;
  const last = display.split(' › ').pop();
  const base = path.basename(last.replace(/\\/g, '/')).replace(/\.[^.]*$/, '').replace(/[<>:"/\\|?*\x00-\x1f]/g, '_') || 'midi';
  if (p.includes('{name}')) return p.replace(/\{name\}/g, base);
  const n = (csvUsed.get(p) || 0) + 1; csvUsed.set(p, n);
  if (n === 1) return p;
  const ext = path.extname(p);
  return p.slice(0, p.length - ext.length) + '_' + n + ext;
}

function analyzePlain(file, base, length, display, opts, P) {
  const t0 = process.hrtime.bigint();
  const mf = openMidiFile(file, base, length);
  try {
    const total = mf.tracks.reduce((s, t) => s + t.length, 0);
    const budget = 64 << 20;
    const readers = mf.tracks.map((t) => new FileTrackReader(mf.fd, t.start, t.length,
      Math.max(EV_BUF_MIN, Math.min(1 << 20, Math.floor(budget * t.length / Math.max(1, total)) || EV_BUF_MIN)), opts.chanMask));
    return finishAnalysis(display, { format: mf.format, tracks: mf.tracks.length }, readers, mf.division, P, opts, t0);
  } finally { fs.closeSync(mf.fd); }
}

function analyzeSeq(seq, limit, display, opts, P, t0) {
  t0 = t0 || process.hrtime.bigint();
  const parsed = parseMidiStream(seq, limit, opts);
  try {
    const total = parsed.tracks.reduce((s, t) => s + t.length, 0);
    const budget = 64 << 20;
    const readers = parsed.tracks.map((t) => new CompactTrackReader(parsed.arena, t.start, t.length,
      Math.max(4096, Math.min(1 << 20, Math.floor(budget * t.length / Math.max(1, total)) || 4096))));
    return finishAnalysis(display, { format: parsed.format, tracks: parsed.tracks.length }, readers, parsed.division, P, opts, t0,
      { compactMB: parsed.arena.size / 1048576, spilledMB: parsed.arena.spilled / 1048576 });
  } finally { parsed.arena.close(); }
}

// decoded stream that is either one MIDI or a tar of MIDIs
function analyzeDecoded(src, display, opts, P, out) {
  const t0 = process.hrtime.bigint();
  const seq = new SeqReader(src);
  try {
    seq.fill(512);
    const head = seq.buf.subarray(seq.pos, seq.end);
    if (isTarHeader(head, 0)) {
      tarWalk(seq, (name, size) => {
        if (wantEntry(name, opts)) out.push(analyzeSeq(seq, size, display + ' › ' + name, opts, P));
        return true;
      });
    } else out.push(analyzeSeq(seq, Infinity, display, opts, P, t0));
  } finally { src.close(); }
}

function externalTool(kind) {
  if (kind === 'xz' && hasTool('xz')) return { cmd: 'xz', args: ['-dc', '--'] };
  if (kind === 'bzip2' && hasTool('bzip2')) return { cmd: 'bzip2', args: ['-dc', '--'] };
  if (kind === 'zstd' && hasTool('zstd')) return { cmd: 'zstd', args: ['-dc', '--'] };
  const z = find7z();
  if (z) return { cmd: z, args: ['x', '-so', '-bd', '--'] };
  return null;
}

function analyzeInput(file, opts, P) {
  const out = [];
  const kind = sniff(file);
  const size = fs.statSync(file).size;
  const dec = opts.decoder || 'auto';
  const name = path.basename(file);
  const ringMB = opts.ringMB || 32;
  if (kind === 'midi') { out.push(analyzePlain(file, 0, size, file, opts, P)); return out; }
  if (kind === 'unknown') throw new Error('unknown file type (not MIDI / zip / 7z / tar / xz / gz / bz2 / zst / lzma / rar)');

  if (['xz', 'gzip', 'bzip2', 'zstd', 'lzma', 'brotli'].includes(kind)) {
    let job;
    const ext = (dec !== 'builtin' && ['xz', 'bzip2', 'lzma'].includes(kind)) || (kind === 'zstd' && typeof zlib.createZstdDecompress !== 'function')
      ? externalTool(kind === 'lzma' ? 'lzma' : kind) : null;
    if (ext && dec !== 'builtin') job = { kind: 'external', cmd: ext.cmd, args: ext.args.concat([file]) };
    else if (kind === 'gzip') job = { kind: 'zlib', method: 'gunzip', file, start: 0, length: size };
    else if (kind === 'zstd') job = { kind: 'zlib', method: 'zstd', file, start: 0, length: size };
    else if (kind === 'brotli') job = { kind: 'zlib', method: 'brotli', file, start: 0, length: size };
    else job = { kind, file, start: 0, length: size };
    analyzeDecoded(startJob(job, ringMB), file, opts, P, out);
    return out;
  }

  if (kind === 'tar') {
    const src = new FileSource(fs.openSync(file, 'r'), 0, size);
    const seq = new SeqReader(src, 1 << 16);
    const items = [];
    tarWalk(seq, (n, sz) => { if (wantEntry(n, opts)) items.push({ name: n, start: seq.tell(), size: sz }); return true; });
    fs.closeSync(src.fd);
    for (const it of items) out.push(analyzePlain(file, it.start, it.size, file + ' › ' + it.name, opts, P));
    if (!items.length) throw new Error('no MIDI files in this tar');
    return out;
  }

  if (kind === 'zip') {
    const entries = readZip(file).filter((e) => !e.name.endsWith('/'));
    let pick = entries.filter((e) => wantEntry(e.name, opts));
    if (!pick.length && entries.length === 1 && !opts.entry) pick = entries;
    if (!pick.length) throw new Error('no MIDI files in this zip' + (opts.entry ? ` matching "${opts.entry}"` : ''));
    for (const e of pick) {
      const disp = file + ' › ' + e.name;
      if (e.flags & 1) throw new Error(`${e.name}: encrypted zip entries are not supported`);
      if (e.dataStart < 0) throw new Error(`${e.name}: corrupt local header`);
      if (e.method === 0) { out.push(analyzePlain(file, e.dataStart, e.csize, disp, opts, P)); continue; }
      let job = null;
      if (e.method === 8) job = { kind: 'zlib', method: 'inflateRaw', file, start: e.dataStart, length: e.csize };
      else if (e.method === 12) job = { kind: 'bzip2', file, start: e.dataStart, length: e.csize };
      else if (e.method === 14) job = { kind: 'zip-lzma', file, start: e.dataStart, length: e.csize, size: e.size, eos: !!(e.flags & 2) };
      else if (e.method === 93 && typeof zlib.createZstdDecompress === 'function') job = { kind: 'zlib', method: 'zstd', file, start: e.dataStart, length: e.csize };
      else if (e.method === 95) job = { kind: 'xz', file, start: e.dataStart, length: e.csize };
      else {
        const z = find7z();
        if (!z) throw new Error(`${e.name}: zip method ${e.method} needs 7-Zip installed`);
        job = { kind: 'external', cmd: z, args: ['x', '-so', '-bd', '--', file, e.name] };
      }
      const t0 = process.hrtime.bigint();
      const src = startJob(job, ringMB);
      try { out.push(analyzeSeq(new SeqReader(src), e.size, disp, opts, P, t0)); } finally { src.close(); }
    }
    return out;
  }

  if (kind === '7z') {
    let info = null, err = null;
    try { info = read7z(file); } catch (e) { err = e; }
    const z = dec !== 'builtin' ? find7z() : null;
    if (info) {
      let pick = info.entries.filter((e) => wantEntry(e.name, opts));
      if (!pick.length && info.entries.length === 1 && !opts.entry) pick = info.entries;
      if (!pick.length) throw new Error('no MIDI files in this 7z' + (opts.entry ? ` matching "${opts.entry}"` : ''));
      // group by folder (solid blocks are decoded once)
      const byFolder = new Map();
      for (const e of pick) { if (!byFolder.has(e.folder)) byFolder.set(e.folder, []); byFolder.get(e.folder).push(e); }
      for (const [fi, list] of byFolder) {
        list.sort((a, b) => a.offset - b.offset);
        const f = info.folders[fi];
        let chain = null;
        try { chain = folderChain(f); for (const c of chain) if (!['21', '030101', '00', '03', '040202', '040108'].includes(c.id)) throw new Error('unsupported 7z method ' + c.id); } catch (e) { chain = null; err = e; }
        if (z && (dec === 'external' || dec === 'auto' || !chain)) {
          for (const e of list) {
            const t0 = process.hrtime.bigint();
            const src = startJob({ kind: 'external', cmd: z, args: ['x', '-so', '-bd', '--', file, e.name] }, ringMB);
            try { out.push(analyzeSeq(new SeqReader(src), e.size, file + ' › ' + e.name, opts, P, t0)); } finally { src.close(); }
          }
          continue;
        }
        if (!chain) throw new Error(`${err.message} — install 7-Zip to read this archive`);
        const t0 = process.hrtime.bigint();
        const src = startJob({ kind: '7z', file, start: f.packStart, length: f.packSize, chain, ranges: list.map((e) => [e.offset, e.size]) }, ringMB);
        try {
          const seq = new SeqReader(src);
          for (const e of list) out.push(analyzeSeq(seq, e.size, file + ' › ' + e.name, opts, P, out.length ? undefined : t0));
        } finally { src.close(); }
      }
      return out;
    }
    if (!z) throw new Error(`${err ? err.message : '7z'} — install 7-Zip to read this archive`);
    // fall through to external listing
  }

  // rar / unreadable 7z: external 7-Zip
  const z = find7z();
  if (!z) throw new Error(`${kind} archives need 7-Zip installed (7z / 7zz / 7za in PATH)`);
  const entries = list7zExternal(z, file);
  let pick = entries.filter((e) => wantEntry(e.name, opts));
  if (!pick.length && entries.length === 1 && !opts.entry) pick = entries;
  if (!pick.length) throw new Error('no MIDI files in this archive');
  for (const e of pick) {
    const t0 = process.hrtime.bigint();
    const src = startJob({ kind: 'external', cmd: z, args: ['x', '-so', '-bd', '--', file, e.name] }, ringMB);
    try { out.push(analyzeSeq(new SeqReader(src), e.size >= 0 ? e.size : Infinity, file + ' › ' + e.name, opts, P, t0)); } finally { src.close(); }
  }
  return out;
}

module.exports = { analyzeInput, analyzeReaders, openMidiFile, FileTrackReader, DEFAULTS, phiExact, phiV, PH, HandModel, Accumulator };

// =====================================================================================
// CLI
// =====================================================================================
function parseArgs(argv) {
  const o = { files: [], csv: null, fps: 60, format: '{p},{ps},{mps}', decimals: 3, set: {}, json: false,
    timingOnly: false, overlaps: true, chanMask: 0xFFFF, entry: null, quiet: false, dumpInputs: null };
  for (let i = 0; i < argv.length; i++) {
    const a = argv[i];
    const val = () => { if (i + 1 >= argv.length) throw new Error('missing value for ' + a); return argv[++i]; };
    switch (a) {
      case '-h': case '--help': o.help = true; break;
      case '-c': case '--csv': o.csv = val(); break;
      case '--fps': o.fps = +val(); break;
      case '-f': case '--format': o.format = val(); break;
      case '-d': case '--decimals': o.decimals = +val(); break;
      case '--hours': o.set.hours = +val(); break;
      case '--restart': o.set.restart = +val(); break;
      case '--window': o.set.window = +val() / 1000; break;
      case '--timing-only': o.timingOnly = true; break;
      case '--overlaps': case '--swift': o.overlaps = true; break;
      case '--no-overlaps': o.overlaps = false; break;
      case '--no-velocity': o.velocity = false; break;
      case '--skip-channel': { for (const c of val().split(',')) o.chanMask &= ~(1 << ((+c - 1) & 15)); break; }
      case '--set': { const [k, v] = val().split('='); if (!(k in DEFAULTS)) throw new Error('unknown constant ' + k); o.set[k] = +v; break; }
      case '--json': o.json = true; break;
      case '-q': case '--quiet': o.quiet = true; break;
      case '-e': case '--entry': o.entry = val(); break;
      case '--decoder': o.decoder = val(); if (!['auto', 'builtin', 'external'].includes(o.decoder)) throw new Error('--decoder must be auto, builtin or external'); break;
      case '--ram': o.ramMB = +val(); break;
      case '--temp': o.tmpDir = val(); break;
      case '--dump-inputs': o.dumpInputs = val(); break;
      case '--version': o.version = true; break;
      default:
        if (a.startsWith('-') && a.length > 1) throw new Error('unknown option ' + a);
        o.files.push(a);
    }
  }
  return o;
}

const HELP = `midiprec ${VERSION} — raw precision (24h) of MIDI files on a 128-key piano

usage: node midiprec.js <file> [more files...] [options]

  -c, --csv <out.csv>       export precision over time (60 lines per second of MIDI)
  -f, --format "<tpl>"      CSV line template (default "{p},{ps},{mps}")
                            {p} precision  {ps} precision gained in the last second
                            {mps} max of {ps} so far  {t} time (s)  {f} frame
                            {n} notes so far  {nps} notes in the last second  {mnps} max NPS
  -d, --decimals <n>        decimals in the CSV (default 3)
      --fps <n>             CSV lines per second (default 60)
      --no-overlaps         merge overlapped notes (same key, same tick) into one
                            press. Default: each extra copy is a swift re-strike
                            of that key (like a GD swift/multi click)
      --no-velocity         ignore note velocity (default: window ÷ velocity/127,
                            so soft notes need less precision; vel 0 = note-off)
      --timing-only         NaNDL-style timing only (no piano physics multipliers)
      --hours <h>           completion target (default 24)
      --restart <s>         overhead per attempt in seconds (default 1)
      --window <ms>         max ± timing window for isolated notes (default 50)
      --skip-channel <list> ignore channels, e.g. 10 (drums)
      --set name=value      override a model constant (see README)
      --json                machine-readable output
  -q, --quiet               only print the result line

archives: .zip .7z .xz .gz .bz2 .zst .lzma .br .tar(.gz/.xz/...) .rar(needs 7-Zip)
  -e, --entry <text>        only analyse archive entries whose name contains <text>
      --decoder <mode>      auto (default: use installed 7-Zip/xz when faster),
                            builtin (pure JS), external (prefer 7-Zip)
      --ram <MB>            RAM for notes of compressed MIDIs before spilling to a
                            temp file (default 1024)
      --temp <dir>          folder for that temp file (default: system temp)
`;

function report(r, opts) {
  const { res, fin, mf } = r;
  const c = res.ctx;
  const tempoTxt = c.smpte ? 'SMPTE timing' : (c.tempoCount
    ? (c.tempoMin === c.tempoMax ? `${fmtBpm(c.tempoMin)} BPM` : `${fmtBpm(c.tempoMax)}–${fmtBpm(c.tempoMin)} BPM (${fmtNum(c.tempoCount)} tempo events)`)
    : '120 BPM (default)');
  if (opts.json) {
    return JSON.stringify({ file: r.file, precision: fin.L === Infinity ? 'inf' : fin.L, sigmaSeconds: fin.L > 0 ? 1 / fin.L : null,
      notes: res.notes, presses: res.presses, onsets: res.onsets, inputs: fin.inputs, duration: res.endTime,
      ppq: c.smpte ? null : c.ppq, tracks: mf.tracks, maxNps: res.acc.npsMax, overlapped: res.notes - res.presses, seconds: r.secs, reason: fin.reason || undefined });
  }
  if (opts.quiet) return `${fmtPrec(fin.L)}\t${r.file}`;
  const lines = [];
  const title = r.file.split(' › ').map((x, i) => i === 0 ? path.basename(x) : x).join(' › ');
  lines.push(`── ${title} ${'─'.repeat(Math.max(3, 60 - title.length))}`);
  lines.push(`  format ${mf.format} · ${mf.tracks} tracks · ${c.smpte ? 'SMPTE' : 'PPQ ' + c.ppq} · ${tempoTxt}`);
  lines.push(`  length ${fmtTime(res.endTime)} · ${fmtNum(res.notes)} notes · max ${fmtNum(res.acc.npsMax)} NPS`);
  lines.push(`  ${fmtNum(res.presses)} key presses in ${fmtNum(res.onsets)} chords/ticks → ${fmtNum(fin.inputs)} hand inputs`);
  const dup = res.notes - res.presses;
  if (dup > 0) lines.push(res.overlaps
    ? `  ${fmtNum(dup)} overlapped notes (same key & tick) counted as swift re-strikes`
    : `  ${fmtNum(dup)} overlapped notes (same key & tick) merged (--no-overlaps)`);
  lines.push('');
  lines.push(`  RAW PRECISION (24h): ${fmtPrec(fin.L)} σ/s${opts.timingOnly ? '  [timing only]' : ''}`);
  lines.push(`  = needs a timing consistency (1σ) of ${fmtSigma(fin.L)}`);
  if (fin.reason) lines.push(`  (${fin.reason})`);
  if (r.csvInfo) lines.push(`  CSV: ${r.csvPath} (${fmtNum(r.csvInfo.lines)} lines, max P/S ${fmtPrec(r.csvInfo.maxPS)})`);
  if (r.compactMB !== undefined && r.compactMB > 64) lines.push(`  notes buffer ${r.compactMB.toFixed(0)} MB${r.spilledMB ? ` (${r.spilledMB.toFixed(0)} MB spilled to temp file)` : ''}`);
  const mnps = res.notes / Math.max(r.secs, 1e-9);
  lines.push(`  analysed in ${r.secs.toFixed(2)} s (${(mnps / 1e6).toFixed(1)} M notes/s)${r.csvInfo ? `, CSV ${r.csvSecs.toFixed(2)} s` : ''}`);
  return lines.join('\n');
}

function main() {
  let opts;
  try { opts = parseArgs(process.argv.slice(2)); } catch (e) { console.error(e.message); console.error(HELP); process.exit(2); }
  if (opts.version) { console.log(VERSION); return; }
  if (opts.help || opts.files.length === 0) { console.log(HELP); return; }
  const P = Object.assign({}, DEFAULTS, opts.set);
  const results = [];
  for (const f of opts.files) {
    try {
      if (!fs.existsSync(f)) throw new Error('file not found');
      if (opts.dumpInputs) dumpInputs(f, opts, P);
      const rs = analyzeInput(f, opts, P);
      for (const r of rs) { results.push(r); console.log(report(r, opts)); }
    } catch (e) {
      console.error(`${f}: ${e && e.message ? e.message : e}`);
      if (process.env.MIDIPREC_DEBUG) console.error(e.stack);
      process.exitCode = 1;
    }
  }
  if (process.env.MIDIPREC_DEBUG) console.error(`[debug] peak RSS ${(process.resourceUsage().maxRSS / 1024).toFixed(0)} MB`);
  if (results.length > 1 && !opts.json && !opts.quiet) {
    console.log('\nRanking (hardest first):');
    results.sort((a, b) => b.fin.L - a.fin.L);
    results.forEach((r, i) => console.log(`  ${String(i + 1).padStart(2)}. ${fmtPrec(r.fin.L).padStart(22)} σ/s  ${r.file.split(' › ').map((x, j) => j === 0 ? path.basename(x) : x).join(' › ')}`));
  }
}

// debug: re-run and dump (t, a) of every input as float64 pairs
function dumpInputs(file, opts, P) {
  const out = [];
  const mf = openMidiFile(file);
  const readers = mf.tracks.map((t) => new FileTrackReader(mf.fd, t.start, t.length, 1 << 16, opts.chanMask));
  const o1 = Accumulator.prototype.addInput, o2 = Accumulator.prototype.addInputN;
  Accumulator.prototype.addInput = function (t, a) { out.push(t, a); return o1.call(this, t, a); };
  Accumulator.prototype.addInputN = function (t, a, n) { for (let i = 0; i < n; i++) out.push(t, a); return o2.call(this, t, a, n); };
  try { analyzeReaders(readers, mf.division, P, Object.assign({}, opts, { csv: null })); }
  finally { Accumulator.prototype.addInput = o1; Accumulator.prototype.addInputN = o2; fs.closeSync(mf.fd); }
  fs.writeFileSync(opts.dumpInputs, Buffer.from(new Float64Array(out).buffer));
}

if (!WT.isMainThread) { if (WT.workerData && WT.workerData.midiprecJob) workerMain(); }
else if (require.main === module) main();
