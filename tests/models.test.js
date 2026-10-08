import { test } from 'node:test';
import assert from 'node:assert/strict';
import * as M from '../web/models.js';

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg ?? ''} expected ${b}, got ${a}`);

test('special functions: logGamma, normal CDF and inverse, bivariate normal', () => {
  close(Math.exp(M.logGamma(5)), 24, 1e-9, 'Gamma(5)');
  close(M.normCdf(1.96), 0.9750021, 2e-6, 'Phi(1.96)');
  close(M.normInv(0.975), 1.959964, 1e-5, 'Phi^-1(0.975)');
  // P(X<0, Y<0) = 1/4 + asin(rho)/(2 pi)
  close(M.bivariateNormCdf(0, 0, 0.5), 0.25 + Math.asin(0.5) / (2 * Math.PI), 1e-4, 'bivariate at origin');
});

test('Poisson partition: size bias doubles the multi-founder risk at low density', () => {
  const r = M.poissonPartition({ f0: 0.7, phi: 0.05 });
  close(r.lambda, 0.35667, 1e-4); close(r.pMulti, 0.16776, 1e-4); close(r.pMultiSelected, 0.29374, 1e-4);
  const low = M.poissonPartition({ f0: 0.999, phi: 1e-4 });
  close(low.ratio, 2, 0.01, 'ratio tends to 2');
});

test('MOI: multi-copy share among scored cells exceeds the share among transduced cells', () => {
  const r = M.moi({ m: 1, fi: 0.5, phi: 0.05 });
  close(r.p0, Math.exp(-0.5), 1e-12); assert.ok(r.pMultiAmongScored > r.pMultiAmongTransduced);
  close(r.pk.reduce((a, b) => a + b, 0), 1, 1e-9, 'pmf sums to 1');
});

test('Kärber design: SE shrinks as 1/sqrt(n-1) and grows with the dilution step', () => {
  const a = M.karberDesign({ n: 4 }), b = M.karberDesign({ n: 13 });
  close(a.seLog10 / b.seLog10, Math.sqrt(12 / 3), 1e-9);
  assert.ok(M.karberDesign({ step: 10 }).seLog10 > M.karberDesign({ step: 2, levels: 20, logStart: 2 }).seLog10);
});

test('Low-copy detection: LoD95 is three copies before any chemistry', () => {
  const r = M.lowCopy({ copiesPerReaction: 3 });
  close(r.pDetect, 1 - Math.exp(-3), 1e-12); close(r.lod95, 2.9957, 1e-3);
});

test('Bottlenecks: neutral resampling CV matches sqrt(k / coverage) and loss matches e^-c for one stage', () => {
  const r = M.bottlenecks({ coverage: 30, k: 3 });
  close(r.cvAfterK, Math.sqrt(3 / 30), 0.01);
  const one = M.bottlenecks({ coverage: 3, k: 1 });
  close(one.lostAfterK, Math.exp(-3), 1e-6);
});

test('PCR branching: closed-form CV matches a Galton–Watson simulation', () => {
  const e = 0.7, n = 12, R = M.rng(7); const ys = [];
  for (let s = 0; s < 4000; s++) { let z = 1; for (let c = 0; c < n; c++) z += R.binomial(z, e); ys.push(z); }
  const mean = ys.reduce((a, b) => a + b) / ys.length, v = ys.reduce((a, b) => a + (b - mean) ** 2, 0) / (ys.length - 1);
  const r = M.pcrBranching({ efficiency: e, cycles: n, templates: 1 });
  close(mean, r.meanAmplification, 0.05 * r.meanAmplification, 'mean'); close(Math.sqrt(v) / mean, r.cvSingleTemplate, 0.03, 'cv');
});

test('UMI collisions: expected distinct labels follow U(1 - e^{-N/U})', () => {
  const r = M.pcrBranching({ efficiency: 0.9, cycles: 20, umis: 1000, molecules: 500 });
  close(r.umiCollisionFraction, 1 - 1000 * (1 - Math.exp(-0.5)) / 500, 1e-12);
});

test('Allele dropout and capture dropout closed forms', () => {
  const a = M.alleleDropout({ q: 0.8 }); close(a.adoAmongAmplified, 0.32 / 0.96, 1e-12);
  const c = M.captureDropout({ trueCount: 5, capture: 0.1, trueCv: 0 }); close(c.pZero, 0.9 ** 5, 1e-12); close(c.cvObserved, Math.sqrt(0.9 / 0.5), 1e-12);
});

test('Shot noise and bleaching steps', () => {
  close(M.shotNoise({ photons: 100, readNoise: 0 }).snr, 10, 1e-12);
  const b = M.bleachingSteps({ subunits: 4, pLabel: 0.75 }); close(b.dist.reduce((x, y) => x + y, 0), 1, 1e-12); close(b.pFullCount, 0.75 ** 4, 1e-12); assert.equal(b.modeSteps, 3);
});

test('Number fluctuations and Berg–Purcell scale as expected', () => {
  const f = M.numberFluctuations({ concentration_nM: 1, volume_fL: 1 }); close(f.meanMolecules, 0.602, 1e-3);
  const a = M.bergPurcell({ radius_um: 1, concentration_nM: 1, D: 100, T: 1 }), b = M.bergPurcell({ radius_um: 1, concentration_nM: 1, D: 100, T: 4 });
  close(a.relativeError / b.relativeError, 2, 1e-9, 'halves with 4x time');
});

test('Channel noise: variance is a parabola in the mean with apex at p = 1/2', () => {
  const r = M.channelNoise({ N: 100, i: 1, pOpen: 0.5 }); close(r.variance, 25, 1e-9);
  const best = r.parabola.reduce((m, x) => x.variance > m.variance ? x : m); close(best.mean, 50, 1.5);
});

test('Bursting: negative binomial has mean fb and Fano 1 + b; pmf sums to one', () => {
  let s = 0, m = 0, m2 = 0; for (let k = 0; k < 3000; k++) { const p = M.nbPmf(k, 2, 20); s += p; m += k * p; m2 += k * k * p; }
  close(s, 1, 1e-9); close(m, 40, 1e-6); close((m2 - m * m) / m, 21, 1e-5);
});

test('Sort relaxation: enrichment decays as e^{-t/tau} and the sorted mean matches the truncated normal', () => {
  const r = M.sortRelaxation({ sigma: 0.5, topFraction: 0.05, tau: 2, t: 2 });
  close(r.enrichmentRemaining, Math.exp(-1), 1e-12); close(r.sortedMeanSd, 2.0627, 1e-3);
});

test('Fractional killing: two-dose survival matches a bivariate simulation', () => {
  const R = M.rng(3), sigma = 1, d50 = 1, dose = 3, rho = 0.8, n = 60000; let s1 = 0, both = 0;
  for (let i = 0; i < n; i++) { const z1 = R.normal(), z2 = rho * z1 + Math.sqrt(1 - rho * rho) * R.normal(); const zc = Math.log(dose / d50) / sigma; if (z1 > zc) { s1++; if (z2 > zc) both++; } }
  const r = M.fractionalKilling({ sigma, d50, dose, rho });
  close(s1 / n, r.survival1, 0.006); close(both / s1, r.survival2Given, 0.012); assert.ok(r.survival2Given > r.survival1);
});

test('Biphasic kill and switching asymptotes', () => {
  const k = M.biphasicKill({ f: 0.001, ks: 2, kp: 0.01, t: 6 }); close(k.survivingFraction, 0.999 * Math.exp(-12) + 0.001 * Math.exp(-0.06), 1e-12); assert.ok(k.timeTo999 > k.timeTo999NoPersisters);
  const s = M.switching({ kOff: 0.1, kOn: 0.02, generations: 1e4 }); close(s.fractionOn, 0.02 / 0.12, 1e-9);
});

test('Luria–Delbrück: p0 = e^-m, p1 = m e^-m / 2, and the distribution is heavier-tailed than Poisson', () => {
  const r = M.luriaDelbruck({ m: 2, rmax: 400 });
  close(r.p0, Math.exp(-2), 1e-12); close(r.p1, Math.exp(-2), 1e-12); // m e^-m / 2 = e^-2 for m = 2
  close(r.mass, 1, 0.02, 'mass with truncation');
  assert.ok(r.tailAbove(20) > 1 - M.poissonCdf(20, r.meanTruncated));
});

test('Clonal drift: without selection the variance follows p(1-p)(1 - (1 - 1/N)^k)', () => {
  const r = M.clonalDrift({ p0: 0.2, N: 500, s: 0, passages: 10 });
  close(r.finalSd ** 2, r.driftOnlyVar, 0.02 * r.driftOnlyVar); close(r.final, 0.2, 1e-12);
  assert.ok(M.clonalDrift({ p0: 0.2, N: 500, s: 0.1, passages: 10 }).final > 0.2);
});

test('Mitochondrial bottleneck and initiating-cell frequency information', () => {
  close(M.mtBottleneck({ h: 0.3, n: 10 }).variance, 0.021, 1e-12);
  const a = M.initiatingFrequency({ frequency: 1e-3, doses: [100, 300, 1000, 3000], n: 6 }), b = M.initiatingFrequency({ frequency: 1e-3, doses: [100, 300, 1000, 3000], n: 24 });
  close(a.seLogF / b.seLogF, 2, 1e-9, 'SE halves with 4x replicates');
});

test('Hierarchical SE: components sum and a single day is dominated by day variance', () => {
  const r = M.hierarchicalSe({ sdWell: 0.1, sdPlate: 0.1, sdDay: 0.15, nWells: 8, nPlates: 2, nDays: 3 });
  close(r.shareDay + r.sharePlate + r.shareWell, 1, 1e-12); assert.ok(r.se > r.seNaive);
});

test('Spheroid hook: secondary limitation alone gives a plateau; carryover, crowding and quenching give a fall; the direct label is monotonic', () => {
  const C = M.logspace(0.1, 1000, 25);
  const mono = (ys) => ys.every((y, i) => i === 0 || y >= ys[i - 1] - 1e-9);
  const ideal = M.hookCurve({ gamma: 0, rhoQ: 1e12, eps: 0, mu: 0, limitSecondary: false }, C).map(x => x.signal);
  const limited = M.hookCurve({ gamma: 0, rhoQ: 1e12, eps: 0, mu: 0, limitSecondary: true }, C).map(x => x.signal);
  const carry = M.hookCurve({ gamma: 0, rhoQ: 1e12, eps: 0, mu: 0, limitSecondary: true, carry: 0.02 }, C).map(x => x.signal);
  const crowd = M.hookCurve({ gamma: 0.75, rhoQ: 1e12, eps: 0, mu: 0, limitSecondary: false }, C).map(x => x.signal);
  const quench = M.hookCurve({ gamma: 0, rhoQ: 2000, eps: 0, mu: 0, limitSecondary: false }, C).map(x => x.signal);
  const direct = M.hookCurve({}, C).map(x => x.directLabel);
  assert.ok(mono(ideal) && mono(limited) && mono(direct), 'ideal, limited and direct are monotonic');
  assert.ok(limited[24] < ideal[24], 'limitation caps the signal');
  for (const ys of [carry, crowd, quench]) assert.ok(Math.max(...ys) > ys[24] * 1.02, 'hook present');
});

test('Wicksell: mean section radius is pi R / 4 and the fraction of fully stained-looking sections is t / R', () => {
  const r = M.wicksell({ R: 250, t: 50, n: 20000 });
  close(r.meanApparentRadius, Math.PI * 250 / 4, 0.05); close(r.fractionLookFullyStained, 0.2, 1e-3); assert.ok(r.biasFactor > 1);
  // Monte Carlo check of the apparent thickness
  const R = M.rng(11); let s = 0; const n = 200000;
  for (let i = 0; i < n; i++) { const z = R() * 250, Rp = Math.sqrt(250 ** 2 - z * z), inner = 200 ** 2 - z * z; s += inner > 0 ? Rp - Math.sqrt(inner) : Rp; }
  close(s / n, r.meanApparentThickness, 0.3);
});

test('Seeding size: radius CV is one third of the count CV', () => {
  const r = M.seedingSize({ meanCells: 400 }); close(r.cvRadius, 0.05 / 3, 1e-12);
});
