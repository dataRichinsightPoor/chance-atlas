// Chance Atlas · model engine
// Closed-form and small numerical models for stochastic processes that shape bioassay outcomes.
// Dependency-free ES module. Every function here is covered by tests/models.test.js.
// All examples are synthetic. Nothing here is validated against any experimental system.

export const VERSION = '0.1.0-alpha';

/* ----------------------------------------------------------------- numerics */

export function logGamma(x) {
  // Lanczos approximation (g = 7, n = 9), relative error < 1e-13 for x > 0.
  const c = [0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7];
  if (x < 0.5) return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1 - x);
  x -= 1;
  let a = c[0];
  const t = x + 7.5;
  for (let i = 1; i < 9; i++) a += c[i] / (x + i);
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}
export const logChoose = (n, k) => logGamma(n + 1) - logGamma(k + 1) - logGamma(n - k + 1);
export const poissonPmf = (k, lam) => lam === 0 ? (k === 0 ? 1 : 0) : Math.exp(k * Math.log(lam) - lam - logGamma(k + 1));
export const binomPmf = (k, n, p) => (p <= 0) ? (k === 0 ? 1 : 0) : (p >= 1) ? (k === n ? 1 : 0) : Math.exp(logChoose(n, k) + k * Math.log(p) + (n - k) * Math.log(1 - p));
export function poissonCdf(k, lam) { let s = 0; for (let i = 0; i <= k; i++) s += poissonPmf(i, lam); return Math.min(1, s); }
export function binomCdf(k, n, p) { let s = 0; for (let i = 0; i <= k; i++) s += binomPmf(i, n, p); return Math.min(1, s); }

export function erf(x) {
  // Abramowitz & Stegun 7.1.26 refined by one Newton step on the complementary error function series is unnecessary here:
  // use a high-precision rational approximation (W. J. Cody style), |error| < 1.2e-7.
  const t = 1 / (1 + 0.3275911 * Math.abs(x));
  const y = 1 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t - 0.284496736) * t + 0.254829592) * t * Math.exp(-x * x);
  return x >= 0 ? y : -y;
}
export const normCdf = (z) => 0.5 * (1 + erf(z / Math.SQRT2));
export const normPdf = (z) => Math.exp(-0.5 * z * z) / Math.sqrt(2 * Math.PI);
export function normInv(p) {
  // Acklam's algorithm with one Newton refinement.
  if (p <= 0) return -Infinity; if (p >= 1) return Infinity;
  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02, 1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02, 6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00, -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00, 3.754408661907416e+00];
  const pl = 0.02425, ph = 1 - pl; let q, r, x;
  if (p < pl) { q = Math.sqrt(-2 * Math.log(p)); x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  else if (p <= ph) { q = p - 0.5; r = q * q; x = (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1); }
  else { q = Math.sqrt(-2 * Math.log(1 - p)); x = -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1); }
  const e = normCdf(x) - p; x -= e / normPdf(x);
  return x;
}
export function bivariateNormCdf(a, b, rho, n = 400) {
  // P(X < a, Y < b) for standard bivariate normal with correlation rho, by integrating over X.
  if (Math.abs(rho) >= 1) rho = Math.sign(rho) * 0.999999;
  const lo = -8, hi = Math.min(a, 8);
  if (hi <= lo) return 0;
  const h = (hi - lo) / n, s = Math.sqrt(1 - rho * rho);
  let sum = 0;
  for (let i = 0; i <= n; i++) {
    const x = lo + i * h, w = (i === 0 || i === n) ? 0.5 : 1;
    sum += w * normPdf(x) * normCdf((b - rho * x) / s);
  }
  return sum * h;
}
export function linspace(a, b, n) { const o = []; for (let i = 0; i < n; i++) o.push(a + (b - a) * i / (n - 1)); return o; }
export function logspace(a, b, n) { return linspace(Math.log10(a), Math.log10(b), n).map(x => 10 ** x); }

// Small deterministic PRNG (mulberry32) so simulations in tests are reproducible.
export function rng(seed = 1) {
  let a = seed >>> 0;
  const u = () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; };
  u.normal = () => { const r = Math.sqrt(-2 * Math.log(1 - u())); return r * Math.cos(2 * Math.PI * u()); };
  u.poisson = (lam) => { if (lam > 50) return Math.max(0, Math.round(lam + Math.sqrt(lam) * u.normal())); const L = Math.exp(-lam); let k = 0, p = 1; do { k++; p *= u(); } while (p > L); return k - 1; };
  u.binomial = (n, p) => { if (n > 200) return Math.max(0, Math.min(n, Math.round(n * p + Math.sqrt(n * p * (1 - p)) * u.normal()))); let k = 0; for (let i = 0; i < n; i++) if (u() < p) k++; return k; };
  return u;
}

/* ------------------------------------------------- A. sampling and partitioning */

// Poisson partitions. f0 = empty fraction. Returns occupancy statistics and selection size bias.
export function poissonPartition({ f0, phi = 1 }) {
  const lam = -Math.log(f0);
  const pMulti = 1 - lam * Math.exp(-lam) / (1 - Math.exp(-lam));           // P(k >= 2 | k >= 1)
  const pMonoSel = phi * lam * Math.exp(-lam) / (1 - Math.exp(-phi * lam));  // P(k = 1 | passes screen)
  return { lambda: lam, pMulti, pMultiSelected: 1 - pMonoSel, rareLimit: 1 - Math.exp(-lam), ratio: (1 - pMonoSel) / pMulti };
}

// Multiplicity of infection. m = particles per cell, fi = infectious fraction, phi = per-copy chance a copy scores.
export function moi({ m, fi = 1, phi = 0.05, kmax = 12 }) {
  const lam = m * fi; const pk = [];
  for (let k = 0; k <= kmax; k++) pk.push(poissonPmf(k, lam));
  const p0 = pk[0], p1 = pk[1], pMulti = 1 - p0 - p1;
  const pMultiAmongTransduced = pMulti / (1 - p0);
  let num = 0, den = 0, meanScored = 0;
  for (let k = 1; k <= kmax; k++) { const w = 1 - (1 - phi) ** k; den += w * pk[k]; if (k >= 2) num += w * pk[k]; meanScored += k * w * pk[k]; }
  return { lambda: lam, p0, p1, pMulti, pMultiAmongTransduced, pMultiAmongScored: num / den, meanCopiesTransduced: lam / (1 - p0), meanCopiesScored: meanScored / den, pk };
}

// Endpoint titration (Spearman–Kärber). dilutions: array of doses (infectious units per well); n wells each.
// Expected response p_i = 1 - exp(-dose_i). Returns the Kärber variance of log10 titer under the design.
export function karberDesign({ logStart = 2, step = 10, levels = 8, n = 8, titerLog10 = 4.3 }) {
  // doses in TCID-like units: dose_i = 10^(titerLog10 - level_i) with level_i = logStart..., infectious units per well in natural log terms: convert to Poisson mean by ln2 (one TCID50 = ln 2 units)
  const rows = [];
  let varSum = 0;
  for (let i = 0; i < levels; i++) {
    const dil = logStart + i * Math.log10(step);
    const mean = Math.LN2 * 10 ** (titerLog10 - dil);
    const p = 1 - Math.exp(-mean);
    varSum += p * (1 - p) / Math.max(1, n - 1);
    rows.push({ dilutionLog10: dil, pExpected: p });
  }
  const se = Math.log10(step) * Math.sqrt(varSum);
  return { rows, seLog10: se, ci95FactorLow: 10 ** (-1.96 * se), ci95FactorHigh: 10 ** (1.96 * se) };
}

// Low-copy detection. c = mean copies per reaction after extraction efficiency.
export function lowCopy({ copiesPerReaction, efficiency = 1, replicates = 1 }) {
  const c = copiesPerReaction * efficiency;
  const pDetect = 1 - Math.exp(-c);
  return { effectiveCopies: c, pDetect, pDetectAny: 1 - (1 - pDetect) ** replicates, pAllReplicates: pDetect ** replicates, lod95: -Math.log(0.05) / efficiency, lod99: -Math.log(0.01) / efficiency };
}

// Counting statistics: CV from Poisson counts; events needed for a target CV.
export function countingCv({ events, targetCv = 0.05 }) {
  return { cv: 1 / Math.sqrt(events), eventsForTarget: Math.ceil(1 / targetCv ** 2), ci95Low: events - 1.96 * Math.sqrt(events), ci95High: events + 1.96 * Math.sqrt(events) };
}

// Library bottlenecks. L elements, coverage c cells per element at each of k bottlenecks.
// Each bottleneck resamples counts ~ Poisson(c * current relative abundance) approximately.
export function bottlenecks({ coverage, k = 3, foldThreshold = 2, elements = 1 }) {
  // Fraction of elements lost at one Poisson bottleneck of mean c: e^-c. After k independent bottlenecks of a neutral element,
  // survival requires count >= 1 at each stage; the count process is a Galton–Watson with Poisson offspring (mean 1 after renormalization).
  // Compute by iterating the pgf on a truncated support.
  const K = Math.max(60, Math.ceil(coverage * 8));
  let dist = new Array(K + 1).fill(0); dist[Math.min(K, Math.round(coverage))] = 1; // start at nominal count
  const mean0 = coverage;
  for (let s = 0; s < k; s++) {
    const next = new Array(K + 1).fill(0);
    for (let i = 0; i <= K; i++) if (dist[i] > 0) {
      const lam = i; // neutral: expected count after resampling equals current count (library kept at fixed coverage)
      if (lam === 0) { next[0] += dist[i]; continue; }
      for (let j = 0; j <= K; j++) next[j] += dist[i] * poissonPmf(j, lam);
    }
    dist = next;
  }
  const lost = dist[0];
  let depleted = 0, enriched = 0, mean = 0, m2 = 0;
  for (let j = 0; j <= K; j++) { mean += j * dist[j]; m2 += j * j * dist[j]; if (j <= mean0 / foldThreshold) depleted += dist[j]; if (j >= mean0 * foldThreshold) enriched += dist[j]; }
  const varc = m2 - mean * mean;
  return { lostAfterK: lost, lostOne: Math.exp(-coverage), falseDepleted: depleted, falseEnriched: enriched, cvAfterK: Math.sqrt(varc) / mean, cvOneTheory: Math.sqrt(k / coverage), expectedLostElements: lost * elements, dist };
}

// Undersampling: molecules observed at least once at mean depth per molecule m (Lander–Waterman style).
export function undersampling({ depthPerMolecule, clonotypeFrequency = 1e-4, reads = 1e5 }) {
  const seen = 1 - Math.exp(-depthPerMolecule);
  return { fractionSeen: seen, duplicateRate: 1 - seen / depthPerMolecule, pClonotypeSeen: 1 - Math.exp(-clonotypeFrequency * reads), readsFor95: -Math.log(0.05) / clonotypeFrequency };
}

/* ----------------------------------------------------------- B. amplification */

// PCR as a Galton–Watson branching process with per-cycle efficiency e.
export function pcrBranching({ efficiency, cycles = 30, templates = 1, umis = 4 ** 10, molecules = 1e5 }) {
  const e = efficiency, n = cycles;
  const mean = (1 + e) ** n;
  // Variance of the yield from a single template: Var_n = (1-e)/(1+e) * mean^2 * (1 - mean^-1)  [exact for Bernoulli duplication]
  const cv2one = (1 - e) / (1 + e) * (1 - 1 / mean);
  const cv = Math.sqrt(cv2one / templates);
  const expectedDistinct = umis * (1 - Math.exp(-molecules / umis));
  return { meanAmplification: mean, cvSingleTemplate: Math.sqrt(cv2one), cv, umiCollisionFraction: 1 - expectedDistinct / molecules, umisForOnePct: molecules / 0.0201 };
}

// Allele dropout from two templates each captured with probability q.
export function alleleDropout({ q }) {
  const both = q * q, one = 2 * q * (1 - q), none = (1 - q) ** 2;
  return { pBoth: both, pOne: one, pNone: none, adoAmongAmplified: one / (1 - none) };
}

// Capture and dropout in single-cell counting. n true molecules, capture probability cap.
export function captureDropout({ trueCount, capture, trueCv = 0.3 }) {
  const pZero = (1 - capture) ** trueCount;
  const meanObs = trueCount * capture;
  // CV^2 of observed = CV^2_true + (1 - cap) / (cap * mean_true)  (binomial thinning of a variable count)
  const cv2 = trueCv ** 2 + (1 - capture) / (capture * trueCount);
  return { pZero, meanObserved: meanObs, cvObserved: Math.sqrt(cv2), cvTrue: trueCv, noiseShare: ((1 - capture) / (capture * trueCount)) / cv2 };
}

/* ------------------------------------------------------- C. detection physics */

export function shotNoise({ photons, readNoise = 1.5, threshold = null }) {
  const snr = photons / Math.sqrt(photons + readNoise ** 2);
  const T = threshold ?? photons / 2;
  const pBelow = poissonCdf(Math.floor(T), photons);
  return { snr, cv: 1 / snr, pBelowThreshold: pBelow, photonsForSnr10: 100 + readNoise ** 2 * 0 + 10 * Math.sqrt(100 + readNoise ** 2) };
}

// Photobleaching step counts with labeling/maturation probability pl.
export function bleachingSteps({ subunits, pLabel, nmax = 6 }) {
  const dist = []; for (let k = 0; k <= subunits; k++) dist.push(binomPmf(k, subunits, pLabel));
  const meanSteps = subunits * pLabel;
  const pFull = pLabel ** subunits;
  // naive estimate of subunit number from the modal step count
  let mode = 0; for (let k = 1; k <= subunits; k++) if (dist[k] > dist[mode]) mode = k;
  return { dist, meanSteps, pFullCount: pFull, modeSteps: mode, pZero: dist[0], pLabelFromMean: (m) => m / subunits };
}

// Number fluctuations in a detection volume (FCS). c in nM, volume in femtoliters.
export function numberFluctuations({ concentration_nM, volume_fL }) {
  const N = concentration_nM * 1e-9 * volume_fL * 1e-15 * 6.02214076e23;
  return { meanMolecules: N, relativeVariance: 1 / N, cv: 1 / Math.sqrt(N), g0: 1 / N };
}

// Berg–Purcell limit for a perfectly absorbing sphere of radius a (um) in concentration c (nM), diffusion D (um^2/s), time T (s).
export function bergPurcell({ radius_um, concentration_nM, D = 100, T }) {
  const cPerUm3 = concentration_nM * 1e-9 * 6.02214076e23 * 1e-15; // molecules per um^3 (1 um^3 = 1e-15 L)
  const rel = 1 / Math.sqrt(D * radius_um * cPerUm3 * T);
  const encounters = 4 * Math.PI * D * radius_um * cPerUm3 * T;
  return { relativeError: rel, encounters, timeFor5pct: 1 / (D * radius_um * cPerUm3 * 0.05 ** 2) };
}

// Single-channel binomial fluctuation (non-stationary noise analysis). N channels, unitary current i (pA), open probability p.
export function channelNoise({ N, i, pOpen }) {
  const mean = N * i * pOpen, variance = i * i * N * pOpen * (1 - pOpen);
  const ps = linspace(0.001, 0.999, 60);
  return { meanCurrent: mean, variance, cv: Math.sqrt(variance) / mean, parabola: ps.map(p => ({ mean: N * i * p, variance: i * i * N * p * (1 - p) })) };
}

/* ------------------------------------------------------ D. single-cell biology */

// Negative binomial pmf for a bursty telegraph model: mean = f*b, Fano = 1 + b, r = f (bursts per lifetime), p = b/(1+b).
export function nbPmf(k, r, b) { const p = b / (1 + b); return Math.exp(logGamma(k + r) - logGamma(k + 1) - logGamma(r) + r * Math.log(1 - p) + k * Math.log(p)); }
export function bursting({ frequency, size, threshold, kmax = 2000 }) {
  const mean = frequency * size, fano = 1 + size;
  const above = (f, b) => { let s = 0; for (let k = 0; k < Math.ceil(threshold); k++) s += nbPmf(k, f, b); return 1 - s; };
  return { mean, fano, cv: Math.sqrt(fano / mean), pctPositive: above(frequency, size), pctPositiveDoubleFreq: above(2 * frequency, size), pctPositiveDoubleSize: above(frequency, 2 * size) };
}

// Regression to the mean after sorting. Log-protein state ~ N(0, sigma^2); sort the top fraction q; OU relaxation with mixing time tau.
export function sortRelaxation({ sigma, topFraction, tau, t }) {
  const z = normInv(1 - topFraction);
  const meanSortedZ = normPdf(z) / topFraction;             // E[Z | Z > z]
  const varSortedZ = 1 + z * meanSortedZ - meanSortedZ ** 2;  // Var[Z | Z > z]
  const decay = Math.exp(-t / tau);
  const meanT = meanSortedZ * decay;                         // in sd units
  const varT = varSortedZ * decay ** 2 + (1 - decay ** 2);
  return { sortedMeanSd: meanSortedZ, meanAtT_sd: meanT, enrichmentRemaining: decay, foldMeanLinear: Math.exp(sigma * meanT + 0.5 * sigma ** 2 * varT) / Math.exp(0.5 * sigma ** 2), sdAtT: Math.sqrt(varT), zCut: z };
}

// Fractional killing from a lognormal protective state. A cell dies if dose * exp(sigma Z) < threshold... parameterize as:
// survive if Z > zc(dose) where zc = (ln(dose/d50)) / sigma. Second identical dose: states correlated rho between doses.
export function fractionalKilling({ sigma, d50, dose, rho }) {
  const zc = Math.log(dose / d50) / sigma;
  const s1 = 1 - normCdf(zc);
  // P(survive both) = P(Z1 > zc, Z2 > zc) = 1 - 2 Phi(zc) + Phi2(zc, zc, rho)
  const both = 1 - 2 * normCdf(zc) + bivariateNormCdf(zc, zc, rho);
  const s2given = both / s1;
  return { survival1: s1, survival2Given: s2given, survivalTwoDoses: both, naiveTwoDoses: s1 * s1, foldLessKilling: (1 - s1) / (1 - s2given), zc };
}

// Persisters: biphasic kill with switching. Fractions normal (1-f) and persister f; kill rates ks, kp; switching ignored over short exposures unless a > 0.
export function biphasicKill({ f, ks, kp, t }) {
  const N = (1 - f) * Math.exp(-ks * t) + f * Math.exp(-kp * t);
  // time to 99.9% kill
  let t999 = NaN; const lo = 0, hi = 1e4; let a = lo, b = hi;
  const g = (x) => (1 - f) * Math.exp(-ks * x) + f * Math.exp(-kp * x) - 1e-3;
  if (g(hi) < 0) { for (let i = 0; i < 80; i++) { const m = 0.5 * (a + b); if (g(m) > 0) a = m; else b = m; } t999 = 0.5 * (a + b); }
  return { survivingFraction: N, log10Kill: -Math.log10(N), timeTo999: t999, timeTo999NoPersisters: Math.log(1000) / ks, persisterShareOfSurvivors: f * Math.exp(-kp * t) / N };
}

// Two-state epigenetic switching (silencing). kOff = silencing rate per generation, kOn = reactivation rate.
export function switching({ kOff, kOn, generations, initialOn = 1 }) {
  const s = kOff + kOn, eq = kOn / s;
  const on = eq + (initialOn - eq) * Math.exp(-s * generations);
  return { fractionOn: on, equilibriumOn: eq, relaxationGenerations: 1 / s, halfLifeGenerations: Math.log(2) / s };
}

// Growth variability: clone size after T doublings-worth of time when per-division time has CV cvT. Lognormal approximation.
export function cloneSize({ meanDivisions, cvDivisionTime }) {
  // ln(size) = ln2 * (number of divisions). Number of divisions in fixed time ~ Normal(meanDivisions, var = meanDivisions * cvT^2) (renewal CLT).
  const sdDiv = Math.sqrt(meanDivisions) * cvDivisionTime;
  const sdLn = Math.LN2 * sdDiv, meanLn = Math.LN2 * meanDivisions;
  const cv = Math.sqrt(Math.exp(sdLn ** 2) - 1);
  return { medianSize: Math.exp(meanLn), cvSize: cv, p10: Math.exp(meanLn - 1.2816 * sdLn), p90: Math.exp(meanLn + 1.2816 * sdLn), foldP90P10: Math.exp(2 * 1.2816 * sdLn) };
}

/* ------------------------------------------------- E. population and evolution */

// Luria–Delbrück (Lea–Coulson) distribution by the Ma–Sandri–Sarkar recursion. m = expected mutations per culture.
export function luriaDelbruck({ m, rmax = 200 }) {
  const p = [Math.exp(-m)];
  for (let r = 1; r <= rmax; r++) { let s = 0; for (let i = 0; i < r; i++) s += p[i] / (r - i + 1); p.push(m / r * s); }
  const mass = p.reduce((a, b) => a + b, 0);
  let median = 0, c = 0; for (let r = 0; r <= rmax; r++) { c += p[r]; if (c >= 0.5) { median = r; break; } }
  // Poisson with the same median for comparison; and mean of the truncated distribution
  let mean = 0; for (let r = 0; r <= rmax; r++) mean += r * p[r];
  const p0Estimate = (f0) => -Math.log(f0);
  return { pmf: p, mass, median, meanTruncated: mean, p0: p[0], p1: p[1], tailAbove: (r) => 1 - p.slice(0, r + 1).reduce((a, b) => a + b, 0), p0Estimate };
}

// Clonal drift with selection across passages. p0 initial subclone frequency, N effective cells carried per passage, s selection coefficient per passage.
export function clonalDrift({ p0, N, s, passages }) {
  const traj = [], sd = [];
  let p = p0, v = 0;
  for (let k = 0; k <= passages; k++) {
    traj.push(p); sd.push(Math.sqrt(v));
    const pn = p * Math.exp(s) / (p * Math.exp(s) + (1 - p));
    // variance propagation: Var_{k+1} ≈ Var_k * (dp'/dp)^2 + p'(1-p')/N
    const dp = Math.exp(s) / (p * Math.exp(s) + (1 - p)) ** 2;
    v = v * dp * dp + pn * (1 - pn) / N;
    p = pn;
  }
  return { trajectory: traj, sd, final: p, finalSd: sd[passages], driftOnlyVar: p0 * (1 - p0) * (1 - (1 - 1 / N) ** passages) };
}

// Mitochondrial bottleneck. Heteroplasmy h, effective segregating units n.
export function mtBottleneck({ h, n, generationsOfDrift = 1 }) {
  const varOne = h * (1 - h) / n;
  const varK = h * (1 - h) * (1 - (1 - 1 / n) ** generationsOfDrift);
  const sd = Math.sqrt(varK);
  return { variance: varK, sd, pAbove: (thr) => 1 - normCdf((thr - h) / sd), pBelow: (thr) => normCdf((thr - h) / sd), varianceOneGeneration: varOne };
}

// Limiting-dilution transplantation (ELDA-type). Frequency f of initiating cells; design doses with n replicates each.
export function initiatingFrequency({ frequency, doses, n }) {
  // Fisher information for log f under single-hit Poisson: I = sum n x^2 f^2 e^{-xf} / (1 - e^{-xf})
  let I = 0; const rows = [];
  for (const x of doses) { const mu = x * frequency, pTake = 1 - Math.exp(-mu); if (pTake > 0 && pTake < 1) I += n * mu * mu * Math.exp(-mu) / pTake; rows.push({ dose: x, pTake }); }
  const seLog = 1 / Math.sqrt(I);
  return { rows, seLogF: seLog, ci95Low: frequency * Math.exp(-1.96 * seLog), ci95High: frequency * Math.exp(1.96 * seLog), foldWidth: Math.exp(2 * 1.96 * seLog) };
}

/* ----------------------------------------------------- F. apparent stochasticity */

// Hierarchical variance of a condition mean with wells nested in plates nested in days.
export function hierarchicalSe({ sdWell, sdPlate, sdDay, nWells, nPlates, nDays }) {
  const v = sdDay ** 2 / nDays + sdPlate ** 2 / (nDays * nPlates) + sdWell ** 2 / (nDays * nPlates * nWells);
  const naive = (sdWell ** 2 + sdPlate ** 2 + sdDay ** 2) / (nDays * nPlates * nWells);
  return { se: Math.sqrt(v), seNaive: Math.sqrt(naive), understatement: Math.sqrt(v / naive), shareDay: (sdDay ** 2 / nDays) / v, sharePlate: (sdPlate ** 2 / (nDays * nPlates)) / v, shareWell: (sdWell ** 2 / (nDays * nPlates * nWells)) / v };
}

/* ------------------------------------------ G. worked case: spheroid penetration */

// Shell model of antibody penetration with a detection stage.
// Primary: saturated shell depth L = min(R, sqrt(2 D t C0 / BT)) (binding-site barrier regime), occupancy theta = C0/(C0+Kd).
// Secondary per primary ns = nmax (1 - gamma theta) (crowding). Dye density rho = BT theta ns DOL. Brightness per dye q = exp(-rho/rhoQ), so total brightness rho q peaks at rho = rhoQ.
// Attenuation: measured signal at depth z weighted by exp(-(mu + eps rho) z) for whole mounts (eps = 0, mu = 0 for sections).
// Secondary limitation: effective secondary Cs_eff = max(0, Cs - ns * carry * C0); secondary front Ls = sqrt(2 D ts Cs_eff / (BT theta ns)).
export function spheroidHook(P) {
  const { C0, R = 250, D = 10, t = 3600 * 24, BT = 1000, Kd = 1, nmax = 2, gamma = 0.75, DOL = 4, rhoQ = 6000, mu = 0.004, eps = 2e-6, Cs = 20, ts = 3600 * 24, carry = 0, limitSecondary = true, zRef = 100 } = P;
  const theta = C0 / (C0 + Kd);
  const L = Math.min(R, Math.sqrt(2 * D * t * C0 / BT));
  const rhoP = BT * theta;                       // bound primary per volume (nM-equivalent)
  const ns = nmax * (1 - gamma * theta);
  let Ls = Infinity;
  if (limitSecondary) { const csEff = Math.max(0, Cs - ns * carry * C0); Ls = csEff > 0 ? Math.sqrt(2 * D * ts * csEff / Math.max(1e-12, rhoP * ns)) : 0; }
  const Ldet = Math.min(L, Ls);
  const rhoDye = rhoP * ns * DOL;
  const q = Math.exp(-rhoDye / rhoQ);               // phenomenological self-quenching; brightness per dye falls with local dye density
  // integrate over the detected shell: volume-weighted, with attenuation by depth z = R - r
  const n = 200, dr = Ldet / n; let sig = 0, sigNoAtt = 0, direct = 0;
  for (let i = 0; i < n; i++) {
    const z = (i + 0.5) * dr, r = R - z, w = 4 * Math.PI * r * r * dr;
    const att = Math.exp(-(mu + eps * rhoDye) * z);
    sig += rhoDye * q * att * w; sigNoAtt += rhoDye * q * w;
  }
  for (let i = 0; i < n; i++) { const z = (i + 0.5) * (L / n), r = R - z; direct += rhoP * DOL * 4 * Math.PI * r * r * (L / n); }
  const vol = 4 / 3 * Math.PI * R ** 3;
  const rim = rhoDye * q, deep = Ldet >= zRef ? rhoDye * q * Math.exp(-(mu + eps * rhoDye) * zRef) : 0;
  const deepDirect = L >= zRef ? rhoP * DOL : 0;
  return { theta, L, Ls, Ldet, rhoP, ns, rhoDye, q, signal: sig / vol, signalNoAttenuation: sigNoAtt / vol, directLabel: direct / vol, rimIntensity: rim, deepIntensity: deep, deepToRim: rim > 0 ? deep / rim : 0, deepDirect, penetrationFraction: L / R, detectedFraction: Ldet / R };
}
export function hookCurve(P, concentrations) {
  return concentrations.map(C0 => ({ C0, ...spheroidHook({ ...P, C0 }) }));
}

// Wicksell: a sphere of radius R with a stained shell of thickness tTrue, sectioned at random offset z ~ U(0, R).
export function wicksell({ R, t, n = 2000 }) {
  let sumR = 0, sumT = 0, full = 0, sumRatio = 0;
  for (let i = 0; i < n; i++) {
    const z = (i + 0.5) / n * R;
    const Rp = Math.sqrt(R * R - z * z);
    const inner = (R - t) ** 2 - z * z;
    const tp = inner > 0 ? Rp - Math.sqrt(inner) : Rp;
    if (inner <= 0) full++;
    sumR += Rp; sumT += tp; sumRatio += tp / Rp;
  }
  return { meanApparentRadius: sumR / n, meanApparentRadiusTheory: Math.PI * R / 4, meanApparentThickness: sumT / n, biasFactor: (sumT / n) / t, fractionLookFullyStained: full / n, fractionTheory: t / R, meanApparentRelativeDepth: sumRatio / n, trueRelativeDepth: t / R };
}

// Spheroid size from Poisson seeding: cells per drop ~ Poisson(mean), radius ∝ cells^(1/3).
export function seedingSize({ meanCells, cvExtra = 0 }) {
  const cvCount = Math.sqrt(1 / meanCells + cvExtra ** 2);
  return { cvCount, cvRadius: cvCount / 3, p10Radius: (1 - 1.2816 * cvCount / 3), p90Radius: (1 + 1.2816 * cvCount / 3) };
}

export function runAllExamples() {
  return {
    poisson: poissonPartition({ f0: 0.7, phi: 0.05 }), moi: moi({ m: 1, fi: 0.5, phi: 0.05 }), karber: karberDesign({}), lowCopy: lowCopy({ copiesPerReaction: 3 }),
    bottlenecks: bottlenecks({ coverage: 200, k: 3 }), pcr: pcrBranching({ efficiency: 0.9, cycles: 30, templates: 10 }), bursting: bursting({ frequency: 2, size: 20, threshold: 60 }),
    killing: fractionalKilling({ sigma: 1, d50: 1, dose: 3, rho: 0.8 }), ld: (() => { const r = luriaDelbruck({ m: 2 }); return { p0: r.p0, median: r.median, tail20: r.tailAbove(20) }; })(),
    hook: spheroidHook({ C0: 100 }), wicksell: wicksell({ R: 250, t: 50 }),
  };
}
