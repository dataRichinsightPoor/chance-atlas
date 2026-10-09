// Chance Atlas · calculator specifications
// Each calculator declares its parameters and a compute() that returns outputs and a chart from models.js.
import * as M from './models.js';

const pct = (x, dp = 1) => Number.isFinite(x) ? (100 * x).toFixed(dp) + '%' : 'n/a';
const f = (x, dp = 2) => Number.isFinite(x) ? x.toFixed(dp) : 'n/a';
const sci = (x) => Number.isFinite(x) ? (Math.abs(x) >= 1e4 || Math.abs(x) < 1e-3 ? x.toExponential(2) : x.toPrecision(3)) : 'n/a';
const int = (x) => Number.isFinite(x) ? Math.round(x).toLocaleString('en-US') : 'n/a';
const P = (key, label, def, min, max, step, unit = '') => ({ key, label, def, min, max, step, unit });

export const CALCS = {
  poisson: {
    params: [P('f0', 'Empty fraction', 0.7, 0.05, 0.99, 0.01), P('phi', 'Lineages passing the screen', 0.05, 0.001, 1, 0.001)],
    compute: (p) => {
      const r = M.poissonPartition(p);
      const xs = M.linspace(0.01, 0.75, 75);
      const rand = xs.map(g => M.poissonPartition({ f0: 1 - g, phi: 1 }).pMulti), sel = xs.map(g => M.poissonPartition({ f0: 1 - g, phi: p.phi }).pMultiSelected);
      return { outputs: [['Units per partition λ', f(r.lambda, 3)], ['Occupied partition holds ≥2 units', pct(r.pMulti)], ['Selected partition holds ≥2 units', pct(r.pMultiSelected)], ['Risk ratio, selected vs random', f(r.ratio, 2) + ' (→ 2 at low density)']],
        chart: { x: xs.map(g => 100 * g), series: [{ name: 'Random occupied partition', y: rand.map(v => 100 * v) }, { name: 'Partition selected by a screen', y: sel.map(v => 100 * v) }], xLabel: 'Occupied partitions (%)', yLabel: 'Holds two or more units (%)', mark: 100 * (1 - p.f0) } };
    },
  },
  moi: {
    params: [P('m', 'Particles per cell', 1, 0.05, 10, 0.05), P('fi', 'Infectious fraction', 0.5, 0.01, 1, 0.01), P('phi', 'Chance a copy scores', 0.05, 0.001, 1, 0.001)],
    compute: (p) => {
      const r = M.moi(p), xs = M.linspace(0.05, 5, 100);
      return { outputs: [['Functional copies per cell λ', f(r.lambda, 2)], ['Untransduced cells', pct(r.p0)], ['Multi-copy among transduced', pct(r.pMultiAmongTransduced)], ['Multi-copy among cells that score', pct(r.pMultiAmongScored)], ['Mean copies, transduced vs scoring', f(r.meanCopiesTransduced, 2) + ' vs ' + f(r.meanCopiesScored, 2)]],
        chart: { x: xs, series: [{ name: 'Multi-copy among transduced', y: xs.map(l => 100 * M.moi({ m: l, fi: 1, phi: p.phi }).pMultiAmongTransduced) }, { name: 'Multi-copy among scoring cells', y: xs.map(l => 100 * M.moi({ m: l, fi: 1, phi: p.phi }).pMultiAmongScored) }, { name: 'Untransduced', y: xs.map(l => 100 * Math.exp(-l)) }], xLabel: 'Functional copies per cell λ', yLabel: 'Percent of cells', mark: r.lambda } };
    },
  },
  karber: {
    params: [P('n', 'Wells per dilution', 8, 2, 48, 1), P('step', 'Dilution factor', 10, 2, 10, 1), P('levels', 'Dilution levels', 8, 3, 16, 1), P('titerLog10', 'True titer, log10 per well at level 0', 4.3, 1, 8, 0.1)],
    compute: (p) => {
      const r = M.karberDesign({ ...p, logStart: 1 }), ns = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48];
      return { outputs: [['SE of log10 titer', f(r.seLog10, 3)], ['95% interval, fold below and above', f(r.ci95FactorLow, 2) + '× to ' + f(r.ci95FactorHigh, 2) + '×'], ['Smallest reliably resolved ratio (≈ 2.8 SE)', f(10 ** (2.8 * r.seLog10), 1) + '×']],
        chart: { x: ns, series: [{ name: `Dilution factor ${p.step}`, y: ns.map(n => M.karberDesign({ ...p, n, logStart: 1 }).seLog10) }, { name: 'Dilution factor 2 (same range)', y: ns.map(n => M.karberDesign({ ...p, n, step: 2, levels: Math.ceil(p.levels * Math.log2(p.step)), logStart: 1 }).seLog10) }], xLabel: 'Wells per dilution', yLabel: 'SE of log10 titer', mark: p.n } };
    },
  },
  lowcopy: {
    params: [P('copiesPerReaction', 'Mean copies per reaction', 3, 0.1, 30, 0.1), P('efficiency', 'Extraction × amplification efficiency', 1, 0.05, 1, 0.05), P('replicates', 'Replicates', 3, 1, 12, 1)],
    compute: (p) => {
      const r = M.lowCopy(p), xs = M.linspace(0.1, 15, 150);
      return { outputs: [['Effective copies', f(r.effectiveCopies, 2)], ['Detection per replicate', pct(r.pDetect)], ['All replicates positive', pct(r.pAllReplicates)], ['At least one positive', pct(r.pDetectAny)], ['LoD95, nominal copies', f(r.lod95, 2)]],
        chart: { x: xs, series: [{ name: 'One replicate', y: xs.map(c => 100 * (1 - Math.exp(-c * p.efficiency))) }, { name: `All of ${p.replicates} replicates`, y: xs.map(c => 100 * (1 - Math.exp(-c * p.efficiency)) ** p.replicates) }], xLabel: 'Nominal copies per reaction', yLabel: 'Detection probability (%)', mark: p.copiesPerReaction, hline: 95 } };
    },
  },
  counting: {
    params: [P('events', 'Events counted', 50, 5, 5000, 5), P('targetCv', 'Target CV', 0.05, 0.01, 0.5, 0.01)],
    compute: (p) => {
      const r = M.countingCv(p), xs = M.logspace(5, 1e5, 60);
      return { outputs: [['CV from counting alone', pct(r.cv)], ['95% interval on the count', int(r.ci95Low) + ' to ' + int(r.ci95High)], ['Events needed for target CV', int(r.eventsForTarget)]],
        chart: { x: xs, series: [{ name: 'Poisson CV', y: xs.map(n => 100 / Math.sqrt(n)) }], xLabel: 'Events counted', yLabel: 'CV (%)', logX: true, logY: true, mark: p.events } };
    },
  },
  bottleneck: {
    params: [P('coverage', 'Cells per element at each bottleneck', 100, 3, 1000, 1), P('k', 'Bottlenecks', 3, 1, 8, 1), P('foldThreshold', 'Fold-change call threshold', 2, 1.2, 8, 0.1)],
    compute: (p) => {
      const r = M.bottlenecks(p), xs = M.logspace(3, 1000, 40);
      return { outputs: [['Neutral CV after all bottlenecks', pct(r.cvAfterK)], ['Theory √(k / c)', pct(r.cvOneTheory)], [`Neutral elements called ≥${p.foldThreshold}× depleted`, pct(r.falseDepleted, 2)], [`Called ≥${p.foldThreshold}× enriched`, pct(r.falseEnriched, 2)], ['Elements lost entirely', pct(r.lostAfterK, 3)]],
        chart: { x: xs, series: [{ name: 'Neutral CV', y: xs.map(c => 100 * M.bottlenecks({ coverage: c, k: p.k, foldThreshold: p.foldThreshold }).cvAfterK) }, { name: 'False depletion calls', y: xs.map(c => 100 * M.bottlenecks({ coverage: c, k: p.k, foldThreshold: p.foldThreshold }).falseDepleted) }], xLabel: 'Coverage, cells per element', yLabel: 'Percent', logX: true, logY: true, mark: p.coverage } };
    },
  },
  undersampling: {
    params: [P('depthPerMolecule', 'Reads per molecule', 1, 0.05, 10, 0.05), P('clonotypeFrequency', 'Clonotype frequency', 1e-4, 1e-6, 1e-2, 1e-6), P('reads', 'Reads sequenced', 1e5, 1e3, 1e7, 1e3)],
    compute: (p) => {
      const r = M.undersampling(p), xs = M.logspace(0.05, 20, 60);
      return { outputs: [['Molecules seen at least once', pct(r.fractionSeen)], ['Duplicate reads', pct(r.duplicateRate)], ['This clonotype seen', pct(r.pClonotypeSeen)], ['Reads for 95% chance of seeing it', int(r.readsFor95)]],
        chart: { x: xs, series: [{ name: 'Molecules seen', y: xs.map(m => 100 * (1 - Math.exp(-m))) }, { name: 'Duplicate reads', y: xs.map(m => 100 * (1 - (1 - Math.exp(-m)) / m)) }], xLabel: 'Reads per molecule', yLabel: 'Percent', logX: true, mark: p.depthPerMolecule } };
    },
  },
  pcr: {
    params: [P('efficiency', 'Per-cycle efficiency', 0.9, 0.3, 1, 0.01), P('cycles', 'Cycles', 30, 5, 40, 1), P('templates', 'Starting templates', 10, 1, 10000, 1), P('umis', 'Possible UMIs', 1048576, 1024, 1.7e7, 1024), P('molecules', 'Molecules labeled', 1e5, 1e3, 1e7, 1e3)],
    compute: (p) => {
      const r = M.pcrBranching(p), xs = M.linspace(0.3, 1, 71);
      return { outputs: [['Mean amplification', sci(r.meanAmplification)], ['CV of yield from one template', pct(r.cvSingleTemplate)], [`CV from ${int(p.templates)} templates`, pct(r.cv, 2)], ['UMI collisions', pct(r.umiCollisionFraction, 2)]],
        chart: { x: xs, series: [{ name: '1 template', y: xs.map(e => 100 * M.pcrBranching({ efficiency: e, cycles: p.cycles, templates: 1 }).cv) }, { name: '10 templates', y: xs.map(e => 100 * M.pcrBranching({ efficiency: e, cycles: p.cycles, templates: 10 }).cv) }, { name: '100 templates', y: xs.map(e => 100 * M.pcrBranching({ efficiency: e, cycles: p.cycles, templates: 100 }).cv) }], xLabel: 'Per-cycle efficiency', yLabel: 'CV of product (%)', mark: p.efficiency } };
    },
  },
  ado: {
    params: [P('q', 'Per-allele capture probability', 0.8, 0.05, 1, 0.01)],
    compute: (p) => {
      const r = M.alleleDropout(p), xs = M.linspace(0.05, 1, 96);
      return { outputs: [['Both alleles', pct(r.pBoth)], ['One allele (dropout)', pct(r.pOne)], ['No amplification', pct(r.pNone)], ['Dropout among amplified', pct(r.adoAmongAmplified)]],
        chart: { x: xs, series: [{ name: 'Dropout among amplified', y: xs.map(q => 100 * M.alleleDropout({ q }).adoAmongAmplified) }, { name: 'No amplification', y: xs.map(q => 100 * (1 - q) ** 2) }], xLabel: 'Per-allele capture probability', yLabel: 'Percent of cells', mark: p.q } };
    },
  },
  capture: {
    params: [P('trueCount', 'True molecules per cell', 5, 1, 200, 1), P('capture', 'Capture probability', 0.1, 0.01, 1, 0.01), P('trueCv', 'True biological CV', 0.3, 0, 2, 0.05)],
    compute: (p) => {
      const r = M.captureDropout(p), xs = M.logspace(1, 500, 60);
      return { outputs: [['Cells with zero observed', pct(r.pZero)], ['Mean observed count', f(r.meanObserved, 2)], ['Observed CV vs true CV', pct(r.cvObserved) + ' vs ' + pct(r.cvTrue)], ['Share of observed variance that is technical', pct(r.noiseShare)]],
        chart: { x: xs, series: [{ name: 'Zero observed', y: xs.map(n => 100 * (1 - p.capture) ** n) }, { name: 'Technical share of variance', y: xs.map(n => 100 * M.captureDropout({ trueCount: n, capture: p.capture, trueCv: p.trueCv }).noiseShare) }], xLabel: 'True molecules per cell', yLabel: 'Percent', logX: true, mark: p.trueCount } };
    },
  },
  shot: {
    params: [P('photons', 'Photons per measurement', 100, 1, 10000, 1), P('readNoise', 'Read noise (electrons)', 1.5, 0, 20, 0.1), P('threshold', 'Threshold (photons)', 50, 1, 10000, 1)],
    compute: (p) => {
      const r = M.shotNoise(p), xs = M.logspace(1, 1e4, 60);
      return { outputs: [['Signal-to-noise', f(r.snr, 1)], ['CV', pct(r.cv)], ['Falls below threshold', pct(r.pBelowThreshold, 2)]],
        chart: { x: xs, series: [{ name: 'SNR, shot noise only', y: xs.map(n => Math.sqrt(n)) }, { name: `SNR with read noise ${p.readNoise}`, y: xs.map(n => n / Math.sqrt(n + p.readNoise ** 2)) }], xLabel: 'Photons', yLabel: 'SNR', logX: true, logY: true, mark: p.photons } };
    },
  },
  bleach: {
    params: [P('subunits', 'True subunits', 4, 1, 8, 1), P('pLabel', 'Effective labeling probability', 0.75, 0.1, 1, 0.01)],
    compute: (p) => {
      const r = M.bleachingSteps(p);
      return { outputs: [['Modal step count', r.modeSteps], ['Complexes showing all steps', pct(r.pFullCount)], ['Mean steps', f(r.meanSteps, 2)], ['Invisible complexes (zero steps)', pct(r.pZero, 2)]],
        chart: { x: r.dist.map((_, k) => k), series: [{ name: 'Observed step count', y: r.dist.map(v => 100 * v) }], xLabel: 'Bleaching steps observed', yLabel: 'Percent of complexes', bars: true } };
    },
  },
  fcs: {
    params: [P('concentration_nM', 'Concentration (nM)', 1, 0.01, 1000, 0.01), P('volume_fL', 'Detection volume (fL)', 1, 0.1, 1000, 0.1)],
    compute: (p) => {
      const r = M.numberFluctuations(p), xs = M.logspace(0.01, 1000, 60);
      return { outputs: [['Mean molecules in volume', f(r.meanMolecules, 2)], ['Relative variance 1/⟨N⟩', f(r.relativeVariance, 3)], ['CV of instantaneous signal', pct(r.cv)]],
        chart: { x: xs, series: [{ name: `CV at ${p.volume_fL} fL`, y: xs.map(c => 100 * M.numberFluctuations({ concentration_nM: c, volume_fL: p.volume_fL }).cv) }], xLabel: 'Concentration (nM)', yLabel: 'CV (%)', logX: true, logY: true, mark: p.concentration_nM } };
    },
  },
  bp: {
    params: [P('radius_um', 'Sensor radius (µm)', 5, 0.01, 20, 0.01), P('concentration_nM', 'Ligand (nM)', 0.01, 1e-4, 100, 1e-4), P('D', 'Diffusion coefficient (µm²/s)', 100, 1, 500, 1), P('T', 'Integration time (s)', 60, 0.1, 3600, 0.1)],
    compute: (p) => {
      const r = M.bergPurcell(p), xs = M.logspace(1e-4, 100, 60);
      return { outputs: [['Error floor, perfectly monitoring sphere', pct(r.relativeError)], ['Error floor, perfect absorber', pct(r.relAbsorber)], ['Molecular encounters in T', sci(r.encounters)], ['Time for 5% precision (monitor)', sci(r.timeFor5pct) + ' s']],
        chart: { x: xs, series: [{ name: 'Perfectly monitoring sphere √(3/(5πDacT))', y: xs.map(c => 100 * M.bergPurcell({ ...p, concentration_nM: c }).relativeError) }, { name: 'Perfect absorber 1/√(4πDacT)', y: xs.map(c => 100 * M.bergPurcell({ ...p, concentration_nM: c }).relAbsorber) }], xLabel: 'Ligand (nM)', yLabel: 'Relative error (%)', logX: true, logY: true, mark: p.concentration_nM } };
    },
  },
  channel: {
    params: [P('N', 'Channels in patch', 100, 1, 2000, 1), P('i', 'Unitary current (pA)', 1, 0.05, 20, 0.05), P('pOpen', 'Open probability', 0.5, 0.01, 0.99, 0.01)],
    compute: (p) => {
      const r = M.channelNoise(p);
      return { outputs: [['Mean current (pA)', f(r.meanCurrent, 1)], ['Variance (pA²)', f(r.variance, 2)], ['CV of current', pct(r.cv)]],
        chart: { x: r.parabola.map(x => x.mean), series: [{ name: 'Variance vs mean', y: r.parabola.map(x => x.variance) }], xLabel: 'Mean current (pA)', yLabel: 'Variance (pA²)', mark: r.meanCurrent } };
    },
  },
  burst: {
    params: [P('frequency', 'Bursts per mRNA lifetime', 2, 0.1, 20, 0.1), P('size', 'Burst size (mRNA)', 20, 1, 200, 1), P('threshold', 'Positive threshold (mRNA)', 60, 1, 500, 1)],
    compute: (p) => {
      const r = M.bursting(p), xs = M.linspace(1, Math.max(10, 4 * r.mean), 60);
      return { outputs: [['Mean mRNA', f(r.mean, 1)], ['Fano factor 1 + b', f(r.fano, 1)], ['Percent positive', pct(r.pctPositive)], ['After doubling frequency', pct(r.pctPositiveDoubleFreq)], ['After doubling size', pct(r.pctPositiveDoubleSize)]],
        chart: { x: xs, series: [{ name: 'Baseline', y: xs.map(T => 100 * M.bursting({ ...p, threshold: T }).pctPositive) }, { name: 'Double frequency', y: xs.map(T => 100 * M.bursting({ ...p, frequency: 2 * p.frequency, threshold: T }).pctPositive) }, { name: 'Double size', y: xs.map(T => 100 * M.bursting({ ...p, size: 2 * p.size, threshold: T }).pctPositive) }], xLabel: 'Threshold (mRNA per cell)', yLabel: 'Percent positive', mark: p.threshold } };
    },
  },
  sort: {
    params: [P('sigma', 'SD of log protein level', 0.5, 0.05, 2, 0.05), P('topFraction', 'Fraction sorted', 0.05, 0.005, 0.5, 0.005), P('tau', 'Mixing time (generations)', 2, 0.2, 20, 0.1), P('t', 'Time since sort (generations)', 4, 0, 40, 0.5)],
    compute: (p) => {
      const r = M.sortRelaxation(p), xs = M.linspace(0, Math.max(10, 5 * p.tau), 80);
      return { outputs: [['Enrichment at sort (SD units)', f(r.sortedMeanSd, 2)], ['Enrichment remaining', pct(r.enrichmentRemaining)], ['Mean now (SD units)', f(r.meanAtT_sd, 2)], ['Fold above parent, linear scale', f(r.foldMeanLinear, 2) + '×']],
        chart: { x: xs, series: [{ name: 'Mean of sorted population (SD units)', y: xs.map(t => M.sortRelaxation({ ...p, t }).meanAtT_sd) }, { name: 'SD of sorted population', y: xs.map(t => M.sortRelaxation({ ...p, t }).sdAtT) }], xLabel: 'Generations since sort', yLabel: 'Standard-deviation units', mark: p.t } };
    },
  },
  kill: {
    params: [P('sigma', 'SD of log protective state', 1, 0.1, 3, 0.05), P('d50', 'Dose killing half the cells', 1, 0.1, 100, 0.1), P('dose', 'Dose', 3, 0.1, 100, 0.1), P('rho', 'State correlation between doses', 0.8, 0, 0.99, 0.01)],
    compute: (p) => {
      const r = M.fractionalKilling(p), xs = M.logspace(p.d50 / 30, p.d50 * 30, 60);
      return { outputs: [['Survival, first dose', pct(r.survival1)], ['Survival of survivors, second dose', pct(r.survival2Given)], ['Two doses, actual vs independent', pct(r.survivalTwoDoses, 2) + ' vs ' + pct(r.naiveTwoDoses, 2)], ['Second dose kills this much less', f(r.foldLessKilling, 2) + '× smaller fraction']],
        chart: { x: xs, series: [{ name: 'Survival, first dose', y: xs.map(d => 100 * M.fractionalKilling({ ...p, dose: d }).survival1) }, { name: 'Survival of survivors, second dose', y: xs.map(d => 100 * M.fractionalKilling({ ...p, dose: d }).survival2Given) }], xLabel: 'Dose', yLabel: 'Survival (%)', logX: true, mark: p.dose } };
    },
  },
  persist: {
    params: [P('f', 'Persister fraction', 0.001, 1e-5, 0.1, 1e-5), P('ks', 'Kill rate, main population (per h)', 2, 0.1, 10, 0.1), P('kp', 'Kill rate, persisters (per h)', 0.01, 0, 1, 0.005), P('t', 'Exposure (h)', 6, 0.5, 48, 0.5)],
    compute: (p) => {
      const r = M.biphasicKill(p), xs = M.linspace(0, 24, 97);
      return { outputs: [['Surviving fraction', sci(r.survivingFraction)], ['Log10 kill', f(r.log10Kill, 2)], ['Time to 99.9% kill (h)', f(r.timeTo999, 1)], ['Same without persisters (h)', f(r.timeTo999NoPersisters, 2)], ['Persisters among survivors', pct(r.persisterShareOfSurvivors)]],
        chart: { x: xs, series: [{ name: 'With persisters', y: xs.map(t => Math.log10(M.biphasicKill({ ...p, t }).survivingFraction)) }, { name: 'Without persisters', y: xs.map(t => -p.ks * t / Math.LN10) }], xLabel: 'Exposure (h)', yLabel: 'log10 surviving fraction', mark: p.t, yMin: -8 } };
    },
  },
  switch: {
    params: [P('kOff', 'Silencing rate per generation', 0.05, 0, 1, 0.005), P('kOn', 'Reactivation rate per generation', 0.01, 0, 1, 0.005), P('generations', 'Generations', 20, 0, 200, 1)],
    compute: (p) => {
      const r = M.switching(p), xs = M.linspace(0, Math.max(40, 4 * r.relaxationGenerations), 80);
      return { outputs: [['Fraction expressing now', pct(r.fractionOn)], ['Equilibrium fraction', pct(r.equilibriumOn)], ['Relaxation time (generations)', f(r.relaxationGenerations, 1)]],
        chart: { x: xs, series: [{ name: 'Fraction expressing', y: xs.map(g => 100 * M.switching({ ...p, generations: g }).fractionOn) }], xLabel: 'Generations', yLabel: 'Percent expressing', mark: p.generations } };
    },
  },
  clone: {
    params: [P('meanDivisions', 'Mean divisions', 10, 1, 30, 1), P('cvDivisionTime', 'CV of division time', 0.2, 0.02, 0.6, 0.01)],
    compute: (p) => {
      const r = M.cloneSize(p), xs = M.linspace(0.02, 0.6, 59);
      return { outputs: [['Median clone size', int(r.medianSize)], ['CV of clone size', pct(r.cvSize)], ['10th to 90th percentile', int(r.p10) + ' to ' + int(r.p90)], ['Fold range', f(r.foldP90P10, 1) + '×']],
        chart: { x: xs, series: [{ name: 'Fold range, 90th/10th percentile', y: xs.map(c => M.cloneSize({ ...p, cvDivisionTime: c }).foldP90P10) }], xLabel: 'CV of division time', yLabel: 'Fold range of clone size', mark: p.cvDivisionTime } };
    },
  },
  ld: {
    params: [P('m', 'Expected mutations per culture', 2, 0.1, 20, 0.1)],
    compute: (p) => {
      const r = M.luriaDelbruck({ m: p.m, rmax: 400 }), K = Math.min(60, Math.max(15, Math.round(12 * p.m)));
      const xs = Array.from({ length: K + 1 }, (_, i) => i);
      const pois = xs.map(k => 100 * M.poissonPmf(k, r.meanTruncated));
      return { outputs: [['P(no mutants) = e^−m', pct(r.p0, 2)], ['Median mutants', r.median], ['Cultures with more than 10× the median', pct(r.tailAbove(10 * Math.max(1, r.median)), 2)], ['Estimate m from zero class', 'm = −ln(fraction of cultures with no mutants)']],
        chart: { x: xs, series: [{ name: 'Luria–Delbrück', y: r.pmf.slice(0, K + 1).map(v => 100 * v) }, { name: 'Poisson with the same mean', y: pois }], xLabel: 'Mutants per culture', yLabel: 'Percent of cultures', logY: true, yMin: 0.01 } };
    },
  },
  drift: {
    params: [P('p0', 'Initial subclone frequency', 0.2, 0.001, 0.99, 0.001), P('N', 'Effective cells per passage', 500, 10, 100000, 10), P('s', 'Selection per passage', 0.05, -0.5, 0.5, 0.01), P('passages', 'Passages', 20, 1, 100, 1)],
    compute: (p) => {
      const r = M.clonalDrift(p), xs = r.trajectory.map((_, i) => i);
      return { outputs: [['Final frequency', pct(r.final)], ['SD from drift', pct(r.finalSd)], ['Drift-only SD, no selection', pct(Math.sqrt(r.driftOnlyVar))]],
        chart: { x: xs, series: [{ name: 'Expected frequency', y: r.trajectory.map(v => 100 * v) }, { name: '+1 SD', y: r.trajectory.map((v, i) => 100 * Math.min(1, v + r.sd[i])) }, { name: '−1 SD', y: r.trajectory.map((v, i) => 100 * Math.max(0, v - r.sd[i])) }], xLabel: 'Passage', yLabel: 'Subclone frequency (%)' } };
    },
  },
  mt: {
    params: [P('h', 'Parental heteroplasmy', 0.3, 0.01, 0.99, 0.01), P('n', 'Segregating units', 10, 1, 1000, 1), P('generationsOfDrift', 'Generations of drift', 1, 1, 50, 1)],
    compute: (p) => {
      const r = M.mtBottleneck(p), xs = M.linspace(0, 1, 101);
      return { outputs: [['SD of heteroplasmy in progeny', pct(r.sd)], ['Progeny above 60% mutant', pct(r.pAbove(0.6))], ['Progeny below 10% mutant', pct(r.pBelow(0.1))]],
        chart: { x: xs.map(v => 100 * v), series: [{ name: 'Heteroplasmy distribution (normal approximation)', y: xs.map(v => M.normPdf((v - p.h) / r.sd) / r.sd / 100) }], xLabel: 'Heteroplasmy in progeny (%)', yLabel: 'Density', mark: 100 * p.h } };
    },
  },
  elda: {
    params: [P('frequency', 'True frequency (1 in …)', 1000, 10, 1e6, 10), P('n', 'Replicates per dose', 6, 1, 48, 1), P('lowDose', 'Lowest dose (cells)', 100, 1, 1e6, 1), P('doses', 'Doses, threefold apart', 4, 2, 8, 1)],
    compute: (p) => {
      const f0 = 1 / p.frequency, doses = Array.from({ length: p.doses }, (_, i) => p.lowDose * 3 ** i);
      const r = M.initiatingFrequency({ frequency: f0, doses, n: p.n }), xs = M.logspace(doses[0] / 10, doses[doses.length - 1] * 10, 60);
      return { outputs: [['Doses', doses.map(int).join(', ')], ['Expected take rates', r.rows.map(x => pct(x.pTake, 0)).join(', ')], ['SE of ln f', f(r.seLogF, 2)], ['95% interval', '1 in ' + int(1 / r.ci95High) + ' to 1 in ' + int(1 / r.ci95Low)], ['Fold width of interval', f(r.foldWidth, 1) + '×']],
        chart: { x: xs, series: [{ name: 'Expected take rate', y: xs.map(x => 100 * (1 - Math.exp(-x * f0))) }, { name: 'Information per animal (scaled)', y: xs.map(x => { const mu = x * f0; return 100 * mu * mu * Math.exp(-mu) / (1 - Math.exp(-mu)) / 0.65; }) }], xLabel: 'Cells transplanted', yLabel: 'Percent', logX: true } };
    },
  },
  hier: {
    params: [P('sdWell', 'SD between wells', 0.1, 0, 1, 0.01), P('sdPlate', 'SD between plates', 0.1, 0, 1, 0.01), P('sdDay', 'SD between days', 0.15, 0, 1, 0.01), P('nWells', 'Wells per plate', 8, 1, 96, 1), P('nPlates', 'Plates per day', 2, 1, 20, 1), P('nDays', 'Days', 3, 1, 20, 1)],
    compute: (p) => {
      const r = M.hierarchicalSe(p), xs = Array.from({ length: 12 }, (_, i) => i + 1);
      return { outputs: [['True SE of the condition mean', f(r.se, 3)], ['Naive SE (all wells independent)', f(r.seNaive, 3)], ['Understatement factor', f(r.understatement, 2) + '×'], ['Variance share: day / plate / well', pct(r.shareDay, 0) + ' / ' + pct(r.sharePlate, 0) + ' / ' + pct(r.shareWell, 0)]],
        chart: { x: xs, series: [{ name: 'Add days', y: xs.map(n => M.hierarchicalSe({ ...p, nDays: n }).se) }, { name: 'Add plates per day', y: xs.map(n => M.hierarchicalSe({ ...p, nPlates: n }).se) }, { name: 'Add wells per plate (×8)', y: xs.map(n => M.hierarchicalSe({ ...p, nWells: 8 * n }).se) }], xLabel: 'Replicates added at one level', yLabel: 'SE of condition mean' } };
    },
  },
  hook: {
    params: [P('R', 'Spheroid radius (µm)', 250, 50, 600, 10), P('t', 'Primary incubation (h)', 24, 1, 96, 1), P('BT', 'Antigen density (nM)', 1000, 10, 10000, 10), P('Kd', 'Primary Kd (nM)', 1, 0.01, 100, 0.01), P('D', 'Diffusion coefficient (µm²/s)', 10, 1, 50, 0.5),
      P('nmax', 'Secondaries per primary, uncrowded', 2, 1, 3, 0.1), P('gamma', 'Crowding loss at saturation', 0.75, 0, 1, 0.05), P('DOL', 'Dyes per secondary', 4, 1, 10, 0.5), P('rhoQ', 'Quenching density (nM dye)', 6000, 500, 50000, 100),
      P('mu', 'Tissue attenuation (per µm)', 0.004, 0, 0.02, 0.0005), P('eps', 'Dye attenuation (per µm per nM)', 2e-6, 0, 2e-5, 1e-7), P('Cs', 'Secondary concentration (nM)', 20, 1, 500, 1), P('carry', 'Free primary carried into secondary step', 0.01, 0, 0.1, 0.001), P('zRef', 'Reference depth (µm)', 100, 10, 500, 10)],
    compute: (p) => {
      const Q = { ...p, t: p.t * 3600, ts: p.t * 3600, limitSecondary: true }, C = M.logspace(0.03, 10000, 60), here = M.spheroidHook({ ...Q, C0: 10 });
      const all = M.hookCurve(Q, C), ideal = M.hookCurve({ ...Q, gamma: 0, rhoQ: 1e15, eps: 0, mu: 0, limitSecondary: false }, C);
      const norm = (ys) => { const m = Math.max(...ys); return ys.map(y => 100 * y / m); };
      const limitOnly = M.hookCurve({ ...Q, gamma: 0, rhoQ: 1e15, eps: 0, mu: 0, carry: 0 }, C).map(x => x.signal);
      const sectionView = M.hookCurve({ ...Q, eps: 0, mu: 0 }, C).map(x => x.signal);
      const peak = all.reduce((m, x) => x.signal > m.signal ? x : m);
      const last = all[C.length - 1];
      const mech = [['carryover', M.hookCurve({ ...Q, gamma: 0, rhoQ: 1e15, eps: 0, mu: 0 }, C)], ['crowding', M.hookCurve({ ...Q, carry: 0, rhoQ: 1e15, eps: 0, mu: 0, limitSecondary: false }, C)], ['self-quenching', M.hookCurve({ ...Q, carry: 0, gamma: 0, eps: 0, mu: 0, limitSecondary: false }, C)], ['attenuation', M.hookCurve({ ...Q, carry: 0, gamma: 0, rhoQ: 1e15, limitSecondary: false }, C)]]
        .map(([n, cur]) => { const ys = cur.map(x => x.signal), mx = Math.max(...ys); return [n, 1 - ys[ys.length - 1] / mx]; }).sort((a, b) => b[1] - a[1]);
      return { outputs: [['At 10 nM: primary front / detected front', f(here.L, 0) + ' / ' + f(here.Ldet, 0) + ' µm'], ['At 10 nM: secondaries per primary', f(here.ns, 2)], ['Signal peaks at', sci(peak.C0) + ' nM'], ['Signal at 10 µM relative to peak', pct(last.signal / peak.signal, 0)], ['Largest single-mechanism fall at 10 µM', mech[0][0] + ' (' + pct(mech[0][1], 0) + ')'], ['Rim vs ' + p.zRef + ' µm intensity at 10 µM', f(last.rimIntensity, 0) + ' vs ' + f(last.deepIntensity, 0)]],
        chart: { x: C, series: [{ name: 'Your settings, whole mount', y: norm(all.map(x => x.signal)) }, { name: 'Your settings, section (no attenuation)', y: norm(sectionView) }, { name: 'Secondary limitation only', y: norm(limitOnly) }, { name: 'Directly labeled primary', y: norm(all.map(x => x.directLabel)) }], xLabel: 'Primary concentration (nM)', yLabel: 'Integrated signal (% of each curve’s maximum)', logX: true } };
    },
  },
  wicksell: {
    params: [P('R', 'Spheroid radius (µm)', 250, 50, 600, 10), P('t', 'True stained shell (µm)', 50, 5, 300, 5)],
    compute: (p) => {
      const t = Math.min(p.t, p.R), r = M.wicksell({ R: p.R, t }), xs = M.linspace(0.02, 1, 50);
      return { outputs: [['Mean apparent section radius', f(r.meanApparentRadius, 0) + ' µm (πR/4 = ' + f(r.meanApparentRadiusTheory, 0) + ')'], ['Mean apparent shell', f(r.meanApparentThickness, 0) + ' µm, true ' + f(t, 0)], ['Overstatement', pct(r.biasFactor - 1, 0)], ['Sections that look fully stained', pct(r.fractionLookFullyStained, 0)]],
        chart: { x: xs.map(v => 100 * v), series: [{ name: 'Mean apparent shell / true shell', y: xs.map(v => M.wicksell({ R: p.R, t: v * p.R, n: 400 }).biasFactor) }, { name: 'Fraction looking fully stained', y: xs.map(v => v) }], xLabel: 'True shell as % of radius', yLabel: 'Ratio or fraction', mark: 100 * t / p.R } };
    },
  },
  seed: {
    params: [P('meanCells', 'Mean cells per drop', 400, 10, 10000, 10), P('cvExtra', 'Extra dispensing CV', 0.05, 0, 0.5, 0.01)],
    compute: (p) => {
      const r = M.seedingSize(p), xs = M.logspace(10, 10000, 60);
      return { outputs: [['CV of cell count', pct(r.cvCount)], ['CV of radius', pct(r.cvRadius)], ['10th to 90th percentile radius', pct(r.p10Radius, 0) + ' to ' + pct(r.p90Radius, 0) + ' of mean']],
        chart: { x: xs, series: [{ name: 'Radius CV, Poisson only', y: xs.map(n => 100 / Math.sqrt(n) / 3) }, { name: `Radius CV with ${pct(p.cvExtra, 0)} dispensing CV`, y: xs.map(n => 100 * M.seedingSize({ meanCells: n, cvExtra: p.cvExtra }).cvRadius) }], xLabel: 'Mean cells per drop', yLabel: 'CV of radius (%)', logX: true, mark: p.meanCells } };
    },
  },
};
