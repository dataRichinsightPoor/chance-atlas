# Stochastic Processes That Shape Bioassay Outcomes

A working catalog for 99 Small Problems and Data-Rich, Insight-Poor. Each entry lists the process, its usual statistical form, the assays it affects, how it changes the reported result, and the diagnostic that exposes it. Entries marked as candidates share the structure that made limiting-dilution cloning a good installment: a standard calculation that conditions on the wrong event, a closed form that corrects it, and a parameter that a laboratory can measure.

The catalog separates five places where randomness enters an assay: sampling discrete units, amplifying molecules, detecting signal, the biology of individual cells, and the history of the population. A sixth section covers variation that looks stochastic but is not.

## 1. Sampling and partitioning of discrete units

These processes follow from the fact that molecules, cells and virions are counted in whole numbers. At low copy number, Poisson or binomial sampling sets a limit that no instrument improvement removes.

### Poisson loading of partitions

Statistical form: units per partition ~ Poisson(λ), with λ estimated from the negative fraction as \( \lambda = -\ln f_0 \).

Assays: limiting-dilution cloning, droplet and chip digital PCR ([Hindson et al., Anal Chem 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3216358/)), single-molecule digital ELISA ([Rissin et al., Nat Biotechnol 2010](https://pubmed.ncbi.nlm.nih.gov/20495550/)), single-cell encapsulation, microwell and bead-loading assays.

Effect on outcome: precision is best at an intermediate occupancy and collapses near saturation. Multi-occupancy hides inside positive partitions. Selecting partitions by signal size-biases toward multiply occupied ones, the mechanism behind No. 07.

Hidden failure: non-Poisson loading from aggregates, partition volume variation and target clustering on long DNA all inflate the true dispersion. Overdispersion in partition counts is the diagnostic.

Candidate: covered by No. 07 for cloning. A digital-assay sequel (droplet volume variation and linked targets) is feasible.

### Multiplicity of infection and transduction

Statistical form: infectious particles per cell ~ Poisson(MOI × fraction infectious).

Assays: viral infection assays, lentiviral transduction, pooled CRISPR and shRNA screens, reporter assays after transfection.

Effect on outcome: at MOI 1, roughly a third of cells receive nothing and a quarter receive two or more, so a "transduced population" is a mixture of copy numbers. In pooled screens, multiple integrations per cell couple phenotypes of different guides. The infectious fraction is rarely known, so the MOI by particle count and by function differ.

Candidate: strong. The conditional "cells that score in a screen carry more integrations" is the same size bias as cloning, applied to guide calls.

### Endpoint-dilution titration

Statistical form: wells infected at dilution i ~ Binomial(n, 1 − e^(−dose_i)).

Assays: TCID50, LD50 and ID50 endpoints, plaque assays, most-probable-number counts. Reed and Muench published the interpolation method in 1938 ([Am J Epidemiol](https://academic.oup.com/aje/article-abstract/27/3/493/99616)); Spearman–Kärber is the common alternative ([Ramakrishnan, World J Virol 2016](https://www.wjgnet.com/2220-3249/abstract/v5/i2/85.htm)). Plaque purification rests on the same single-particle logic ([Dulbecco and Vogt, J Exp Med 1954](https://pmc.ncbi.nlm.nih.gov/articles/PMC2180341/)).

Effect on outcome: interpolation methods carry no likelihood, so their intervals are often missing or wrong, and a fixed number of wells per dilution limits resolution to about a factor of the dilution step.

Candidate: moderate. Maximum-likelihood titration with intervals is a well-defined tool, but the field is well trodden.

### Low-copy detection limits

Statistical form: target copies in a reaction ~ Poisson(c × V).

Assays: qPCR and RT-qPCR near the limit of detection, rare-variant detection, minimal residual disease, pathogen detection, low-input sequencing. LoD definitions for qPCR handle this sampling explicitly ([Forootan et al., Biomol Detect Quantif 2017](https://pubmed.ncbi.nlm.nih.gov/28702366/)).

Effect on outcome: at a mean of three copies per reaction, about 5 percent of reactions contain none, which sets an LoD floor near three copies regardless of chemistry. Replicate Ct values at low input widen for sampling reasons before amplification adds its own noise.

Candidate: strong for a short installment. "Your LoD is set by Poisson before it is set by the assay" has an inspectable closed form.

### Cell seeding and counting

Statistical form: cells per well ~ Poisson or overdispersed (negative binomial) with clumping and settling.

Assays: viability and proliferation plates, colony-forming assays, cytotoxicity, flow event counts.

Effect on outcome: at low seeding the coefficient of variation of the starting count is \( 1/\sqrt{N} \), which propagates directly into fold-growth readouts. Flow-gated rare populations follow Poisson counting, so a gate with 50 events carries a CV of about 14 percent before any biology.

Candidate: moderate, as a rare-event gating installment.

### Library representation and bottlenecks

Statistical form: multinomial sampling of library elements at each step (transduction, selection, passaging, PCR, sequencing).

Assays: pooled CRISPR, shRNA and ORF screens; phage, yeast and mRNA display; DNA-encoded libraries; barcoded lineage tracing. Coverage requirements are central to screen design ([Hart et al., Cell 2015](http://sites.utoronto.ca/sidhulab/pdf/Hart_2015_free.pdf)).

Effect on outcome: each bottleneck adds variance that looks like a phenotype. Guides lost by sampling read as depleted. Selection rounds amplify early random winners (jackpotting), so enrichment ranks partly reflect who got lucky in round one.

Candidate: very strong. Neutral drift through serial bottlenecks has a closed-form variance, and the false-depletion rate at a given coverage is easy to show.

### Sequencing and mass spectrometry undersampling

Statistical form: reads ~ multinomial; data-dependent MS2 selection behaves as stochastic sampling of precursors.

Assays: RNA-seq, amplicon and barcode counting, shotgun proteomics. In one study, more than 100,000 peptide features eluted per run, but only about 16 percent were targeted for MS/MS ([Michalski et al., J Proteome Res 2011](https://pubmed.ncbi.nlm.nih.gov/21309581/)).

Effect on outcome: missing values are not missing at random. They concentrate in low-abundance species, which biases fold changes when they are imputed or dropped.

Candidate: moderate. It overlaps with the censoring logic of Censor Check.

## 2. Molecular amplification

### Stochastic early-cycle PCR

Statistical form: branching process. Each template copies with probability equal to the per-cycle efficiency, so early-cycle luck is multiplied exponentially.

Assays: qPCR at low input, single-cell PCR, amplicon sequencing, library preparation. Kebschull and Zador separated bias, stochasticity, template switching and polymerase error as distinct distortions ([Nucleic Acids Res 2015](https://pubmed.ncbi.nlm.nih.gov/26187991/)).

Effect on outcome: relative abundances drift away from their input proportions, and the drift is largest for the rarest templates.

Mitigation: unique molecular identifiers count molecules rather than reads ([Islam et al., Nat Methods 2014](https://pubmed.ncbi.nlm.nih.gov/24363023/)).

Candidate: strong. The Galton–Watson variance of a PCR product has a clean closed form, and the cost of UMI collisions at high depth is a natural second panel.

### Allele dropout and preferential amplification

Statistical form: Bernoulli failure of one template among few.

Assays: single-cell genotyping, preimplantation genetic testing, low-input forensic and clinical genotyping. Walsh and colleagues described the mechanisms in 1992 ([PCR Methods Appl](https://pubmed.ncbi.nlm.nih.gov/1477658/)).

Effect on outcome: heterozygotes are called homozygous.

### Capture and dropout in single-cell sequencing

Statistical form: binomial or Poisson capture of each transcript, compounded with amplification.

Assays: scRNA-seq and other single-cell omics. Technical noise in these data is large and has to be modeled explicitly ([Brennecke et al., Nat Methods 2013](https://www.nature.com/articles/nmeth.2645); [Kharchenko et al., Nat Methods 2014](https://www.nature.com/articles/nmeth.2967)).

Effect on outcome: zeros mix true absence with failed capture, which inflates apparent heterogeneity and the bimodality of expression.

## 3. Detection physics

### Photon shot noise and camera noise

Statistical form: photon counts ~ Poisson; read noise is roughly Gaussian and, for sCMOS sensors, specific to each pixel ([Huang et al., Nat Methods 2013](https://pubmed.ncbi.nlm.nih.gov/23708387/)).

Assays: fluorescence plate readers, microscopy, flow cytometry, luminescence at low signal.

Effect on outcome: the variance of dim signals scales with their mean, so dim conditions look noisier and fail thresholds more often. Photon counts also set the localization precision of single molecules.

Candidate: it belongs with the threshold installment ("The threshold is part of the measurement").

### Photobleaching and blinking

Statistical form: exponential waiting times; step counts ~ Binomial(n, labeling × detection efficiency).

Assays: single-molecule step counting ([Ulbrich and Isacoff, Nat Methods 2007](https://www.nature.com/articles/nmeth1024)), single-molecule localization microscopy, FRAP, long time-lapse imaging.

Effect on outcome: incomplete maturation or labeling of fluorophores biases stoichiometry downward in a binomially predictable way. Blinking inflates counts in localization microscopy.

Candidate: strong. The binomial correction for subunit counting is compact and often omitted.

## 4. Single-cell biology

### Transcriptional bursting and intrinsic or extrinsic noise

Statistical form: two-state telegraph model. mRNA counts follow a negative binomial distribution, with burst frequency and size as its parameters. Elowitz and colleagues separated intrinsic from extrinsic noise with dual reporters ([Science 2002](https://pubmed.ncbi.nlm.nih.gov/12183631/)), reviewed by [Raj and van Oudenaarden (Cell 2008)](https://pubmed.ncbi.nlm.nih.gov/18957198/).

Assays: reporter assays, single-cell expression assays, any threshold-based readout of cell fraction (percent positive, percent responding).

Effect on outcome: a population mean can shift because burst frequency changed (more cells on) or because burst size changed (brighter cells), and the two mechanisms imply different pharmacology. Percent-positive readouts depend on where the threshold sits within a skewed distribution.

Candidate: strong. Frequency versus size from mean and Fano factor is inspectable math.

### Variability and memory of protein levels

Statistical form: protein levels fluctuate around a stationary distribution, with a mixing time of one or more generations. Sigal and colleagues found memories of protein level lasting more than two generations ([Nature 2006](https://pubmed.ncbi.nlm.nih.gov/17122776/)).

Assays: drug response, time-course imaging, sorting-then-assay experiments. Cohen and colleagues linked cell-to-cell differences in protein dynamics to which cells survived a drug ([Science 2008](https://pubmed.ncbi.nlm.nih.gov/19023046/)).

Effect on outcome: a sorted "high" population relaxes back to the parental distribution, so measurements depend on the time since sorting.

Candidate: strong. Regression to the mean after sorting has a closed form under an Ornstein–Uhlenbeck model and is a frequent source of false "induced" effects.

### Fractional killing and fate variability

Statistical form: thresholded response to a variable protein state, often treated as a hazard model.

Assays: apoptosis and cytotoxicity assays. Spencer and colleagues traced variability in TRAIL-induced death to differences in protein levels between cells rather than genetic differences ([Nature 2009](https://pubmed.ncbi.nlm.nih.gov/19363473/)).

Effect on outcome: sigmoidal dose-response curves can plateau below 100 percent killing for non-genetic reasons, and a second dose kills less than the first predicts because survivors are enriched for resistant states that later relax.

Candidate: very strong. "The survivors are not resistant, they are sampled" ties to IC50 and Emax interpretation.

### Phenotypic switching, persisters and rare pre-resistant states

Statistical form: two-state Markov switching with low rates. Kill curves are biphasic, with a sum of exponentials.

Assays: antibiotic time-kill and MIC assays, where [Balaban and colleagues (Science 2004)](https://pubmed.ncbi.nlm.nih.gov/15308767/) linked persistence to phenotypic switching that existed before treatment. Cancer drug-tolerance assays, where Shaffer and colleagues showed rare transcriptional states that predicted melanoma resistance ([Nature 2017](https://pubmed.ncbi.nlm.nih.gov/28607484/)).

Effect on outcome: endpoint assays read the switching equilibrium, not the drug. Results depend on the growth phase and inoculum history.

Candidate: strong. A biphasic kill model with the switching rates recovered from the curve is a clean tool.

### Random monoallelic expression

Statistical form: random allelic choice made once per clone and inherited stably. [Gimelbrant and colleagues (Science 2007)](https://pubmed.ncbi.nlm.nih.gov/18006746/) found it widespread on human autosomes in clonal lines.

Assays: clonal cell lines, allele-specific assays, knockouts.

Effect on outcome: clones from the same parent differ in which allele is expressed, so a heterozygous mutation can look absent or dominant depending on the clone.

### Growth and cell-cycle stochasticity

Statistical form: variable division times; asynchronous cycle phase. Expression fluctuations of metabolic enzymes can propagate into growth fluctuations ([Kiviet et al., Nature 2014](https://www.nature.com/articles/nature13582)).

Assays: proliferation endpoints, phase-specific drug responses, outgrowth after single-cell deposition (the p in No. 07).

Effect on outcome: endpoint counts mix growth rate with lag-time variability. Clone-size distributions are broad even for identical genotypes.

## 5. Population and evolutionary history

### Mutation timing and jackpots (Luria–Delbrück)

Statistical form: mutants arise at random times during expansion. The resulting Luria–Delbrück distribution has a heavy right tail ([Luria and Delbrück, Genetics 1943](https://pmc.ncbi.nlm.nih.gov/articles/PMC1209226/)).

Assays: fluctuation tests, Ames-type mutagenesis, resistance frequency, emergence of escape mutants in viral and cancer cultures.

Effect on outcome: means are dominated by rare jackpot cultures. The variance between replicate cultures carries more information than the mean.

Candidate: very strong. It is a classic with an inspectable estimator (the p0 method and maximum likelihood), and it applies directly to how resistance frequency is reported in drug discovery.

### Clonal drift and genetic evolution of cell lines

Statistical form: neutral drift plus selection during passaging and bottlenecks.

Assays: any cell-based assay run across passages or laboratories. When 27 strains of MCF7 were tested against 321 compounds, drug responses differed considerably ([Ben-David et al., Nature 2018](https://www.nature.com/articles/s41586-018-0409-3)).

Effect on outcome: "the same cell line" is a population whose composition is a random walk. This adds variation between laboratories that replication within a laboratory cannot reveal.

Candidate: strong. It pairs naturally with No. 07's closing point that clonal origin does not certify a homogeneous population.

### Integration position effects and transgene silencing

Statistical form: integration site is random, and expression and silencing rates depend on the site. Silencing can be progressive at a rate characteristic of the site ([a study of transgene position effects in PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC88803/)), and position effects at the vector insertion site contribute to silencing soon after integration ([a review of lentiviral insulation strategies](https://pmc.ncbi.nlm.nih.gov/articles/PMC3279150/)).

Assays: stable reporter lines, engineered producer lines, transduced populations.

Effect on outcome: pools drift toward low expressers over time. Clone-to-clone differences in expression often reflect the site rather than the construct.

## 6. Apparent stochasticity

Some variability between cells that looks random is deterministic in variables the assay did not record. [Snijder and colleagues (Nature 2009)](https://pelkmanslab.org/wp-content/uploads/2019/02/Snijder2009.pdf) showed that population context, including local density, cell size and position in the colony, predicts much of the variability in endocytosis and virus infection. Plate position, edge evaporation, time since passaging and dispensing order behave the same way. Treating these as noise widens confidence intervals that a covariate would narrow, and treating true intrinsic noise as a covariate effect overfits. Before modeling, the useful question is whether a dual-reporter, sibling-cell or position-randomized design can split the variance.

## 7. Worked case: antibody penetration into spheroids with a hooked titration curve

The assay. Tumor spheroids are incubated with a primary antibody over a concentration series, fixed, stained with a fluorescently labeled secondary, and imaged either as cryosections or as cleared whole mounts. The readout is signal versus depth from the rim, summarized as a penetration depth or as the core-to-rim ratio, and plotted against primary concentration. The curve rises, reaches a maximum, and then falls at the highest primary concentrations. The fall is the hook, and the first job is to decide which of several mechanisms produced it, because they carry different stochastic signatures.

### Mechanisms that produce the hook

Secondary limitation. The secondary is used at one fixed concentration while the amount of bound primary grows with the primary titration. Once bound primary exceeds the secondary available in the incubation, the secondary becomes the limiting reagent, and it is now the secondary that meets a binding-site barrier: it saturates the outer shell of primary and fails to reach the core. The primary distribution continues to deepen while the detected signal retreats to the rim. The classic binding-site barrier was described for the primary by [Fujimori and colleagues (J Nucl Med 1990)](https://pubmed.ncbi.nlm.nih.gov/2362198/) and generalized as transport opposed by binding and clearance by [Thurber, Schmidt and Wittrup (Adv Drug Deliv Rev 2008)](https://pubmed.ncbi.nlm.nih.gov/18541331/); in an indirect stain, the same mathematics applies a second time to the detection reagent. Diagnostic: the hook moves when the secondary concentration or incubation time changes, and a directly labeled primary at the same concentrations shows no hook.

Fluorophore self-quenching. At high local label density, dyes on adjacent secondaries transfer energy to each other and emit less per molecule. Signal per bound antibody falls exactly where occupancy is highest, which is the rim at high primary concentration. Multiply labeled antibodies exhibit self-quenching, and brightness per dye falls as the degree of labeling rises ([a FRET microscopy analysis in Anal Chem](https://pubs.acs.org/doi/10.1021/acs.analchem.9b01504); [a review of fluorescent tracers in PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC6527780/)). Diagnostic: the hook flattens when a lower-degree-of-labeling secondary is used, and rim intensity saturates before core intensity.

Optical attenuation. In whole mounts the signal from deep voxels passes through densely stained outer shells. When the rim is bright, excitation and emission are attenuated by absorption and scattering on their way in and out, and the apparent core signal falls as the rim signal rises. Light penetration into uncleared samples is usually limited to about 50 to 70 µm ([Nürnberg et al., Front Mol Biosci 2020](https://pmc.ncbi.nlm.nih.gov/articles/PMC7046628/)). Diagnostic: the hook is weaker in cryosections than in whole mounts, and it weakens after clearing.

Prozone-type effects. In solution-phase immunoassays, high analyte concentrations produce lower signal because excess analyte occupies capture and detection reagents separately ([Tate and Ward, Clin Biochem Rev 2004](https://pubmed.ncbi.nlm.nih.gov/18458713/)). In a washed indirect stain, the analogous route is incomplete washing: free primary left in the interstitium binds secondary and removes it before it can reach immobilized primary. Diagnostic: additional wash steps shift the hook to higher concentrations.

Antigen loss at high occupancy. Antibody-induced internalization or shedding reduces surface antigen at saturating concentrations. This is biology, not detection, but it produces the same shape. Diagnostic: it persists with a directly labeled primary and is absent on fixed spheroids stained post-fixation.

### Stochastic processes inside the assay

Seeding and spheroid size. Spheroids form from a dispensed cell count that follows Poisson or overdispersed statistics, so diameter varies between replicates even in a single plate. The observed polydispersity arises from stochastic variation in the initial seeding count, which is inherent to every formation method ([a 2025 screening study in PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC11884609/)); one hanging-drop array reported a Poisson-type size distribution with about 14 percent deviation in diameter ([a KIT hanging-drop study](https://publikationen.bibliothek.kit.edu/1000095611/34071369)). Penetration depth measured in micrometers is independent of size only in the planar limit. The core-to-rim ratio is not: a larger spheroid has more core, so a dose that saturates a 300 µm spheroid leaves a 500 µm spheroid hollow, and the hook position depends on the size distribution of the replicates.

Random section plane. A cryosection cuts a sphere at a random distance from its center. The section radius is \( \sqrt{R^2 - z^2} \), and a section at random z shows an apparent radius smaller than R with a distribution worked out by Wicksell in 1925, the corpuscle problem ([a CWI treatment of Wicksell's problem](https://ir.cwi.nl/pub/5952/5952D.pdf)). An off-equator section passes through the stained shell at a glancing angle, so a shell of true thickness t appears thicker, and the apparent core shrinks. Penetration depth read from sections is biased upward unless the section is equatorial. Choosing the largest section from a serial set is a maximum statistic and corrects most of this, but it also selects, which couples section choice to staining quality.

Antigen density per cell. Expression per cell follows a skewed distribution from transcriptional bursting and extrinsic noise (Section 4). The binding-site barrier depth depends on total antigen per unit volume along the diffusion path, so cell-to-cell variation in antigen produces patchy penetration fronts. A rim of high expressers consumes more antibody than the mean would predict. In spheroids there is also a radial gradient: proliferation, hypoxia and necrosis change expression from rim to core ([Minchinton and Tannock, Nat Rev Cancer 2006](https://www.nature.com/articles/nrc1893)), so the "antigen barrier" is not uniform.

Binding kinetics at low concentration. Association and dissociation are stochastic at the molecular level, and the mean-field diffusion-reaction equation is accurate only when many molecules are present per voxel. At the lowest primary concentrations, where the front is shallow, the number of antibodies in a 10 µm shell can be small enough that the front position fluctuates between spheroids for reasons unrelated to biology.

Secondary stoichiometry. A polyclonal secondary binds a variable number of molecules per primary, typically zero to two or three depending on epitope access and crowding. The number of secondaries per primary is a random variable whose mean decreases with crowding, which is a second route to the hook independent of depletion.

Degree of labeling. Dyes are conjugated to lysines at random. The number of dyes per secondary is approximately Poisson around the nominal degree of labeling, so a batch with mean 4 contains unlabeled, singly labeled and heavily labeled molecules. Heavily labeled molecules quench and may bind differently, so the effective brightness per antibody is a mixture, and lot-to-lot differences in this distribution shift the hook.

Wash dissociation. Unbinding during washes is a first-order stochastic process with rate k_off. For a fast-dissociating primary, the fraction retained after a wash of duration τ is \( e^{-k_{\text{off}} \tau} \) per molecule, and the core, which was last to be reached, is first to be lost.

Photon statistics and depth. Signal from deep voxels is dimmer, so its relative shot noise is larger, and depth-dependent attenuation makes the measured intensity profile a convolution of the true profile with an exponential loss. Threshold-based penetration depths, defined as the depth where intensity falls below some fraction of the rim, therefore depend on the exposure and the attenuation length as much as on the antibody.

Selection of spheroids and fields. Operators image spheroids that look round, intact and well stained. Rejecting irregular spheroids is sensible, but it removes the replicates that penetrated worst, and it couples the sample to the outcome in the same way a screen couples kept wells to founder number in No. 07.

### What the hook does to inference

A non-monotonic curve maps one signal to two concentrations. A single-point assay at a high dose can report a penetration depth that is lower than at a tenfold lower dose, and a comparison of two antibodies at one concentration can rank them in reverse of their rank at another. The remedy is to measure the curve rather than a point, to use a directly labeled primary or a secondary in large excess to remove the detection hook, and to report penetration in absolute micrometers against the spheroid radius distribution rather than as a ratio.

### Candidate installment

Strong. The piece would model antibody penetration as diffusion with binding, add a second binding-site barrier for the secondary at a fixed concentration, and show the hook appearing from detection alone, with a directly labeled primary as the control. Two closed forms are inspectable: the Thiele-modulus condition for penetration versus binding, and the Wicksell section-plane bias on apparent penetration depth. A tool would take the primary concentration series, secondary concentration, antigen density, spheroid radius distribution and degree of labeling, and return the expected curve for each mechanism, so the user can see which hook shape matches theirs. Everything can be synthetic.

## Shortlist for 99 Small Problems

Ranked by fit to the series: a standard calculation that conditions on the wrong event, a closed form, a measurable parameter, and an outcome that changes decisions.

1. Fractional killing and survivor enrichment. Dose-response plateaus and the failure of second doses, from protein-state variability.
2. Library bottlenecks and jackpotting in pooled screens and display selections. Coverage, drift variance and false depletion.
3. Luria–Delbrück jackpots in resistance-frequency assays. Why the mean misleads and what the zero-class estimator recovers.
4. Regression to the mean after sorting. Protein memory sets how quickly a sorted population returns, which mimics an induced effect.
5. Poisson limits of detection in qPCR and rare-variant assays. The copy-number floor beneath any chemistry.
6. MOI size bias in pooled screens. Cells that score carry more integrations.
7. Burst frequency versus burst size. Which pharmacology a shift in the mean implies.
8. Photobleaching step counts under incomplete labeling. The binomial correction for stoichiometry.
9. Biphasic kill curves and persister switching rates.
10. PCR as a branching process. Amplification variance and the cost of UMI collisions.
11. The hooked penetration curve. Secondary-antibody binding-site barrier, self-quenching and optical attenuation produce the same non-monotonic shape, with the Wicksell section bias on apparent depth (Section 7).

Items 1, 2, 3 and 11 are the closest structural siblings to No. 07, because each involves selection acting on a random variable that the standard readout treats as fixed.
