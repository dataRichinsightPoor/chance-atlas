# Verification

`npm test` runs 23 tests with Node's built-in runner against `web/models.js`. Simulations use a seeded generator and are reproducible.

| Test | Checks |
|---|---|
| Special functions | Γ(5) = 24; Φ(1.96); Φ⁻¹(0.975); bivariate normal at the origin against 1/4 + asin(ρ)/2π |
| Poisson partitions | λ, random and selected multi-unit risks at f₀ = 0.7; ratio → 2 at low density |
| MOI | p₀ = e^−λ; multi-copy share among scoring cells exceeds share among transduced; pmf sums to 1 |
| Kärber design | SE ∝ 1/√(n−1); larger dilution step gives larger SE |
| Low copy | P(detect) = 1 − e^−c; LoD95 = 2.996 |
| Bottlenecks | CV ≈ √(k/c); one-stage loss = e^−c |
| PCR branching | closed-form mean and CV against 4,000 Galton–Watson realizations |
| UMI collisions | U(1 − e^−N/U) |
| Allele and capture dropout | closed forms |
| Shot noise, bleaching | SNR = √N; binomial normalization, pⁿ, mode |
| Number fluctuations, Berg–Purcell | ⟨N⟩ at 1 nM in 1 fL; error ∝ 1/√T; monitor/absorber variance ratio 12/5; absorber error = 1/√arrivals |
| Channel noise | variance at p = ½; parabola apex |
| Bursting | negative binomial normalization, mean rb, Fano 1 + b |
| Sort relaxation | truncated-normal mean; e^−t/τ decay |
| Fractional killing | two-dose survival against 60,000 correlated-normal pairs |
| Persisters, switching | closed forms; asymptote |
| Luria–Delbrück | p₀, p₁, truncated mass, heavier tail than Poisson |
| Clonal drift | variance against p(1−p)(1 − (1 − 1/N)^k) at s = 0 |
| Mitochondrial bottleneck, ELDA | closed form; SE halves with 4× replicates |
| Hierarchical SE | shares sum to 1; true SE > naive |
| Spheroid hook | ideal, secondary-limited and direct-label curves monotonic; limitation caps signal; carryover, crowding and quenching each produce a fall |
| Wicksell | mean section radius = πR/4; fully stained fraction = t/R; apparent thickness against 200,000-sample Monte Carlo |
| Seeding size | CV_R = CV_count/3 |

The browser build was checked in headless Chromium: no console errors, all 29 entries render, filters and deep links work, calculators mount on open, and the layout holds at 390 px width.

These tests show that the code implements the equations in `web/MATH.md`. They are not biological validation.
