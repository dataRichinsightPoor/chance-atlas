# Chance Atlas

An interactive catalog of the stochastic processes known to shape bioassay outcomes. Each of the 29 entries gives the statistical form, the assays it affects, the effect on the reported result, the diagnostic that exposes it, the mitigation, and the primary literature, with a closed-form calculator evaluated in the browser.

Companion to the **99 Small Problems** series of Data-Rich, Insight-Poor. Research-use alpha, v0.1.0-alpha.

- Atlas: https://datarichinsightpoor.github.io/chance-atlas/
- Math and tests: https://datarichinsightpoor.github.io/chance-atlas/methods.html
- Engine: [`web/models.js`](web/models.js) · Entries: [`web/catalog.js`](web/catalog.js) · Calculators: [`web/calcs.js`](web/calcs.js)
- Prose catalog with candidate installments: [`docs/CATALOG.md`](docs/CATALOG.md)

## What it is for

- **Designing an experiment.** Before seeding a plate, transducing a library or titrating an antibody, read off the sampling floor, the selection bias or the detection artifact the design will carry.
- **Diagnosing a result.** Start from the assay in the index and find the processes that touch it, each with the measurement that separates it from the others.
- **Teaching.** Every entry pairs a mechanism with a closed form and a chart that responds to its parameters.

## Sections

A. Sampling and partitioning · B. Amplification · C. Detection physics · D. Single-cell biology · E. Population history · F. Apparent stochasticity · G. Worked case: spheroid antibody penetration with a hooked titration curve.

## Inputs and outputs

Each calculator takes the parameters named in its entry (for example the empty-well fraction and screen pass rate; the MOI and infectious fraction; burst frequency, size and threshold; spheroid radius, antigen density, secondary concentration and carryover). Outputs are the closed-form quantities, a chart, a per-calculator JSON copy and PNG export, and a whole-atlas JSON export. Deep links: `index.html#<entry-id>` opens an entry; `methods.html#<entry-id>` opens its derivation.

## Running locally

```
npm test                      # 23 numerical tests (closed forms, simulations, known values)
python3 tools-build-methods.py # rebuild web/methods.html from web/MATH.md
node analysis/examples.mjs    # evaluate every calculator at its defaults
cd web && python3 -m http.server 8080
```

No build step and no dependencies; the site is static ES modules.

## Boundaries

The calculators evaluate idealized models with user-supplied parameters; none is estimated from data. The catalog lists well-documented processes with primary references, not every source of randomness in biology. Passing tests shows that the code implements the stated equations, not that any parameter describes a real system. The implementation is original code written from the public literature cited in each entry. No employer-provided data, constructs, methods or internal results were used.

## License

MIT. See [LICENSE](LICENSE).
