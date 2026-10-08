// Evaluates every calculator at its default parameters and writes analysis/examples-results.json.
import { CALCS } from '../web/calcs.js';
import { VERSION } from '../web/models.js';
import { writeFileSync } from 'node:fs';
const out = { version: VERSION, calculators: {} };
for (const [id, C] of Object.entries(CALCS)) {
  const p = Object.fromEntries(C.params.map(q => [q.key, q.def]));
  const r = C.compute(p);
  out.calculators[id] = { params: p, outputs: Object.fromEntries(r.outputs) };
}
writeFileSync(new URL('./examples-results.json', import.meta.url), JSON.stringify(out, null, 2));
console.log(Object.keys(out.calculators).length, 'calculators evaluated');
