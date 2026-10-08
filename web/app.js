import { SECTIONS, ENTRIES, ASSAYS } from './catalog.js';
import { CALCS } from './calcs.js';
import { VERSION } from './models.js';

const $ = (id) => document.getElementById(id);
const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
const state = { q: '', section: 'all', assay: 'all', calcOnly: false, open: new Set() };
const values = {}; // calc id -> current params

/* ------------------------------------------------------------- plotting */
function plot(canvas, chart) {
  const dpr = window.devicePixelRatio || 1, W = canvas.clientWidth || 640, H = canvas.clientHeight || 300;
  canvas.width = W * dpr; canvas.height = H * dpr;
  const ctx = canvas.getContext('2d'); ctx.setTransform(dpr, 0, 0, dpr, 0, 0); ctx.clearRect(0, 0, W, H);
  const colors = [css('--teal'), css('--amber'), css('--floor'), css('--muted')], ink = css('--ink'), muted = css('--muted'), line = css('--line');
  ctx.font = '12px ' + css('--sans');
  const legendW = chart.series.reduce((a, s) => a + ctx.measureText(s.name).width + 36, 0), legendRows = Math.max(1, Math.ceil(legendW / (W - 78)));
  const pad = { l: 62, r: 34, t: 10 + 18 * legendRows, b: 44 }, pw = W - pad.l - pad.r, ph = H - pad.t - pad.b;
  const xs = chart.x, logX = !!chart.logX, logY = !!chart.logY;
  const fx = (v) => logX ? Math.log10(v) : v, fy = (v) => logY ? Math.log10(v) : v;
  let ys = chart.series.flatMap(s => s.y).filter(v => Number.isFinite(v) && (!logY || v > 0));
  if (chart.hline != null) ys.push(chart.hline);
  let ymin = chart.yMin ?? Math.min(...ys), ymax = Math.max(...ys);
  if (!logY) { if (ymin > 0 && ymin < 0.25 * ymax) ymin = 0; const span = ymax - ymin || 1; ymax += 0.06 * span; if (ymin < 0) ymin -= 0.04 * span; }
  else { ymin = Math.max(ymin, Math.min(...ys)); ymin = fy(ymin); ymax = fy(ymax) + 0.1; }
  const xv = xs.filter(v => Number.isFinite(v) && (!logX || v > 0)); let xmin = fx(Math.min(...xv)), xmax = fx(Math.max(...xv));
  if (chart.bars) { xmin -= 0.6; xmax += 0.6; }
  const X = (v) => pad.l + (fx(v) - xmin) / (xmax - xmin || 1) * pw, Y = (v) => pad.t + ph - ((logY ? fy(v) : v) - ymin) / (ymax - ymin || 1) * ph;
  ctx.font = '12px ' + css('--sans'); ctx.fillStyle = muted; ctx.strokeStyle = line; ctx.lineWidth = 1;
  // y ticks
  const yticks = ticks(ymin, ymax, logY), fmt = (v) => Math.abs(v) >= 1000 ? v.toLocaleString('en-US', { maximumFractionDigits: 0 }) : (+v.toPrecision(3)).toString();
  ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
  for (const t of yticks) { const y = pad.t + ph - (t - ymin) / (ymax - ymin || 1) * ph; ctx.beginPath(); ctx.moveTo(pad.l, y); ctx.lineTo(W - pad.r, y); ctx.globalAlpha = 0.5; ctx.stroke(); ctx.globalAlpha = 1; ctx.fillText(fmt(logY ? 10 ** t : t), pad.l - 8, y); }
  const xticks = ticks(xmin, xmax, logX); ctx.textAlign = 'center'; ctx.textBaseline = 'top';
  for (const t of xticks) { const x = pad.l + (t - xmin) / (xmax - xmin || 1) * pw; ctx.fillText(fmt(logX ? 10 ** t : t), x, pad.t + ph + 6); }
  ctx.fillStyle = ink; ctx.fillText(chart.xLabel, pad.l + pw / 2, pad.t + ph + 24);
  ctx.save(); ctx.translate(14, pad.t + ph / 2); ctx.rotate(-Math.PI / 2); ctx.textBaseline = 'top'; ctx.fillText(chart.yLabel, 0, 0); ctx.restore();
  if (chart.hline != null) { ctx.setLineDash([3, 4]); ctx.strokeStyle = muted; ctx.beginPath(); ctx.moveTo(pad.l, Y(chart.hline)); ctx.lineTo(W - pad.r, Y(chart.hline)); ctx.stroke(); ctx.setLineDash([]); }
  if (chart.mark != null && Number.isFinite(chart.mark) && (!logX || chart.mark > 0)) { ctx.setLineDash([4, 4]); ctx.strokeStyle = muted; ctx.beginPath(); ctx.moveTo(X(chart.mark), pad.t); ctx.lineTo(X(chart.mark), pad.t + ph); ctx.stroke(); ctx.setLineDash([]); }
  chart.series.forEach((s, i) => {
    ctx.strokeStyle = ctx.fillStyle = colors[i % colors.length]; ctx.lineWidth = i === 0 ? 2.4 : 1.8;
    if (chart.bars) { const bw = pw / (xs.length + 1) * 0.7; s.y.forEach((v, j) => { if (Number.isFinite(v)) ctx.fillRect(X(xs[j]) - bw / 2, Y(v), bw, Y(ymin < 0 ? 0 : ymin) - Y(v)); }); return; }
    ctx.beginPath(); let pen = false;
    s.y.forEach((v, j) => { const ok = Number.isFinite(v) && (!logY || v > 0) && (!logX || xs[j] > 0); if (!ok) { pen = false; return; } const x = X(xs[j]), y = Math.max(pad.t, Math.min(pad.t + ph, Y(v))); if (!pen) { ctx.moveTo(x, y); pen = true; } else ctx.lineTo(x, y); });
    ctx.stroke();
  });
  // legend
  ctx.textAlign = 'left'; ctx.textBaseline = 'middle'; let lx = pad.l, ly = 12;
  chart.series.forEach((s, i) => { const w = ctx.measureText(s.name).width + 36; if (lx + w > W - pad.r && lx > pad.l) { lx = pad.l; ly += 18; } ctx.fillStyle = colors[i % colors.length]; ctx.fillRect(lx, ly - 1, 14, 3); ctx.fillStyle = ink; ctx.fillText(s.name, lx + 20, ly); lx += w; });
}
function ticks(a, b, log) {
  if (log) { const o = []; for (let e = Math.floor(a); e <= Math.ceil(b); e++) if (e >= a - 1e-9 && e <= b + 1e-9) o.push(e); if (o.length < 2) { o.length = 0; for (let e = Math.floor(a); e <= Math.ceil(b); e++) for (const m of [1, 2, 5]) { const v = e + Math.log10(m); if (v >= a && v <= b) o.push(v); } } return o; }
  const span = b - a || 1, raw = span / 5, mag = 10 ** Math.floor(Math.log10(raw)), step = [1, 2, 2.5, 5, 10].map(m => m * mag).find(s => span / s <= 6) || mag;
  const o = []; for (let v = Math.ceil(a / step) * step; v <= b + 1e-9; v += step) o.push(+v.toFixed(10)); return o;
}

/* --------------------------------------------------------------- render */
function matches(e) {
  if (state.section !== 'all' && e.section !== state.section) return false;
  if (state.assay !== 'all' && !e.assays.includes(state.assay)) return false;
  if (state.calcOnly && !e.calc) return false;
  if (state.q) { const q = state.q.toLowerCase(); const hay = [e.title, e.one, e.form, e.effect, e.diagnostic, e.mitigation, ...e.assays].join(' ').toLowerCase(); if (!hay.includes(q)) return false; }
  return true;
}
function renderList() {
  const root = $('entries'); root.innerHTML = '';
  let shown = 0;
  for (const S of SECTIONS) {
    const items = ENTRIES.filter(e => e.section === S.id && matches(e));
    if (!items.length) continue;
    const sec = document.createElement('section'); sec.className = 'atlas-section'; sec.id = 'section-' + S.id;
    sec.innerHTML = `<p class="eyebrow">${S.id} <span>/</span> ${esc(S.name).toUpperCase()}</p><p class="blurb">${esc(S.blurb)}</p>`;
    for (const e of items) { sec.appendChild(card(e)); shown++; }
    root.appendChild(sec);
  }
  $('count').textContent = `${shown} of ${ENTRIES.length} entries`;
}
function card(e) {
  const d = document.createElement('details'); d.className = 'entry'; d.id = e.id; d.open = state.open.has(e.id);
  d.innerHTML = `<summary><span class="t">${esc(e.title)}</span><span class="one">${esc(e.one)}</span><span class="chips">${e.assays.map(a => `<i>${esc(a)}</i>`).join('')}${e.calc ? '<b>calculator</b>' : ''}</span></summary>
  <div class="body">
    <h4>Statistical form</h4><p class="form">${e.form}</p>
    <h4>Effect on the reported result</h4><p>${e.effect}</p>
    <h4>Diagnostic</h4><p>${e.diagnostic}</p>
    <h4>Mitigation</h4><p>${e.mitigation}</p>
    ${e.calc ? `<div class="calc" data-calc="${e.calc}" data-entry="${e.id}"></div>` : ''}
    <h4>References</h4><ol class="refs">${e.refs.map(r => `<li>${esc(r.t)} <a href="${r.u}">${r.u.replace(/^https?:\/\//, '')}</a></li>`).join('')}</ol>
    ${e.tool ? `<p class="hint">Full model in the series: <a href="${e.tool.url}">${esc(e.tool.name)}</a>.</p>` : ''}
    <p class="hint"><a href="#${e.id}" class="permalink">Link to this entry</a> · <a href="methods.html#${e.id}">Derivation</a></p>
  </div>`;
  d.addEventListener('toggle', () => { if (d.open) { state.open.add(e.id); mountCalc(d.querySelector('.calc')); } else state.open.delete(e.id); });
  if (d.open) queueMicrotask(() => mountCalc(d.querySelector('.calc')));
  return d;
}
function mountCalc(box) {
  if (!box || box.dataset.mounted) return; box.dataset.mounted = '1';
  const id = box.dataset.calc, C = CALCS[id]; values[id] ??= Object.fromEntries(C.params.map(p => [p.key, p.def]));
  box.innerHTML = `<div class="calc-head"><h4>Calculator</h4><span class="badge">CLOSED FORM · v${VERSION}</span></div>
  <div class="calc-grid"><form class="params${C.params.length > 6 ? " wide" : ""}" novalidate>${C.params.map(p => `<label>${esc(p.label)}<input type="number" data-key="${p.key}" value="${p.def}" min="${p.min}" max="${p.max}" step="${p.step}"></label>`).join('')}
  <div class="calc-actions"><button type="button" class="reset">Reset</button><button type="button" class="json">Copy JSON</button><button type="button" class="png">Save PNG</button><span class="hint status"></span></div></form>
  <div class="out"><dl class="outputs"></dl><canvas height="300" role="img" aria-label="Calculator chart"></canvas></div></div>`;
  const form = box.querySelector('form'), dl = box.querySelector('.outputs'), cv = box.querySelector('canvas'), st = box.querySelector('.status');
  let last = null;
  const run = () => {
    const p = {}; for (const inp of form.querySelectorAll('input')) { const v = Number(inp.value); p[inp.dataset.key] = Number.isFinite(v) ? v : values[id][inp.dataset.key]; }
    values[id] = p;
    try { last = C.compute(p); dl.innerHTML = last.outputs.map(([k, v]) => `<div><dt>${esc(k)}</dt><dd>${esc(v)}</dd></div>`).join(''); plot(cv, last.chart); st.textContent = ''; }
    catch (err) { st.textContent = err.message; }
  };
  form.addEventListener('input', run);
  form.querySelector('.reset').addEventListener('click', () => { for (const inp of form.querySelectorAll('input')) inp.value = C.params.find(q => q.key === inp.dataset.key).def; run(); });
  form.querySelector('.json').addEventListener('click', async () => { const payload = { tool: 'Chance Atlas', version: VERSION, entry: box.dataset.entry, calculator: id, params: values[id], outputs: Object.fromEntries(last.outputs), chart: last.chart }; try { await navigator.clipboard.writeText(JSON.stringify(payload, null, 2)); st.textContent = 'Copied.'; } catch { st.textContent = 'Clipboard unavailable; use the browser console.'; console.log(payload); } });
  form.querySelector('.png').addEventListener('click', () => { const a = document.createElement('a'); a.download = `chance-atlas-${box.dataset.entry}.png`; a.href = cv.toDataURL('image/png'); a.click(); });
  run();
  new ResizeObserver(() => last && plot(cv, last.chart)).observe(cv);
}
function renderAssayIndex() {
  const root = $('assay-index'); root.innerHTML = '';
  for (const a of ASSAYS) {
    const es = ENTRIES.filter(e => e.assays.includes(a));
    const div = document.createElement('div'); div.innerHTML = `<h4>${esc(a)}</h4><p>${es.map(e => `<a href="#${e.id}" data-entry="${e.id}">${esc(e.title)}</a>`).join(' · ')}</p>`; root.appendChild(div);
  }
}

/* --------------------------------------------------------------- wiring */
function openFromHash() {
  const id = location.hash.slice(1); if (!id) return;
  const e = ENTRIES.find(x => x.id === id); if (!e) return;
  state.section = 'all'; state.assay = 'all'; state.q = ''; state.calcOnly = false; $('q').value = ''; $('section').value = 'all'; $('assay').value = 'all'; $('calcOnly').checked = false;
  state.open.add(id); renderList(); requestAnimationFrame(() => document.getElementById(id)?.scrollIntoView({ block: 'start' }));
}
function init() {
  $('section').innerHTML = '<option value="all">All sections</option>' + SECTIONS.map(s => `<option value="${s.id}">${s.id} · ${esc(s.name)}</option>`).join('');
  $('assay').innerHTML = '<option value="all">All assays</option>' + ASSAYS.map(a => `<option value="${esc(a)}">${esc(a)}</option>`).join('');
  $('q').addEventListener('input', (ev) => { state.q = ev.target.value.trim(); renderList(); });
  $('section').addEventListener('change', (ev) => { state.section = ev.target.value; renderList(); });
  $('assay').addEventListener('change', (ev) => { state.assay = ev.target.value; renderList(); });
  $('calcOnly').addEventListener('change', (ev) => { state.calcOnly = ev.target.checked; renderList(); });
  $('expand').addEventListener('click', () => { ENTRIES.forEach(e => matches(e) && state.open.add(e.id)); renderList(); });
  $('collapse').addEventListener('click', () => { state.open.clear(); renderList(); });
  $('theme').addEventListener('click', () => { const r = document.documentElement, light = r.dataset.theme !== 'light'; r.dataset.theme = light ? 'light' : 'dark'; $('theme').textContent = light ? 'Dark mode' : 'Light mode'; document.querySelectorAll('.calc[data-mounted] canvas').forEach(c => c.dispatchEvent(new Event('resize'))); renderList(); });
  $('export-all').addEventListener('click', async () => {
    const out = {}; for (const [id, C] of Object.entries(CALCS)) { const p = values[id] ?? Object.fromEntries(C.params.map(q => [q.key, q.def])); const r = C.compute(p); out[id] = { params: p, outputs: Object.fromEntries(r.outputs) }; }
    const blob = new Blob([JSON.stringify({ tool: 'Chance Atlas', version: VERSION, generated: new Date().toISOString(), calculators: out }, null, 2)], { type: 'application/json' });
    const a = document.createElement('a'); a.download = `chance-atlas-${VERSION}.json`; a.href = URL.createObjectURL(blob); a.click();
  });
  window.addEventListener('hashchange', openFromHash);
  renderAssayIndex(); renderList(); openFromHash();
  $('version').textContent = 'v' + VERSION;
}
init();
