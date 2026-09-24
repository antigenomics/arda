"""The QC dashboard: one HTML file with the data inside it and no external reference at all.

Never: **no CDN, no bundle, no dependency.** arda has no plotting or templating library and gains
none here. The QC JSON is inlined and a few hundred lines of plain JavaScript draw SVG, so the
page opens on an air-gapped cluster login node, off a USB stick, or as an email attachment, and
it still works after the results directory is deleted. A chart library would be richer and would
cost exactly the property the file exists for. That is the same stance as :mod:`arda.scplot`,
whose gnuplot script is written whether or not gnuplot exists.

Colours are :data:`arda.scplot.DARK2` (ColorBrewer, colour-blind safe) over a transparent
background with one mid-grey ink, so the page reads the same light, dark and printed -- the rule
``scplot`` already states, kept here rather than restated as a second palette.

Takes either a batch document from :func:`arda.qc.write_batch` or one sample's ``.stats.json``:
a single run is a cohort of one, and the panels are the same.
"""

from __future__ import annotations

import json
from pathlib import Path

from .scplot import DARK2, INK

__all__ = ["render", "build"]

#: Distribution scopes, in the order the page offers them. Every one is keyed ``locus:bucket``,
#: which is why one chart draws all of them.
_DISTRIBUTIONS = [
    ("junction_aa_len", "Junction length (aa)"),
    ("read_len", "Aligned read length (nt)"),
    ("clone_size", "Clonotype size (reads)"),
    ("chain_support", "Molecules per chain"),
    ("isotype", "Isotype"),
    ("v_gene", "V gene usage"),
    ("j_gene", "J gene usage"),
]


def _as_batch(doc: dict, name: str) -> dict:
    """Accept a single sample's ``.stats.json`` where a batch document is expected.

    One run is a cohort of one. Rather than a second renderer, the single-sample shape is lifted
    into the batch shape here -- so every panel below has exactly one input format to read.
    """
    if "samples" in doc and "stats" in doc:
        return doc
    return {"samples": [name], "labels": {name: {"project": "", "batch": ""}},
            "stats": {name: doc.get("stats", {})}, "outliers": [], "z_flag": None,
            "min_group": None}


def build(doc: dict, *, title: str = "arda QC") -> str:
    """The whole page as one string. Separated from :func:`render` so a test can skip the disk."""
    # Never: substitution, not `str.format` -- the CSS and the script are full of braces, and
    # `{` in a format string is not a brace. `</` inside the payload is escaped because a `</script`
    # anywhere in a JSON string ends the block early, whatever the JSON says.
    payload = json.dumps(doc, sort_keys=True).replace("</", "<\\/")
    js = (_JS.replace("__PALETTE__", json.dumps(DARK2))
             .replace("__DISTRIBUTIONS__", json.dumps(_DISTRIBUTIONS)))
    return (_PAGE.replace("__TITLE__", _escape(title)).replace("__CSS__", _CSS)
            .replace("__PAYLOAD__", payload).replace("__JS__", js))


def render(input: str | Path, output: str | Path, *, title: str | None = None) -> Path:
    """Render ``input`` (a ``.qc.json`` or a ``.stats.json``) to ``output``. Returns ``output``."""
    input, output = Path(input), Path(output)
    doc = _as_batch(json.loads(input.read_text()), input.name.split(".")[0])
    output.write_text(build(doc, title=title or f"arda QC — {input.name.split('.')[0]}"))
    return output


def _escape(text: str) -> str:
    return (text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


_CSS = """
:root { color-scheme: light dark; }
body { font: 14px/1.5 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
       margin: 0 auto; padding: 1.5rem; max-width: 1100px; background: transparent; }
h1 { font-size: 1.4rem; margin: 0 0 .25rem; }
h2 { font-size: 1.05rem; margin: 2rem 0 .5rem; border-bottom: 1px solid #8884; padding-bottom: .2rem; }
.sub { color: INKCOLOR; margin: 0 0 1rem; }
table { border-collapse: collapse; font-variant-numeric: tabular-nums; font-size: 13px; }
th, td { padding: .25rem .5rem; text-align: right; white-space: nowrap; }
th:first-child, td:first-child { text-align: left; }
th { cursor: pointer; user-select: none; border-bottom: 1px solid #8886; }
th.sorted::after { content: " \\25be"; }
th.sorted.asc::after { content: " \\25b4"; }
tbody tr:hover { background: #8881; }
.scroll { overflow-x: auto; }
.controls { display: flex; flex-wrap: wrap; gap: .5rem; align-items: center; margin: .5rem 0; }
select, input { font: inherit; padding: .15rem .3rem; }
.flag { outline: 2px solid #d95f02; outline-offset: -2px; font-weight: 600; }
.legend { display: flex; flex-wrap: wrap; gap: .75rem; margin: .4rem 0; font-size: 12px; }
.legend span { display: inline-flex; align-items: center; gap: .3rem; cursor: pointer; }
.legend i { width: .75rem; height: .75rem; border-radius: 2px; display: inline-block; }
.legend .off { opacity: .35; }
.empty { color: INKCOLOR; font-style: italic; }
svg text { fill: INKCOLOR; font-size: 11px; }
.tip { position: fixed; pointer-events: none; background: #222; color: #eee; padding: .2rem .45rem;
       border-radius: 3px; font-size: 12px; opacity: 0; transition: opacity .1s; }
dl { display: grid; grid-template-columns: max-content 1fr; gap: .1rem .75rem; margin: 0; }
dt { color: INKCOLOR; }
dd { margin: 0; }
""".replace("INKCOLOR", INK)


# The whole renderer. Plain DOM and hand-drawn SVG: no framework, no build step, and nothing the
# page has to fetch. Kept in one string rather than a package data file so a wheel cannot ship a
# page without its script.
_JS = r"""
const DOC = JSON.parse(document.getElementById('qc-data').textContent);
const PALETTE = __PALETTE__, DISTS = __DISTRIBUTIONS__;
const Z_FLAG = DOC.z_flag ?? 3.5, MIN_GROUP = DOC.min_group ?? 5;

// The columns worth seeing first. Anything else is one checkbox away; showing all ~90 metrics by
// default makes the one that matters unfindable.
const FRONT = ['map.total_reads', 'map.mapped_reads', 'map.mapped_fraction', 'reads',
  'reads_with_junction', 'reads_productive', 'reads_truncated_junction', 'shm_rate',
  'clonotypes', 'clonotype_reads', 'clonotypes_chimeric', 'v_gene_coverage_reads',
  'cells', 'chains', 'pairing_rate', 'doublet_rate', 'molecules_placed_fraction',
  'map.wall_seconds', 'map.peak_rss_mb'];

const $ = (sel) => document.querySelector(sel);
// Never: `createElementNS(null, 'tbody')` is NOT an HTML tbody -- it is an element in the null
// namespace, which the layout engine draws as nothing. The table filled correctly, had six rows
// in the DOM, and rendered as a bare header. HTML goes through createElement; only SVG needs a
// namespace.
const SVG_NS = 'http://www.w3.org/2000/svg';
const SVG_TAGS = new Set(['svg', 'g', 'rect', 'line', 'text', 'path', 'circle']);
const el = (tag, attrs = {}, kids = []) => {
  const n = SVG_TAGS.has(tag) ? document.createElementNS(SVG_NS, tag)
                              : document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === 'class') n.setAttribute('class', v);
    else if (k === 'text') n.textContent = v;
    else if (k.startsWith('on')) n.addEventListener(k.slice(2), v);
    else if (v !== null && v !== undefined) n.setAttribute(k, v);
  }
  for (const kid of [].concat(kids)) n.appendChild(typeof kid === 'string'
    ? document.createTextNode(kid) : kid);
  return n;
};

// ── the table a sample is flattened into ────────────────────────────────────────────────────
function flatten(stats) {
  const row = {};
  for (const [key, metrics] of Object.entries((stats || {}).sample || {}))
    for (const [m, v] of Object.entries(metrics)) row[key ? key + '.' + m : m] = v;
  for (const [key, metrics] of Object.entries((stats || {}).run || {}))
    for (const [m, v] of Object.entries(metrics)) row[key ? key + '.' + m : m] = v;
  return row;
}
const ROWS = Object.fromEntries(DOC.samples.map((s) => [s, flatten(DOC.stats[s])]));
const LABEL = (s) => DOC.labels[s] || { project: '', batch: '' };
const GROUP = (s) => LABEL(s).project + '\u0000' + LABEL(s).batch;

const median = (xs) => {
  const v = xs.slice().sort((a, b) => a - b);
  if (!v.length) return null;
  const h = v.length >> 1;
  return v.length % 2 ? v[h] : (v[h - 1] + v[h]) / 2;
};

// Robust z, computed HERE rather than shipped per cell: it keeps the JSON small and guarantees
// the shading agrees with the numbers on screen. Same rule as `arda.qc.outliers` -- no z for a
// group under MIN_GROUP or a metric with zero MAD, because an invented scale flags whichever
// sample differs by a rounding error.
function zScores(metric) {
  const out = {}, byGroup = {};
  for (const s of DOC.samples) {
    const v = ROWS[s][metric];
    if (typeof v === 'number' && isFinite(v)) (byGroup[GROUP(s)] ||= []).push([s, v]);
  }
  for (const pairs of Object.values(byGroup)) {
    const med = median(pairs.map((p) => p[1]));
    const mad = median(pairs.map((p) => Math.abs(p[1] - med)));
    for (const [s, v] of pairs) {
      out[s] = { median: med, mad: mad,
                 z: (pairs.length >= MIN_GROUP && mad > 0) ? 0.6745 * (v - med) / mad : null };
    }
  }
  return out;
}

const fmt = (v) => typeof v !== 'number' ? (v ?? '')
  : Number.isInteger(v) ? v.toLocaleString('en-US')
  : Math.abs(v) >= 1000 ? v.toFixed(0) : v.toPrecision(4).replace(/\.?0+$/, '');

// ── panel 1: the sample table ───────────────────────────────────────────────────────────────
let sortCol = 'sample', sortAsc = true, showAll = false;

function columns() {
  const seen = new Set();
  for (const s of DOC.samples) for (const k of Object.keys(ROWS[s])) seen.add(k);
  if (showAll) return [...seen].sort();
  const front = FRONT.filter((c) => seen.has(c));
  return front.length ? front : [...seen].sort();
}

function visibleSamples() {
  const proj = $('#f-project').value, bat = $('#f-batch').value;
  const q = $('#f-text').value.trim().toLowerCase();
  return DOC.samples.filter((s) =>
    (!proj || LABEL(s).project === proj) && (!bat || LABEL(s).batch === bat) &&
    (!q || s.toLowerCase().includes(q)));
}

function drawTable() {
  const cols = columns(), zs = Object.fromEntries(cols.map((c) => [c, zScores(c)]));
  const rows = visibleSamples().sort((a, b) => {
    const x = sortCol === 'sample' ? a : ROWS[a][sortCol], y = sortCol === 'sample' ? b : ROWS[b][sortCol];
    if (x === y) return a.localeCompare(b);
    if (x === undefined || x === null) return 1;
    if (y === undefined || y === null) return -1;
    return (x > y ? 1 : -1) * (sortAsc ? 1 : -1);
  });
  const head = el('tr', {}, ['sample', 'project', 'batch', ...cols].map((c) =>
    el('th', { class: c === sortCol ? 'sorted ' + (sortAsc ? 'asc' : '') : '', title: c,
               onclick: () => { sortAsc = c === sortCol ? !sortAsc : true; sortCol = c; drawTable(); } },
      c)));
  const body = rows.map((s) => el('tr', {}, [
    el('td', { text: s }), el('td', { text: LABEL(s).project }), el('td', { text: LABEL(s).batch }),
    ...cols.map((c) => {
      const z = (zs[c][s] || {}).z;
      const flagged = z !== null && z !== undefined && Math.abs(z) >= Z_FLAG;
      return el('td', { class: flagged ? 'flag' : '',
                        title: z == null ? '' : 'robust z ' + z.toFixed(2),
                        text: fmt(ROWS[s][c]) });
    })]));
  $('#table').replaceChildren(el('table', {}, [el('thead', {}, [head]), el('tbody', {}, body)]));
  $('#table-note').textContent = rows.length + ' of ' + DOC.samples.length + ' samples; '
    + cols.length + ' metrics. Outlined cells are |robust z| ≥ ' + Z_FLAG
    + ' within their project/batch.';
}

// ── SVG helpers ─────────────────────────────────────────────────────────────────────────────
const tip = el('div', { class: 'tip' });
document.body.appendChild(tip);
const hover = (node, text) => {
  node.addEventListener('mousemove', (e) => {
    tip.textContent = text; tip.style.opacity = 1;
    tip.style.left = (e.clientX + 12) + 'px'; tip.style.top = (e.clientY + 12) + 'px';
  });
  node.addEventListener('mouseleave', () => { tip.style.opacity = 0; });
  return node;
};

function axes(w, h, pad, xLabel, yMax, ticks) {
  const g = el('g', {});
  g.appendChild(el('line', { x1: pad.l, y1: h - pad.b, x2: w - pad.r, y2: h - pad.b,
                             stroke: 'currentColor', 'stroke-opacity': .35 }));
  g.appendChild(el('line', { x1: pad.l, y1: pad.t, x2: pad.l, y2: h - pad.b,
                             stroke: 'currentColor', 'stroke-opacity': .35 }));
  for (let i = 0; i <= 4; i++) {
    const v = yMax * i / 4, y = h - pad.b - (h - pad.b - pad.t) * i / 4;
    g.appendChild(el('text', { x: pad.l - 6, y: y + 4, 'text-anchor': 'end', text: fmt(v) }));
    if (i) g.appendChild(el('line', { x1: pad.l, y1: y, x2: w - pad.r, y2: y,
                                      stroke: 'currentColor', 'stroke-opacity': .12 }));
  }
  for (const [x, label] of ticks)
    g.appendChild(el('text', { x: x, y: h - pad.b + 14, 'text-anchor': 'middle', text: label }));
  if (xLabel) g.appendChild(el('text', { x: (pad.l + w - pad.r) / 2, y: h - 2,
                                         'text-anchor': 'middle', text: xLabel }));
  return g;
}

// ── panel 2: one metric across the cohort ───────────────────────────────────────────────────
function drawMetric() {
  const metric = $('#m-pick').value, samples = visibleSamples()
    .filter((s) => typeof ROWS[s][metric] === 'number');
  const host = $('#metric');
  if (!samples.length) {
    host.replaceChildren(el('p', { class: 'empty', text: 'No sample carries this metric.' }));
    return;
  }
  const zs = zScores(metric);
  const w = 1040, h = 300, pad = { l: 70, r: 10, t: 10, b: 46 };
  const yMax = Math.max(...samples.map((s) => ROWS[s][metric])) || 1;
  const bw = (w - pad.l - pad.r) / samples.length;
  const y = (v) => h - pad.b - (h - pad.b - pad.t) * (v / yMax);
  const groups = [...new Set(samples.map(GROUP))];
  const svg = el('svg', { viewBox: `0 0 ${w} ${h}`, width: '100%', height: h });
  const ticks = samples.length <= 40
    ? samples.map((s, i) => [pad.l + bw * (i + .5), s.length > 12 ? s.slice(0, 11) + '…' : s])
    : [];
  svg.appendChild(axes(w, h, pad, metric, yMax, ticks));
  samples.forEach((s, i) => {
    const v = ROWS[s][metric], z = (zs[s] || {}).z;
    const flagged = z !== null && z !== undefined && Math.abs(z) >= Z_FLAG;
    svg.appendChild(hover(el('rect', {
      x: pad.l + bw * i + 1, y: y(v), width: Math.max(bw - 2, 1), height: h - pad.b - y(v),
      fill: flagged ? '#d95f02' : PALETTE[groups.indexOf(GROUP(s)) % PALETTE.length],
      'fill-opacity': flagged ? 1 : .8,
    }), `${s}: ${fmt(v)}` + (z == null ? '' : `  (robust z ${z.toFixed(2)})`)));
  });
  // The group median, drawn per group, so "is this sample like its batch" is visible rather
  // than arithmetic.
  for (const g of groups) {
    const idx = samples.map((s, i) => [s, i]).filter(([s]) => GROUP(s) === g);
    const med = median(idx.map(([s]) => ROWS[s][metric]));
    if (med === null) continue;
    svg.appendChild(el('line', {
      x1: pad.l + bw * idx[0][1], x2: pad.l + bw * (idx[idx.length - 1][1] + 1),
      y1: y(med), y2: y(med), stroke: 'currentColor', 'stroke-width': 1.5,
      'stroke-dasharray': '5 3', 'stroke-opacity': .7 }));
  }
  host.replaceChildren(svg);
}

// ── panel 3: the distributions ──────────────────────────────────────────────────────────────
let hidden = new Set();

function distKeys(scope) {
  const loci = new Set(), metrics = new Set();
  for (const s of DOC.samples)
    for (const [key, ms] of Object.entries((DOC.stats[s] || {})[scope] || {})) {
      loci.add(key.includes(':') ? key.split(':')[0] : '');
      for (const m of Object.keys(ms)) metrics.add(m);
    }
  return { loci: [...loci].sort(), metrics: [...metrics].sort() };
}

function drawDist() {
  const scope = $('#d-scope').value, locus = $('#d-locus').value, metric = $('#d-metric').value;
  const host = $('#dist');
  const series = [];
  for (const s of visibleSamples()) {
    const pts = [];
    for (const [key, ms] of Object.entries((DOC.stats[s] || {})[scope] || {})) {
      const [loc, bucket] = key.includes(':') ? [key.split(':')[0], key.slice(key.indexOf(':') + 1)]
                                              : ['', key];
      if (locus && loc !== locus) continue;
      if (!(metric in ms)) continue;
      pts.push([bucket, ms[metric]]);
    }
    if (pts.length) series.push([s, pts]);
  }
  if (!series.length) {
    host.replaceChildren(el('p', { class: 'empty', text: 'Nothing in this scope for these samples.' }));
    return;
  }
  // Buckets are numeric for a histogram and categorical for a gene or an isotype; one axis
  // handles both by falling back to the sorted label order.
  const labels = [...new Set(series.flatMap(([, pts]) => pts.map((p) => p[0])))];
  const numeric = labels.every((b) => b !== '' && isFinite(Number(b)));
  labels.sort(numeric ? (a, b) => Number(a) - Number(b) : (a, b) => a.localeCompare(b));
  const idx = Object.fromEntries(labels.map((b, i) => [b, i]));

  // Shown as a FRACTION of each sample's own total, or a deep sample hides every other curve.
  const shown = series.filter(([s]) => !hidden.has(s));
  const norm = $('#d-norm').checked;
  const scaled = shown.map(([s, pts]) => {
    const total = pts.reduce((a, p) => a + p[1], 0) || 1;
    return [s, pts.map(([b, v]) => [b, norm ? v / total : v])];
  });
  const yMax = Math.max(0, ...scaled.flatMap(([, pts]) => pts.map((p) => p[1]))) || 1;
  const w = 1040, h = 320, pad = { l: 70, r: 10, t: 10, b: 46 };
  const step = (w - pad.l - pad.r) / Math.max(labels.length, 1);
  const x = (b) => pad.l + step * (idx[b] + .5);
  const y = (v) => h - pad.b - (h - pad.b - pad.t) * (v / yMax);
  const every = Math.ceil(labels.length / 24);
  const svg = el('svg', { viewBox: `0 0 ${w} ${h}`, width: '100%', height: h });
  svg.appendChild(axes(w, h, pad, $('#d-scope').selectedOptions[0].textContent, yMax,
    labels.filter((_, i) => i % every === 0).map((b) => [x(b), b])));
  scaled.forEach(([s, pts], i) => {
    const colour = PALETTE[series.findIndex(([n]) => n === s) % PALETTE.length];
    const sorted = pts.slice().sort((a, b) => idx[a[0]] - idx[b[0]]);
    const d = sorted.map((p, k) => (k ? 'L' : 'M') + x(p[0]) + ' ' + y(p[1])).join(' ');
    svg.appendChild(el('path', { d: d, fill: 'none', stroke: colour, 'stroke-width': 1.8,
                                 'stroke-opacity': .85 }));
    for (const p of sorted)
      svg.appendChild(hover(el('circle', { cx: x(p[0]), cy: y(p[1]), r: 2.5, fill: colour }),
        `${s}  ${p[0]}: ${fmt(p[1])}`));
  });
  const legend = el('div', { class: 'legend' }, series.map(([s], i) =>
    el('span', { class: hidden.has(s) ? 'off' : '',
                 onclick: () => { hidden.has(s) ? hidden.delete(s) : hidden.add(s); drawDist(); } },
      [el('i', { style: 'background:' + PALETTE[i % PALETTE.length] }), s])));
  host.replaceChildren(svg, legend);
}

// ── panel 4: provenance ─────────────────────────────────────────────────────────────────────
function drawProvenance() {
  const keys = ['arda_version', 'mmseqs_version', 'reference.path', 'organism', 'map.organism',
    'map.threads', 'map.min_score', 'map.input', 'map.paired', 'map.shards', 'map.read_groups'];
  const rows = visibleSamples().map((s) => {
    const have = keys.filter((k) => ROWS[s][k] !== undefined);
    return el('div', {}, [el('h3', { text: s, style: 'font-size:.95rem;margin:.8rem 0 .2rem' }),
      el('dl', {}, have.flatMap((k) => [el('dt', { text: k }), el('dd', { text: fmt(ROWS[s][k]) })]))]);
  });
  $('#prov').replaceChildren(...rows);
}

// ── wiring ──────────────────────────────────────────────────────────────────────────────────
function fill(select, values, labels) {
  select.replaceChildren(...values.map((v, i) =>
    el('option', { value: v, text: (labels ? labels[i] : v) || '(all)' })));
}

function init() {
  const projects = [...new Set(DOC.samples.map((s) => LABEL(s).project))].filter(Boolean).sort();
  const batches = [...new Set(DOC.samples.map((s) => LABEL(s).batch))].filter(Boolean).sort();
  fill($('#f-project'), ['', ...projects]);
  fill($('#f-batch'), ['', ...batches]);
  $('#f-project').parentElement.style.display = projects.length ? '' : 'none';
  $('#f-batch').parentElement.style.display = batches.length ? '' : 'none';

  const metrics = [...new Set(DOC.samples.flatMap((s) =>
    Object.entries(ROWS[s]).filter(([, v]) => typeof v === 'number').map(([k]) => k)))].sort();
  fill($('#m-pick'), metrics);
  const preferred = ['map.mapped_fraction', 'reads', 'clonotypes'].find((m) => metrics.includes(m));
  if (preferred) $('#m-pick').value = preferred;

  const present = DISTS.filter(([scope]) =>
    DOC.samples.some((s) => Object.keys((DOC.stats[s] || {})[scope] || {}).length));
  fill($('#d-scope'), present.map((d) => d[0]), present.map((d) => d[1]));
  $('#dist-section').style.display = present.length ? '' : 'none';
  syncDist();

  for (const id of ['#f-project', '#f-batch'])
    $(id).addEventListener('change', () => { drawTable(); drawMetric(); drawDist(); drawProvenance(); });
  $('#f-text').addEventListener('input', () => { drawTable(); drawMetric(); drawDist(); drawProvenance(); });
  $('#f-all').addEventListener('change', (e) => { showAll = e.target.checked; drawTable(); });
  $('#m-pick').addEventListener('change', drawMetric);
  $('#d-scope').addEventListener('change', () => { syncDist(); drawDist(); });
  for (const id of ['#d-locus', '#d-metric', '#d-norm']) $(id).addEventListener('change', drawDist);

  drawTable(); drawMetric(); drawDist(); drawProvenance();
}

function syncDist() {
  const { loci, metrics } = distKeys($('#d-scope').value);
  fill($('#d-locus'), loci.length > 1 ? ['', ...loci] : loci);
  fill($('#d-metric'), metrics);
}

init();
"""


_PAGE = """<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>__CSS__</style>
</head>
<body>
<h1>__TITLE__</h1>
<p class="sub">Run quality control. Every figure is arda's own; this page holds its data and
needs nothing from the network.</p>

<div class="controls">
  <label>project <select id="f-project"></select></label>
  <label>batch <select id="f-batch"></select></label>
  <label>sample <input id="f-text" type="search" placeholder="filter…" size="14"></label>
  <label><input id="f-all" type="checkbox"> all metrics</label>
</div>

<h2>Samples</h2>
<div class="scroll" id="table"></div>
<p class="sub" id="table-note"></p>

<h2>One metric across the cohort</h2>
<div class="controls"><label>metric <select id="m-pick"></select></label></div>
<div id="metric"></div>

<section id="dist-section">
<h2>Distributions</h2>
<div class="controls">
  <label>scope <select id="d-scope"></select></label>
  <label>locus <select id="d-locus"></select></label>
  <label>count <select id="d-metric"></select></label>
  <label><input id="d-norm" type="checkbox" checked> as a fraction of each sample</label>
</div>
<div id="dist"></div>
</section>

<h2>Provenance</h2>
<div id="prov"></div>

<script type="application/json" id="qc-data">__PAYLOAD__</script>
<script>__JS__</script>
</body>
</html>
"""
