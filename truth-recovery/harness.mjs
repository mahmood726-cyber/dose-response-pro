// ============================================================
// harness.mjs -- Truth-recovery yardstick for dose-response-pro.
//
// Wires the app's OWN fitLinearModel (copied verbatim into engine.mjs) to the
// known-truth count-level dose-response DGP and measures how often the slope CI
// covers the TRUE log-RR-per-dose slope.
//
// The app builds the slope CI as a z-interval (beta[1] +/- 1.96*se[1]) with DL
// tau^2 by default. We evaluate two candidate calibration fixes from the same
// fitted model: a t critical value (df = nStudies-2) and REML tau^2.
//
// Truth-first: every number printed comes from seeded simulation here.
// Run:  node truth-recovery/harness.mjs --reps 300
// ============================================================

import { fitLinearModel, tCriticalValue, setUseREML } from './engine.mjs';
import { generate, makeRng, SCENARIOS } from './dgp-dr.mjs';

const BASE_SEED = 20260613;
const Z = 1.959963985;

// Returns slope estimate + four CIs from ONE fit per tau2-method.
function evaluate(points) {
  const out = {};
  for (const [tag, reml] of [['DL', false], ['REML', true]]) {
    setUseREML(reml);
    let r; try { r = fitLinearModel(points); } catch { continue; }
    const b = r.beta[1], se = r.se[1];
    const dfT = Math.max(r.nStudies - 2, 1);
    const tc = tCriticalValue(dfT, 0.975);
    out[`${tag}+z`] = { est: b, lo: b - Z * se, hi: b + Z * se };       // shipped path for DL
    out[`${tag}+t`] = { est: b, lo: b - tc * se, hi: b + tc * se };
  }
  setUseREML(false);
  return out;
}

const METHODS = ['DL+z', 'DL+t', 'REML+z', 'REML+t'];
const mean = (a) => a.reduce((x, y) => x + y, 0) / a.length;

export function runCell(slopeTrue, k, scenario, reps, rng) {
  const acc = {};
  for (const m of METHODS) acc[m] = { cov: 0, biasSum: 0, sq: 0, wSum: 0, n: 0 };
  for (let r = 0; r < reps; r++) {
    const { points } = generate(slopeTrue, k, scenario, rng);
    const ev = evaluate(points);
    for (const m of METHODS) {
      const o = ev[m];
      if (!o || !isFinite(o.est) || !isFinite(o.lo) || !isFinite(o.hi)) continue;
      const a = acc[m];
      a.n++;
      a.biasSum += o.est - slopeTrue;
      a.sq += (o.est - slopeTrue) ** 2;
      a.wSum += o.hi - o.lo;
      if (o.lo <= slopeTrue && slopeTrue <= o.hi) a.cov++;
    }
  }
  const res = {};
  for (const m of METHODS) {
    const a = acc[m];
    res[m] = {
      n: a.n,
      coverage: a.n ? +(a.cov / a.n).toFixed(4) : null,
      bias: a.n ? +(a.biasSum / a.n).toFixed(4) : null,
      rmse: a.n ? +Math.sqrt(a.sq / a.n).toFixed(4) : null,
      meanWidth: a.n ? +(a.wSum / a.n).toFixed(4) : null,
    };
  }
  return res;
}

export function runGrid({ reps = 300, ks = [4, 8, 15], slope = 0.30, scenarios = SCENARIOS } = {}) {
  const rng = makeRng(BASE_SEED);
  const grid = [];
  for (const scen of scenarios)
    for (const k of ks) grid.push({ scen, k, results: runCell(slope, k, scen, reps, rng) });
  return grid;
}

export function summarize(grid, filter = () => true) {
  const out = {};
  for (const m of METHODS) {
    const cov = [];
    for (const c of grid) if (filter(c) && c.results[m].coverage != null) cov.push(c.results[m].coverage);
    out[m] = { meanCoverage: +mean(cov).toFixed(4) };
  }
  return out;
}

const isMain = process.argv[1]?.endsWith('harness.mjs');
if (isMain) {
  const i = process.argv.indexOf('--reps');
  const reps = i >= 0 ? Number(process.argv[i + 1]) : 300;
  const t0 = Date.now();
  const grid = runGrid({ reps });
  const s = summarize(grid);
  console.log(`\n# Truth-recovery yardstick -- dose-response-pro`);
  console.log(`reps=${reps}/cell  true slope=0.30 logRR/dose  seed=${BASE_SEED}\n`);
  console.log('## Mean coverage of TRUE slope over all cells\n');
  console.log('method    meanCov');
  for (const m of METHODS) console.log(m.padEnd(9), String(s[m].meanCoverage).padStart(7));
  console.log('\n## Per-cell coverage (shipped DL+z vs candidates)\n');
  console.log('scenario   k    DL+z    DL+t   REML+z  REML+t');
  for (const c of grid) {
    console.log(c.scen.padEnd(10), String(c.k).padStart(2),
      String(c.results['DL+z'].coverage).padStart(7),
      String(c.results['DL+t'].coverage).padStart(7),
      String(c.results['REML+z'].coverage).padStart(7),
      String(c.results['REML+t'].coverage).padStart(7));
  }
  console.log(`\n(${(Date.now() - t0) / 1000}s)`);
}
