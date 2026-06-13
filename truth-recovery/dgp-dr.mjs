// ============================================================
// dgp-dr.mjs -- Known-truth DGP for DOSE-RESPONSE meta-analysis.
//
// dose-response-pro fits a GLS dose-response model (within-study Greenland-
// Longnecker covariance + between-study tau^2). It is dosresmeta/R-validated on
// points; that does not show whether the slope CI RECOVERS the true dose-response
// slope under between-study heterogeneity. This DGP supplies that test by
// generating count-level studies with a KNOWN true log-RR-per-dose slope.
//
// Generative model: study i has slope_i = slope + sqrt(tau2)*Z (heterogeneity on
// the dose-response slope). Reference arm: n0 subjects, risk p0. Dose arm j at
// dose d_j: risk_j = p0 * exp(slope_i * d_j) (log-RR linear in dose). Counts
// a ~ Binomial(n, risk). Empirical log-RR and its SE (with 0.5 cc on zero cells)
// are passed to the app in its own {id,dose,logRR,se,cases} point format, with a
// deterministic reference anchor row (logRR=0, se=0) the app removes.
//
// ESTIMAND: the true slope (log-RR per unit dose). Seeded -> reproducible.
// ============================================================

export const SCENARIOS = ['tau0', 'tau_low', 'tau_mod', 'tau_high'];
const TAU2 = { tau0: 0, tau_low: 0.01, tau_mod: 0.04, tau_high: 0.10 };

export function makeRng(seed) {
  let a = seed >>> 0;
  return function () {
    a |= 0; a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}
function randn(rng) {
  let u1 = rng(), u2 = rng();
  if (u1 < 1e-12) u1 = 1e-12;
  return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
}
function rbinom(rng, n, p) { let x = 0; for (let i = 0; i < n; i++) if (rng() < p) x++; return x; }
function drawN(rng, lo, hi) {
  const a = Math.log(lo), b = Math.log(hi);
  return Math.max(20, Math.round(Math.exp(a + (b - a) * rng())));
}

const DOSES = [0, 1, 2, 3];   // reference + 3 dose levels per study

export function generate(slopeTrue, k, scenario, rng,
                         { p0 = 0.12, nLo = 80, nHi = 600 } = {}) {
  const tau2 = TAU2[scenario] ?? 0;
  const sd = Math.sqrt(tau2);
  const points = [];
  for (let s = 0; s < k; s++) {
    const slope_i = slopeTrue + sd * randn(rng);
    const n0 = drawN(rng, nLo, nHi);
    let a0 = rbinom(rng, n0, p0);
    const a0c = a0 === 0 ? 0.5 : a0;
    points.push({ id: `st${s}`, dose: 0, logRR: 0, se: 0, cases: a0 });   // anchor
    for (let j = 1; j < DOSES.length; j++) {
      const d = DOSES[j];
      const risk = Math.min(0.95, Math.max(1e-4, p0 * Math.exp(slope_i * d)));
      const nj = drawN(rng, nLo, nHi);
      let aj = rbinom(rng, nj, risk);
      const ajc = aj === 0 ? 0.5 : aj;
      const logRR = Math.log((ajc / nj) / (a0c / n0));
      const se = Math.sqrt(1 / ajc - 1 / nj + 1 / a0c - 1 / n0);
      points.push({ id: `st${s}`, dose: d, logRR, se: isFinite(se) && se > 0 ? se : 0.2, cases: aj });
    }
  }
  return { points, slopeTrue, tau2, info: { k } };
}
