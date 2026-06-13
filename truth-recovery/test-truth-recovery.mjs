// node --test truth-recovery/test-truth-recovery.mjs
// Measured invariants for the dose-response-pro truth-recovery yardstick.
// Seeded; no hand-entered numbers.
import assert from 'node:assert/strict';
import { describe, it } from 'node:test';
import { fitLinearModel, setUseREML } from './engine.mjs';
import { generate, makeRng, SCENARIOS } from './dgp-dr.mjs';
import { runCell, runGrid, summarize } from './harness.mjs';

describe('dose-response DGP', () => {
  it('is reproducible for a fixed seed', () => {
    const a = generate(0.3, 8, 'tau_mod', makeRng(7));
    const b = generate(0.3, 8, 'tau_mod', makeRng(7));
    assert.deepEqual(a.points, b.points);
  });
  it('builds k studies each with a reference anchor + 3 dose rows', () => {
    const rng = makeRng(1);
    for (const scen of SCENARIOS) {
      const { points } = generate(0.3, 5, scen, rng);
      assert.equal(points.length, 5 * 4);
      const anchors = points.filter(p => p.dose === 0 && p.logRR === 0 && p.se === 0);
      assert.equal(anchors.length, 5);
    }
  });
});

describe('DL tau^2 bug (measured)', () => {
  it('estimateTau2DL collapses to exactly 0 even under strong heterogeneity', () => {
    // The denominator (sumTrV - df) subtracts an integer df from a sum of small
    // point VARIANCES, so it is almost always negative -> DL returns 0.
    setUseREML(false);
    const rng = makeRng(20260613);
    let allZero = true;
    for (let i = 0; i < 30; i++) {
      const { points } = generate(0.3, 8, 'tau_high', rng);
      if (fitLinearModel(points).tau2 !== 0) { allZero = false; break; }
    }
    assert.ok(allZero, 'DL tau^2 was not identically 0 across reps (bug may be fixed)');
  });
});

describe('Truth-recovery (measured)', () => {
  it('the shipped default (DL + z) UNDER-covers the true slope under heterogeneity', () => {
    const rng = makeRng(20260613);
    const hi = runCell(0.30, 15, 'tau_high', 200, rng);
    assert.ok(hi['DL+z'].coverage < 0.4,
      `DL+z coverage under high het = ${hi['DL+z'].coverage} (expected severe under-coverage)`);
  });

  it('REML produces tau^2>0 and recovers far more truth than the broken DL default', () => {
    const grid = runGrid({ reps: 150 });
    const s = summarize(grid);
    assert.ok(s['REML+z'].meanCoverage > s['DL+z'].meanCoverage + 0.3,
      `REML+z ${s['REML+z'].meanCoverage} not >> DL+z ${s['DL+z'].meanCoverage}`);
    // Honest: REML errs CONSERVATIVE here (it over-covers); it is a safe interim
    // default, not a perfectly calibrated one.
    assert.ok(s['REML+z'].meanCoverage >= 0.95,
      `REML+z mean coverage ${s['REML+z'].meanCoverage} (expected >= nominal / conservative)`);
  });
});
