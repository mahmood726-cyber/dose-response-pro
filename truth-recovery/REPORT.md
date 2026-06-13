# Truth-recovery yardstick — dose-response-pro

**Verdict: BUG found (measured) + a safe interim fix. The default DL τ² estimator
collapses to 0, so the shipped dose-response model is silently fixed-effect and
under-covers the true slope badly under heterogeneity.**

## The bug (`estimateTau2DL`)
```js
const denominator = sumTrV - df;           // sumTrV = Σ point variances ; df = (K-1)*p
if (!isFinite(Q) || denominator <= 0) return 0;
return Math.max(0, (Q - df) / denominator);
```
`sumTrV` is a **sum of point variances** (order ~0.1–1), while `df` is an **integer
count** `(K-1)·p` (order ~10). Subtracting a count from a sum of small variances
makes `denominator` almost always negative, so the guard returns **0**. Measured
example (K=8, high heterogeneity): `sumTrV=1.39`, `df=14`, `denominator=−12.6` → 0.
In 30/30 high-heterogeneity replications `estimateTau2DL` returned exactly 0.

Because the app ships `let useREML = false` (DL is the default), the dose-response
model is in practice **always fixed-effect** regardless of real between-study
heterogeneity.

## Method
- DGP (`dgp-dr.mjs`): count-level studies with a KNOWN true log-RR-per-dose slope
  (0.30) and between-study slope heterogeneity; reference + 3 dose arms,
  `a ~ Binomial(n, p0·exp(slope_i·dose))`. Empirical log-RR + SE passed in the
  app's own `{id,dose,logRR,se,cases}` point format. Seeded → reproducible.
- Engine (`engine.mjs`): `fitLinearModel` + GLS covariance + τ² (DL/REML) +
  `tCriticalValue` copied **verbatim** from `index.html` (lines 1332–3174).
- 300 reps/cell, true slope 0.30, `k∈{4,8,15}`, `τ² ∈ {0,0.01,0.04,0.10}` (slope
  heterogeneity). Coverage = how often the slope CI covers the true slope.

## Results — coverage of the TRUE slope

| method   | mean coverage |
|----------|--------------:|
| DL + z (shipped) | **0.475** |
| DL + t           | 0.819 |
| REML + z         | 0.985 |
| REML + t         | 1.000 |

Per-cell, the shipped DL+z collapses with heterogeneity (nominal 0.95):

| scenario  | k  | DL+z | REML+z |
|-----------|----|-----:|-------:|
| tau0      | 8  | 0.86 | 1.00   |
| tau_mod   | 15 | 0.18 | 0.997  |
| tau_high  | 15 | **0.06** | 0.993 |

## Findings (all measured)
1. **`estimateTau2DL` is broken** — its denominator is dimensionally inconsistent
   and returns 0 for essentially all realistic data. This is DGP-independent.
2. **Consequence:** the default dose-response slope CI under-covers the truth
   catastrophically under heterogeneity (down to 0.06), because it is effectively
   a fixed-effect interval.
3. **Safe interim fix:** set `useREML = true`. REML produces τ²>0 and restores
   coverage to ~0.985. **Honest caveat:** REML here errs *conservative* — it
   over-covers (REML+t ≈ 1.0). It is a safe interim default, not a perfectly
   calibrated one.
4. **Deeper note (honest):** the model's τ² is a single additive scalar on every
   contrast (a random-intercept structure), whereas dose-response heterogeneity is
   naturally a **random slope** (variance that scales with dose). A scalar τ² is a
   misspecified RE structure here; REML happens to over-inflate enough to stay
   conservative. The right long-term fix is a random-slope / one-stage model, not
   just a τ²-estimator swap. (Numeric τ² values are therefore NOT directly
   comparable to the injected slope-τ²; the *coverage* of the true slope is the
   scale-free ground truth this harness reports.)

## Recommendation
- **Short term:** flip the default to `useREML = true` (one line, measured-
  justified) so the app stops silently running fixed-effect. Or fix the DL
  denominator (correct multivariate-DL trace formula) and validate vs R
  `dosresmeta`.
- **Long term:** offer a random-slope heterogeneity structure for dose-response.

## What did NOT transfer
NPE/conformal/SBC/PartialID are estimator-of-μ machinery; this is a GLS dose-
response engine, so only the known-truth harness + a count-level dose-response DGP
transferred. No runtime dependency added.

## Reproduce
```
node truth-recovery/harness.mjs --reps 300
node --test truth-recovery/test-truth-recovery.mjs
```
