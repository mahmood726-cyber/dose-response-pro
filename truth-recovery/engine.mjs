// engine.mjs -- pure dose-response GLS core EXTRACTED VERBATIM from index.html
// (dose-response-pro), lines 1332-3174 (matrix utils, GLS covariance, model
// fitting, tau2 DL/REML, Q, fit stats, chiSqCDF, normCDFinv, tCriticalValue).
// DOM globals it needs are declared here; useREML is the app's own (line 2900),
// re-exposed via setUseREML so the harness can test DL vs REML.

const NUMERICAL_TOLERANCE = 1e-10;
const DETERMINANT_THRESHOLD = 1e-10;
const RIDGE_PENALTY = 1e-10;
let currentCI = 95;
let currentMainModel = 'linear';

function invertMatrix(V, n) {
  // Create augmented matrix [V|I]
  const aug = new Array(n * 2 * n).fill(0);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      aug[i * 2 * n + j] = V[i * n + j];
    }
    aug[i * 2 * n + n + i] = 1;  // Identity matrix on right
  }

  // Gaussian elimination with partial pivoting
  for (let col = 0; col < n; col++) {
    // Find pivot
    let maxRow = col;
    let maxVal = Math.abs(aug[col * 2 * n + col]);
    for (let row = col + 1; row < n; row++) {
      if (Math.abs(aug[row * 2 * n + col]) > maxVal) {
        maxVal = Math.abs(aug[row * 2 * n + col]);
        maxRow = row;
      }
    }

    // Swap rows
    if (maxRow !== col) {
      for (let j = 0; j < 2 * n; j++) {
        [aug[col * 2 * n + j], aug[maxRow * 2 * n + j]] =
          [aug[maxRow * 2 * n + j], aug[col * 2 * n + j]];
      }
    }

    // Check for singular matrix
    if (Math.abs(aug[col * 2 * n + col]) < DETERMINANT_THRESHOLD) {
      console.warn('Matrix is near-singular, adding ridge penalty');
      aug[col * 2 * n + col] += RIDGE_PENALTY;
    }

    // Scale pivot row
    const pivot = aug[col * 2 * n + col];
    for (let j = 0; j < 2 * n; j++) {
      aug[col * 2 * n + j] /= pivot;
    }

    // Eliminate column
    for (let row = 0; row < n; row++) {
      if (row !== col) {
        const factor = aug[row * 2 * n + col];
        for (let j = 0; j < 2 * n; j++) {
          aug[row * 2 * n + j] -= factor * aug[col * 2 * n + j];
        }
      }
    }
  }

  // Extract inverse (right half of augmented matrix)
  const inv = new Array(n * n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      inv[i * n + j] = aug[i * 2 * n + n + j];
    }
  }

  return inv;
}

/**
 * Invert a block-diagonal matrix by inverting each block separately
 * This is the CORRECT implementation for GLS with within-study correlation
 * @param {Array} V_blocks - Array of covariance blocks
 * @returns {Array} - Array of inverted blocks
 */
function invertBlockDiagonal(V_blocks) {
  return V_blocks.map(block => {
    const n = block.n;
    const V = new Array(n * n);

    // Extract block as 2D array for inversion
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        V[i * n + j] = block.V[i * n + j];
      }
    }

    // Invert the full covariance matrix (NOT just diagonal!)
    const V_inv = invertMatrix(V, n);

    return { V_inv, n };
  });
}

/**
 * Build GLS covariance matrix for dose-response meta-analysis
 *
 * Reference: Greenland S, Longnecker MP. Methods for trend estimation from
 * summarized dose-response data. Am J Epidemiol. 1992;135(11):1301-1309.
 *
 * Also: Orsini N, et al. Meta-analysis for linear and nonlinear dose-response
 * relations. Am J Epidemiol. 2012;175(1):66-73.
 *
 * IMPLEMENTATION NOTE: The full Greenland-Longnecker covariance formula
 * requires raw case counts from the reference category:
 *
 *   Cov(log RR_i, log RR_j) = 1/n_0 + 1/N_0
 *
 * where n_0 = cases in reference, N_0 = total in reference.
 *
 * When reference case counts are available (dose = 0), we use the
 * Greenland-Longnecker shared-reference covariance (Cov(i,j) = 1/n0).
 * When only SEs are available, we use a Hamling-style correlation
 * approximation with rho (default 0.5) and fall back to diagonal if
 * the matrix is not positive semi-definite.
 *
 * @param {Array} studyPoints - Array of non-reference study data points with se property
 * @param {number} rho - Correlation fallback when shared-reference counts are unavailable
 * @param {Object|null} referencePoint - Optional reference-category row used for shared-reference covariance
 * @returns {Array} - Covariance matrix as flat array (n x n)
 */
function buildGLSCovariance(studyPoints, rho = 0.5, referencePoint = null) {
  const n = studyPoints.length;
  if (n === 0) return [];

  rho = Math.max(0.1, Math.min(0.9, rho));
  const V = new Array(n * n).fill(0);

  const validSEs = studyPoints.filter(p => typeof p.se === 'number' && p.se > 0).map(p => p.se);
  const medianSE = validSEs.length > 0
    ? validSEs.sort((a, b) => a - b)[Math.floor(validSEs.length / 2)]
    : 0.1;

  const seValues = studyPoints.map((pt, i) => {
    if (typeof pt.se !== 'number' || pt.se <= 0) {
      console.warn(`Invalid SE for point ${i}, using median SE=${medianSE.toFixed(4)}`);
      return medianSE;
    }
    return pt.se;
  });

  const baselineIndex = studyPoints.findIndex(p => p.dose === 0);
  const baselinePointInBlock = baselineIndex >= 0 ? studyPoints[baselineIndex] : null;
  const effectiveReference = referencePoint || baselinePointInBlock;
  const hasExactCounts = Number.isFinite(effectiveReference?.cases) && effectiveReference.cases > 0;

  if (hasExactCounts) {
    const n0 = effectiveReference.cases;
    const refCov = 1 / Math.max(n0, 1);
    const variances = seValues.map(se => se * se);

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (i === j) {
          V[i * n + j] = variances[i];
        } else {
          if (i === baselineIndex || j === baselineIndex) {
            // The reference category is a fixed contrast (logRR=0); do not
            // induce correlation between reference and non-reference rows.
            V[i * n + j] = 0;
            continue;
          }
          const covBound = Math.sqrt(variances[i] * variances[j]) * 0.9;
          V[i * n + j] = Math.min(refCov, covBound);
        }
      }
    }

    if (!isPositiveSemiDefinite(V, n)) {
      console.warn('Covariance matrix not positive semi-definite, using diagonal');
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          V[i * n + j] = (i === j) ? variances[i] : 0;
        }
      }
    }
  } else {
    const variances = seValues.map(se => se * se);

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        if (i === j) {
          V[i * n + j] = variances[i];
        } else {
          const covBound = Math.min(variances[i], variances[j]);
          V[i * n + j] = Math.min(rho * seValues[i] * seValues[j], covBound * 0.9);
        }
      }
    }

    if (!isPositiveSemiDefinite(V, n)) {
      console.warn('Covariance matrix not positive semi-definite, using diagonal');
      for (let i = 0; i < n; i++) {
        for (let j = 0; j < n; j++) {
          V[i * n + j] = (i === j) ? variances[i] : 0;
        }
      }
    }
  }

  return V;
}

function detectStudyReferencePoint(studyPoints) {
  // Treat only deterministic anchors (logRR=0 with non-positive/unknown SE)
  // as removable reference rows. Valid dose-0 rows with positive SE are kept.
  const fixedAnchor = studyPoints.find(p =>
    Number.isFinite(p.logRR) &&
    p.logRR === 0 &&
    (!Number.isFinite(p.se) || p.se <= 0)
  );

  return fixedAnchor || null;
}

function prepareStudyBlocksForModel(points) {
  const studyMap = new Map();
  for (const pt of points) {
    if (!studyMap.has(pt.id)) studyMap.set(pt.id, []);
    studyMap.get(pt.id).push(pt);
  }

  const blocks = [];
  const covarianceSummary = {
    sharedReference: 0,
    approximate: 0,
    totalBlocks: 0,
    referenceRowsExcluded: 0,
    droppedStudies: 0
  };

  for (const [studyId, studyPoints] of studyMap.entries()) {
    const referencePoint = detectStudyReferencePoint(studyPoints);
    const modelPoints = referencePoint
      ? studyPoints.filter(pt => pt !== referencePoint)
      : studyPoints.slice();

    if (referencePoint) {
      covarianceSummary.referenceRowsExcluded += 1;
    }

    if (modelPoints.length === 0) {
      covarianceSummary.droppedStudies += 1;
      continue;
    }

    const hasSharedReference = (Number.isFinite(referencePoint?.cases) && referencePoint.cases > 0) ||
      modelPoints.some(pt => pt.dose === 0 && Number.isFinite(pt.cases) && pt.cases > 0);
    if (hasSharedReference) {
      covarianceSummary.sharedReference += 1;
    } else {
      covarianceSummary.approximate += 1;
    }

    const V = buildGLSCovariance(modelPoints, 0.5, referencePoint);
    blocks.push({
      studyId,
      modelPoints,
      originalPoints: studyPoints.slice(),
      referencePoint,
      V
    });
  }

  covarianceSummary.totalBlocks = blocks.length;
  return { blocks, covarianceSummary };
}

function isPositiveSemiDefinite(V, n) {
  for (let i = 0; i < n; i++) {
    if (V[i * n + i] <= 0) return false;
  }

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      if (i !== j) {
        const bound = Math.sqrt(V[i * n + i] * V[j * n + j]);
        if (Math.abs(V[i * n + j]) > bound + 1e-10) return false;
      }
    }
  }

  if (n === 2) {
    const det = V[0] * V[3] - V[1] * V[2];
    return det >= -1e-10;
  }

  return true;
}

// ================================================================
// ADDITIONAL MODEL FAMILIES
// ================================================================

/**
 * Compute fit statistics shared across model families.
 */
function computeModelFitStats(result, p, nPoints) {
  if (!result || !Number.isFinite(result.WRSS) || !Number.isFinite(result.detV) || !Number.isFinite(nPoints)) {
    return { logLik: NaN, AIC: NaN, BIC: NaN };
  }

  const logLik = -0.5 * (nPoints * Math.log(2 * Math.PI) + result.WRSS + result.detV);
  const AIC = 2 * p - 2 * logLik;
  const BIC = Math.log(Math.max(nPoints, 1)) * p - 2 * logLik;

  return { logLik, AIC, BIC };
}

/**
 * Build restricted cubic spline basis matrix
 *
 * Reference: Durrleman, S., & Simon, R. (1989). Flexible regression models
 * with cubic splines. Statistics in medicine, 8(5), 551-561.
 *
 * @param {Array} doses - Array of dose values
 * @param {number} numKnots - Number of knots (default: 4, must be ≥ 3)
 * @returns {Object} - {basis, knots} where basis is the design matrix
 */
function buildSplineBasis(doses, numKnots = 4) {
  const n = doses.length;
  const knots = [];

  // Place knots at specified percentiles of dose distribution
  const sortedDoses = [...doses].sort((a, b) => a - b);
  for (let k = 0; k < numKnots; k++) {
    const idx = Math.floor((k + 1) * (n - 1) / (numKnots + 1));
    knots.push(sortedDoses[Math.min(idx, n - 1)]);
  }

  // Build basis matrix: [intercept, dose, spline1, spline2, ...]
  const basis = [];

  for (const dose of doses) {
    const row = [1, dose];  // Intercept and linear term

    // Add spline terms for interior knots (knots 1 to K-2)
    for (let k = 1; k < numKnots - 1; k++) {
      const knot = knots[k];
      const term = computeSplineTerm(dose, knot, knots[0], knots[numKnots - 1]);
      row.push(term);
    }

    basis.push(row);
  }

  return { basis, knots };
}

/**
 * Compute a single restricted cubic spline term
 *
 * Formula: S_k(d) = (d - t_k)^3_+ - (d - t_max)^3_+ * (t_max - t_k)/(t_max - t_min)
 *                   + (d - t_min)^3_+ * (t_k - t_min)/(t_max - t_min)
 *
 * where (x)_+ = max(0, x)
 */
function computeSplineTerm(dose, knot, t_min, t_max) {
  const x = dose - knot;
  const x_min = dose - t_min;
  const x_max = dose - t_max;

  // Positive part function
  const pos = (val) => Math.max(0, val);
  const cube = (val) => val * val * val;

  const numerator = cube(pos(x)) -
                   cube(pos(x_max)) * (t_max - knot) / (t_max - t_min) +
                   cube(pos(x_min)) * (knot - t_min) / (t_max - t_min);

  const denominator = 6 * (t_max - t_min);

  return denominator !== 0 ? numerator / denominator : 0;
}

/**
 * Fit restricted cubic spline dose-response model
 *
 * @param {Array} points - Study data points
 * @param {number} tau2Override - Optional tau² override
 * @param {number} numKnots - Number of spline knots (default: 4)
 * @returns {Object} - Fitted model results
 */
function fitSplineModel(points, tau2Override = null, numKnots = 4) {
  if (numKnots < 3 || numKnots > 7) {
    throw new Error('Number of knots must be between 3 and 7');
  }

  const prepared = prepareStudyBlocksForModel(points);
  const studyBlocks = prepared.blocks;
  const covarianceSummary = prepared.covarianceSummary;

  if (studyBlocks.length < 1) {
    throw new Error('Need at least 1 study with non-reference dose contrasts');
  }

  const allDoses = [];
  for (const block of studyBlocks) {
    for (const pt of block.modelPoints) {
      allDoses.push(pt.dose);
    }
  }

  if (allDoses.length < 3) {
    throw new Error('Need at least 3 non-reference dose points for spline model');
  }

  const uniqueDoses = [...new Set(allDoses)];
  if (uniqueDoses.length < 3) {
    throw new Error('Need at least 3 unique non-reference dose levels for spline model');
  }

  if (uniqueDoses.length < numKnots) {
    numKnots = Math.min(uniqueDoses.length, 4);
  }

  const { basis: splineBasis, knots } = buildSplineBasis(allDoses, numKnots);
  const p = splineBasis[0].length;

  const X = [];
  const y = [];
  const V_blocks = [];

  let basisIdx = 0;
  for (const block of studyBlocks) {
    for (const pt of block.modelPoints) {
      X.push(splineBasis[basisIdx]);
      y.push(pt.logRR);
      basisIdx += 1;
    }
    V_blocks.push({ studyId: block.studyId, V: block.V, n: block.modelPoints.length });
  }

  const nStudies = V_blocks.length;
  const nPoints = y.length;

  if (nPoints <= p) {
    throw new Error('Insufficient non-reference points for spline model estimation');
  }

  let tau2 = tau2Override;
  if (tau2 === null) {
    const tauResult = estimateTau2(X, y, V_blocks);
    tau2 = tauResult.tau2;
  }

  const result = solveGLSWithTau2Spline(X, y, V_blocks, tau2, p);

  const Q = computeQStatSpline(y, X, result.beta, V_blocks, tau2, p);
  const df = Math.max((nStudies - 1) * p, 0);
  const I2 = (df > 0 && Q > 0) ? Math.max(0, 100 * (Q - df) / Q) : 0;
  const Qp = df > 0 ? 1 - chiSqCDF(Q, df) : NaN;
  const fitStats = computeModelFitStats(result, p, nPoints);

  return {
    ...result,
    tau2,
    Q,
    df,
    I2,
    Qp,
    nStudies,
    nPoints,
    logLik: fitStats.logLik,
    AIC: fitStats.AIC,
    BIC: fitStats.BIC,
    knots,
    modelType: 'spline',
    numKnots,
    covarianceSummary
  };
}


/**
 * Fit linear dose-response model (intercept + linear term)
 *
 * @param {Array} points - Study data points
 * @param {number} tau2Override - Optional tau² override
 * @returns {Object} - Fitted model results
 */
function fitLinearModel(points, tau2Override = null) {
  const prepared = prepareStudyBlocksForModel(points);
  const studyBlocks = prepared.blocks;
  const covarianceSummary = prepared.covarianceSummary;

  if (studyBlocks.length < 1) {
    throw new Error('Need at least 1 study with non-reference dose contrasts');
  }

  const X = [];
  const y = [];
  const V_blocks = [];

  for (const block of studyBlocks) {
    for (const pt of block.modelPoints) {
      X.push([1, pt.dose]);
      y.push(pt.logRR);
    }
    V_blocks.push({ studyId: block.studyId, V: block.V, n: block.modelPoints.length });
  }

  const p = 2;
  const nStudies = V_blocks.length;
  const nPoints = y.length;

  if (nPoints <= p) {
    throw new Error('Insufficient non-reference points for linear model estimation');
  }

  let tau2 = tau2Override;
  if (tau2 === null) {
    const tauResult = estimateTau2(X, y, V_blocks);
    tau2 = tauResult.tau2;
  }

  const result = solveGLSLinear(X, y, V_blocks, tau2);

  const Q = computeQStatLinear(y, X, result.beta, V_blocks, tau2);
  const df = Math.max((nStudies - 1) * p, 0);
  const I2 = (df > 0 && Q > 0) ? Math.max(0, 100 * (Q - df) / Q) : 0;
  const Qp = df > 0 ? 1 - chiSqCDF(Q, df) : NaN;
  const fitStats = computeModelFitStats(result, p, nPoints);

  return {
    ...result,
    tau2,
    Q,
    df,
    I2,
    Qp,
    nStudies,
    nPoints,
    logLik: fitStats.logLik,
    AIC: fitStats.AIC,
    BIC: fitStats.BIC,
    modelType: 'linear',
    covarianceSummary
  };
}


/**
 * Solve GLS for linear model (2 parameters)
 */
function solveGLSLinear(X, y, V_blocks, tau2) {
  const p = 2;

  const V_inv_blocks = V_blocks.map(block => {
    const blockSize = block.n;
    const V_with_tau2 = new Array(blockSize * blockSize);

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        V_with_tau2[i * blockSize + j] = block.V[i * blockSize + j] + (i === j ? tau2 : 0);
      }
    }

    const V_inv = invertMatrix(V_with_tau2, blockSize);
    return { V_inv, n: blockSize };
  });

  const XtVinvX = new Array(p * p).fill(0);
  const XtVinvY = new Array(p).fill(0);

  let row = 0;
  for (const blockInv of V_inv_blocks) {
    const blockSize = blockInv.n;
    const V_inv = blockInv.V_inv;

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        const w_ij = V_inv[i * blockSize + j];

        for (let k = 0; k < p; k++) {
          XtVinvY[k] += w_ij * X[row + i][k] * y[row + j];
          for (let l = 0; l < p; l++) {
            XtVinvX[k * p + l] += w_ij * X[row + i][k] * X[row + j][l];
          }
        }
      }
    }
    row += blockSize;
  }

  for (let i = 0; i < p; i++) {
    XtVinvX[i * p + i] += RIDGE_PENALTY;
  }

  const det = XtVinvX[0] * XtVinvX[3] - XtVinvX[1] * XtVinvX[2];
  const beta = det !== 0 ? [
    (XtVinvY[0] * XtVinvX[3] - XtVinvY[1] * XtVinvX[1]) / det,
    (XtVinvX[0] * XtVinvY[1] - XtVinvX[2] * XtVinvY[0]) / det
  ] : [0, 0];

  const varMatrix = det !== 0 ? [
    XtVinvX[3] / det, -XtVinvX[1] / det,
    -XtVinvX[2] / det, XtVinvX[0] / det
  ] : [0, 0, 0, 0];

  const se = [
    Math.sqrt(Math.max(varMatrix[0], 0)),
    Math.sqrt(Math.max(varMatrix[3], 0))
  ];

  let WRSS = 0;
  row = 0;
  for (const blockInv of V_inv_blocks) {
    const blockSize = blockInv.n;
    const V_inv = blockInv.V_inv;

    for (let i = 0; i < blockSize; i++) {
      const pred = beta[0] + beta[1] * X[row + i][1];
      const resid_i = y[row + i] - pred;
      for (let j = 0; j < blockSize; j++) {
        const pred_j = beta[0] + beta[1] * X[row + j][1];
        const resid_j = y[row + j] - pred_j;
        WRSS += resid_i * V_inv[i * blockSize + j] * resid_j;
      }
    }
    row += blockSize;
  }

  let detV = 0;
  for (const block of V_blocks) {
    const blockSize = block.n;
    for (let i = 0; i < blockSize; i++) {
      detV += Math.log(Math.max(block.V[i * blockSize + i] + tau2, NUMERICAL_TOLERANCE));
    }
  }

  return { beta, se, WRSS, detV, varMatrix };
}

/**
 * Estimate tau2 for linear model using DL
 */
function estimateTau2Linear(X, y, V_blocks) {
  const K = V_blocks.length;
  const p = 2;

  const feResult = solveGLSLinear(X, y, V_blocks, 0);
  const Q = computeQStatLinear(y, X, feResult.beta, V_blocks, 0);
  const df = (K - 1) * p;

  let sumTrV = 0;
  for (const block of V_blocks) {
    const blockSize = block.n;
    for (let i = 0; i < blockSize; i++) {
      sumTrV += block.V[i * blockSize + i];
    }
  }

  const denom = sumTrV - df;
  if (!Number.isFinite(Q) || !Number.isFinite(denom) || denom <= 0) {
    return 0;
  }

  return Math.max(0, (Q - df) / denom);
}

/**
 * Compute Q statistic for linear model
 */
function computeQStatLinear(y, X, beta, V_blocks, tau2) {
  let Q = 0;
  let row = 0;

  for (const block of V_blocks) {
    const blockSize = block.n;
    const V = new Array(blockSize * blockSize);

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        V[i * blockSize + j] = block.V[i * blockSize + j] + (i === j ? tau2 : 0);
      }
    }

    const V_inv = invertMatrix(V, blockSize);

    for (let i = 0; i < blockSize; i++) {
      const pred = beta[0] + beta[1] * X[row + i][1];
      const resid_i = y[row + i] - pred;
      for (let j = 0; j < blockSize; j++) {
        const pred_j = beta[0] + beta[1] * X[row + j][1];
        const resid_j = y[row + j] - pred_j;
        Q += resid_i * V_inv[i * blockSize + j] * resid_j;
      }
    }
    row += blockSize;
  }

  return Q;
}

/**
 * Solve GLS with spline basis (general p x p case)
 */
function solveGLSWithTau2Spline(X, y, V_blocks, tau2, p) {
  const n = y.length;

  // Invert each covariance block
  const V_inv_blocks = V_blocks.map(block => {
    const blockSize = block.n;
    const V_with_tau2 = new Array(blockSize * blockSize);

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        V_with_tau2[i * blockSize + j] = block.V[i * blockSize + j] + (i === j ? tau2 : 0);
      }
    }

    const V_inv = invertMatrix(V_with_tau2, blockSize);
    return { V_inv, n: blockSize };
  });

  // Compute X'V^(-1)X and X'V^(-1)y
  const XtVinvX = new Array(p * p).fill(0);
  const XtVinvY = new Array(p).fill(0);

  let row = 0;
  for (const blockInv of V_inv_blocks) {
    const blockSize = blockInv.n;
    const V_inv = blockInv.V_inv;

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        const w_ij = V_inv[i * blockSize + j];

        for (let k = 0; k < p; k++) {
          XtVinvY[k] += w_ij * X[row + i][k] * y[row + j];
          for (let l = 0; l < p; l++) {
            XtVinvX[k * p + l] += w_ij * X[row + i][k] * X[row + j][l];
          }
        }
      }
    }
    row += blockSize;
  }

  // Add ridge penalty
  for (let i = 0; i < p; i++) {
    XtVinvX[i * p + i] += RIDGE_PENALTY;
  }

  // Solve for beta using general matrix solver
  const beta = solveGeneralMatrix(XtVinvX, XtVinvY, p);

  // Compute variance-covariance
  const XtVinvX_inv = invertMatrix(XtVinvX, p);
  const se = [];
  for (let i = 0; i < p; i++) {
    se.push(Math.sqrt(Math.max(XtVinvX_inv[i * p + i], 0)));
  }

  // Compute WRSS
  let WRSS = 0;
  row = 0;
  for (const blockInv of V_inv_blocks) {
    const blockSize = blockInv.n;
    const V_inv = blockInv.V_inv;

    // Pre-compute all residuals for this block
    const residuals = [];
    for (let i = 0; i < blockSize; i++) {
      let pred_i = 0;
      for (let k = 0; k < p; k++) {
        pred_i += beta[k] * X[row + i][k];
      }
      residuals.push(y[row + i] - pred_i);
    }

    // Compute quadratic form
    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        WRSS += residuals[i] * V_inv[i * blockSize + j] * residuals[j];
      }
    }
    row += blockSize;
  }

  // Compute log determinant
  let detV = 0;
  for (const block of V_blocks) {
    const blockSize = block.n;
    for (let i = 0; i < blockSize; i++) {
      detV += Math.log(Math.max(block.V[i * blockSize + i] + tau2, NUMERICAL_TOLERANCE));
    }
  }

  return { beta, se, WRSS, detV, varMatrix: XtVinvX_inv };
}

/**
 * Compute Q statistic for spline model
 */
function computeQStatSpline(y, X, beta, V_blocks, tau2, p) {
  let Q = 0;
  let row = 0;

  for (const block of V_blocks) {
    const blockSize = block.n;
    const V_with_tau2 = new Array(blockSize * blockSize);

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        V_with_tau2[i * blockSize + j] = block.V[i * blockSize + j] + (i === j ? tau2 : 0);
      }
    }

    const V_inv = invertMatrix(V_with_tau2, blockSize);

    // Pre-compute all residuals for this block
    const residuals = [];
    for (let i = 0; i < blockSize; i++) {
      let pred_i = 0;
      for (let k = 0; k < p; k++) {
        pred_i += beta[k] * X[row + i][k];
      }
      residuals.push(y[row + i] - pred_i);
    }

    // Compute quadratic form
    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        Q += residuals[i] * V_inv[i * blockSize + j] * residuals[j];
      }
    }
    row += blockSize;
  }

  return Q;
}

/**
 * Solve general linear system Ax = b for any size matrix
 * Uses Gaussian elimination with partial pivoting
 */
function solveGeneralMatrix(A, b, p) {
  // Create augmented matrix [A|b]
  const aug = new Array(p * (p + 1));
  for (let i = 0; i < p; i++) {
    for (let j = 0; j < p; j++) {
      aug[i * (p + 1) + j] = A[i * p + j];
    }
    aug[i * (p + 1) + p] = b[i];
  }

  // Gaussian elimination with partial pivoting
  for (let col = 0; col < p; col++) {
    // Find pivot
    let maxRow = col;
    let maxVal = Math.abs(aug[col * (p + 1) + col]);
    for (let row = col + 1; row < p; row++) {
      if (Math.abs(aug[row * (p + 1) + col]) > maxVal) {
        maxVal = Math.abs(aug[row * (p + 1) + col]);
        maxRow = row;
      }
    }

    // Swap rows
    if (maxRow !== col) {
      for (let j = col; j <= p; j++) {
        [aug[col * (p + 1) + j], aug[maxRow * (p + 1) + j]] =
          [aug[maxRow * (p + 1) + j], aug[col * (p + 1) + j]];
      }
    }

    // Check for singular matrix
    if (Math.abs(aug[col * (p + 1) + col]) < DETERMINANT_THRESHOLD) {
      aug[col * (p + 1) + col] += RIDGE_PENALTY;
    }

    // Eliminate column
    for (let row = col + 1; row < p; row++) {
      const factor = aug[row * (p + 1) + col] / aug[col * (p + 1) + col];
      for (let j = col; j <= p; j++) {
        aug[row * (p + 1) + j] -= factor * aug[col * (p + 1) + j];
      }
    }
  }

  // Back substitution
  const x = new Array(p).fill(0);
  for (let i = p - 1; i >= 0; i--) {
    x[i] = aug[i * (p + 1) + p];
    for (let j = i + 1; j < p; j++) {
      x[i] -= aug[i * (p + 1) + j] * x[j];
    }
    x[i] /= aug[i * (p + 1) + i];
  }

  return x;
}

/**
 * Fit exponential dose-response model
 * Model: log(RR) = β₀ + β₁ × (1 - exp(-α × dose))
 *
 * Reference: Rota, M., et al. (2010). Random-effects dose-response model
 * for pooling non-linear dose-response data from epidemiological studies.
 *
 * @param {Array} points - Study data points
 * @param {number} tau2Override - Optional tau² override
 * @param {number} alpha - Saturation parameter (default: estimated from data)
 * @returns {Object} - Fitted model results
 */
function fitExponentialModel(points, tau2Override = null, alpha = null) {
  const prepared = prepareStudyBlocksForModel(points);
  const studyBlocks = prepared.blocks;
  const covarianceSummary = prepared.covarianceSummary;

  if (studyBlocks.length < 1) {
    throw new Error('Need at least 1 study with non-reference dose contrasts');
  }

  // Estimate alpha from data if not provided
  // Alpha ≈ 1 / (mean of non-zero doses)
  if (alpha === null) {
    const nonZeroDoses = studyBlocks.flatMap(block =>
      block.modelPoints
        .map(p => p.dose)
        .filter(dose => Number.isFinite(dose) && dose > 0)
    );
    if (nonZeroDoses.length > 0) {
      const meanDose = nonZeroDoses.reduce((a, b) => a + b, 0) / nonZeroDoses.length;
      alpha = 1 / Math.max(meanDose, 0.1);
    } else {
      alpha = 0.1;  // Default value
    }
  }

  // Build design matrix X: [intercept, saturation term]
  const X = [];
  const y = [];
  const V_blocks = [];

  for (const block of studyBlocks) {
    for (const pt of block.modelPoints) {
      // Exponential model: 1 - exp(-α × dose)
      const saturation = 1 - Math.exp(-alpha * pt.dose);
      X.push([1, saturation]);
      y.push(pt.logRR);
    }

    V_blocks.push({ studyId: block.studyId, V: block.V, n: block.modelPoints.length });
  }

  const p = 2;
  const nStudies = V_blocks.length;
  const nPoints = y.length;

  if (nPoints <= p) {
    throw new Error('Insufficient non-reference points for exponential model estimation');
  }

  // Estimate tau²
  let tau2 = tau2Override;
  if (tau2 === null) {
    const tauResult = estimateTau2(X, y, V_blocks);
    tau2 = tauResult.tau2;
  }

  // Solve GLS
  const result = solveGLSWithTau2(X, y, V_blocks, tau2);

  // Compute Q and I²
  const Q = computeQStat(y, X, result.beta, V_blocks, tau2);
  const df = Math.max((nStudies - 1) * p, 0);
  const I2 = (df > 0 && Q > 0) ? Math.max(0, 100 * (Q - df) / Q) : 0;
  const Qp = df > 0 ? 1 - chiSqCDF(Q, df) : NaN;
  const fitStats = computeModelFitStats(result, p, nPoints);

  return {
    ...result,
    tau2,
    Q,
    df,
    I2,
    Qp,
    nStudies,
    nPoints,
    logLik: fitStats.logLik,
    AIC: fitStats.AIC,
    BIC: fitStats.BIC,
    alpha,
    modelType: 'exponential',
    covarianceSummary
  };
}

/**
 * Bootstrap confidence intervals for dose-response estimates
 *
 * Reference: Carpenter, J. R., & Bithell, J. (2000). Bootstrap confidence
 * intervals: when, which, what? A practical guide for medical statisticians.
 * Statistics in medicine, 19(9), 1141-1164.
 *
 * @param {Array} points - Study data points
 * @param {number} nBootstrap - Number of bootstrap samples (default: 1000)
 * @param {number} ciLevel - Confidence level (default: 0.95)
 * @param {string} modelType - Model to fit ('quadratic', 'spline', 'exponential')
 * @returns {Object} - Bootstrap results with CIs
 */
function bootstrapDoseResponse(points, nBootstrap = 1000, ciLevel = 0.95, modelType = 'quadratic') {
  showProgress(`Running bootstrap (${nBootstrap} iterations)...`);

  // Store original results
  const originalResult = fitModelByType(points, null, modelType);
  const bootstrapEstimates = [];

  // Get unique study IDs
  const studyIds = [...new Set(points.map(p => p.id))];

  // Bootstrap iterations
  for (let iter = 0; iter < nBootstrap; iter++) {
    // Resample studies with replacement
    const bootStudyIds = [];
    for (let i = 0; i < studyIds.length; i++) {
      const randomIdx = Math.floor(Math.random() * studyIds.length);
      bootStudyIds.push(studyIds[randomIdx]);
    }

    // Create bootstrap dataset
    const bootPoints = [];
    for (const studyId of bootStudyIds) {
      const studyPoints = points.filter(p => p.id === studyId);
      bootPoints.push(...studyPoints);
    }

    // Fit model to bootstrap data
    try {
      const bootResult = fitModelByType(bootPoints, null, modelType);
      bootstrapEstimates.push(bootResult.beta);
    } catch (e) {
      // Skip failed iterations
      continue;
    }

    // Update progress periodically
    if ((iter + 1) % 100 === 0) {
      showProgress(`Bootstrap: ${iter + 1}/${nBootstrap}...`);
    }
  }

  // Calculate bootstrap confidence intervals
  if (bootstrapEstimates.length === 0) {
    hideProgress();
    throw new Error('All bootstrap iterations failed. Check data quality and model choice.');
  }

  const alpha = 1 - ciLevel;
  const nParams = originalResult.beta.length;
  const ciLower = [];
  const ciUpper = [];
  const bootstrapSE = [];

  for (let p = 0; p < nParams; p++) {
    const paramEstimates = bootstrapEstimates
      .map(b => b[p])
      .filter(v => Number.isFinite(v))
      .sort((a, b) => a - b);

    if (paramEstimates.length === 0) {
      ciLower.push(NaN);
      ciUpper.push(NaN);
      bootstrapSE.push(NaN);
      continue;
    }

    const lowerIdx = Math.floor(alpha / 2 * paramEstimates.length);
    const upperIdx = Math.ceil((1 - alpha / 2) * paramEstimates.length) - 1;

    ciLower.push(paramEstimates[Math.max(0, lowerIdx)]);
    ciUpper.push(paramEstimates[Math.max(0, upperIdx)]);

    const mean = paramEstimates.reduce((a, b) => a + b, 0) / paramEstimates.length;
    const variance = paramEstimates.reduce((a, b) => a + (b - mean) ** 2, 0) / paramEstimates.length;
    bootstrapSE.push(Math.sqrt(Math.max(variance, 0)));
  }

  hideProgress();

  return {
    originalBeta: originalResult.beta,
    originalSE: originalResult.se,
    bootstrapSE,
    ciLower,
    ciUpper,
    ciLevel,
    nBootstrap: bootstrapEstimates.length,
    modelType
  };
}

function solveGLS(points, tau2Override = null) {
  const prepared = prepareStudyBlocksForModel(points);
  const studyBlocks = prepared.blocks;
  const covarianceSummary = prepared.covarianceSummary;

  if (studyBlocks.length < 1) {
    throw new Error('Need at least 1 study with non-reference dose contrasts');
  }

  const X = [];
  const y = [];
  const V_blocks = [];

  for (const block of studyBlocks) {
    for (const pt of block.modelPoints) {
      X.push([1, pt.dose, pt.dose * pt.dose]);
      y.push(pt.logRR);
    }
    V_blocks.push({ studyId: block.studyId, V: block.V, n: block.modelPoints.length });
  }

  const p = 3;
  const nStudies = V_blocks.length;
  const nPoints = y.length;

  if (nPoints <= p) {
    throw new Error('Insufficient non-reference points for quadratic model estimation');
  }

  let tau2 = tau2Override;
  if (tau2 === null) {
    const tauResult = estimateTau2(X, y, V_blocks);
    tau2 = tauResult.tau2;
  }

  const result = solveGLSWithTau2(X, y, V_blocks, tau2);

  const Q = computeQStat(y, X, result.beta, V_blocks, tau2);
  const df = Math.max((nStudies - 1) * p, 0);
  const I2 = (df > 0 && Q > 0) ? Math.max(0, 100 * (Q - df) / Q) : 0;
  const Qp = df > 0 ? 1 - chiSqCDF(Q, df) : NaN;

  const logLik = -0.5 * (nPoints * Math.log(2 * Math.PI) + result.WRSS + result.detV);
  const AIC = 2 * p - 2 * logLik;
  const BIC = Math.log(Math.max(nPoints, 1)) * p - 2 * logLik;

  return {
    ...result,
    tau2,
    Q,
    df,
    I2,
    Qp,
    nStudies,
    nPoints,
    AIC,
    BIC,
    covarianceSummary
  };
}


/**
 * Solve GLS with tau2 using PROPER block-diagonal matrix inversion
 * This CORRECTLY accounts for within-study correlation
 * Reference: Greenland & Longnecker (1992); Orsini et al. (2006)
 *
 * @param {Array} X - Design matrix
 * @param {Array} y - Outcome vector
 * @param {Array} V_blocks - Array of covariance blocks {V, n, studyId}
 * @param {number} tau2 - Between-study variance
 * @returns {Object} - Results: beta, se, WRSS, detV, var
 */
function solveGLSWithTau2(X, y, V_blocks, tau2) {
  const n = y.length;
  const p = X[0].length;

  // Invert each block of the covariance matrix (WITHIN-STUDY CORRELATION)
  const V_inv_blocks = V_blocks.map(block => {
    const blockSize = block.n;
    const V = new Array(blockSize * blockSize);
    const V_with_tau2 = new Array(blockSize * blockSize);

    // Build covariance matrix with tau2 added to diagonal
    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        V_with_tau2[i * blockSize + j] = block.V[i * blockSize + j] + (i === j ? tau2 : 0);
      }
    }

    // Properly invert the FULL covariance matrix (not just diagonal!)
    const V_inv = invertMatrix(V_with_tau2, blockSize);

    return { V_inv, n: blockSize };
  });

  // Compute X'V^(-1)X and X'V^(-1)y using block-diagonal structure
  const XtVinvX = new Array(p * p).fill(0);
  const XtVinvY = new Array(p).fill(0);

  let row = 0;
  for (const blockInv of V_inv_blocks) {
    const blockSize = blockInv.n;
    const V_inv = blockInv.V_inv;

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        const w_ij = V_inv[i * blockSize + j];  // Full matrix, not just diagonal!

        for (let k = 0; k < p; k++) {
          XtVinvY[k] += w_ij * X[row + i][k] * y[row + j];
          for (let l = 0; l < p; l++) {
            XtVinvX[k * p + l] += w_ij * X[row + i][k] * X[row + j][l];
          }
        }
      }
    }
    row += blockSize;
  }

  // Add small ridge for numerical stability
  for (let i = 0; i < p; i++) {
    XtVinvX[i * p + i] += RIDGE_PENALTY;
  }

  // Solve for beta using generic inversion (supports any p)
  const XtVinvX_inv = invertMatrix(XtVinvX, p);
  const beta = new Array(p).fill(0);
  for (let i = 0; i < p; i++) {
    for (let j = 0; j < p; j++) {
      beta[i] += XtVinvX_inv[i * p + j] * XtVinvY[j];
    }
  }

  // Compute variance-covariance of beta
  const varMatrix = XtVinvX_inv;  // Full variance-covariance matrix
  const se = [];
  for (let i = 0; i < p; i++) {
    se.push(Math.sqrt(Math.max(varMatrix[i * p + i], 0)));
  }

  // Compute weighted residual sum of squares
  let WRSS = 0;
  row = 0;
  for (const blockInv of V_inv_blocks) {
    const blockSize = blockInv.n;
    const V_inv = blockInv.V_inv;

    // Pre-compute all residuals for this block
    const residuals = [];
    for (let i = 0; i < blockSize; i++) {
      let pred_i = 0;
      for (let k = 0; k < p; k++) {
        pred_i += beta[k] * X[row + i][k];
      }
      residuals.push(y[row + i] - pred_i);
    }

    // Compute quadratic form using full V_inv
    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        WRSS += residuals[i] * V_inv[i * blockSize + j] * residuals[j];
      }
    }
    row += blockSize;
  }

  // Compute log determinant of V
  let detV = 0;
  for (const block of V_blocks) {
    const blockSize = block.n;
    for (let i = 0; i < blockSize; i++) {
      detV += Math.log(Math.max(block.V[i * blockSize + i] + tau2, NUMERICAL_TOLERANCE));
    }
  }

  return { beta, se, WRSS, detV, varMatrix };
}

/**
 * Estimate between-study variance (tau²) using TRUE REML
 * Iterative optimization using Fisher scoring algorithm
 *
 * Reference: van Houwelingen, H. C., Arends, L. R., & Stijnen, T. (2002).
 * Advanced methods in meta-analysis: multivariate approach and meta-regression.
 * Statistics in Medicine, 21(4), 589-624.
 *
 * Also: Viechtbauer, W. (2005). Estimating the mean of a normal distribution
 * with known precision. In R Newsletter, 5(1), 11-13.
 *
 * @param {Array} X - Design matrix
 * @param {Array} y - Outcome vector
 * @param {Array} V_blocks - Array of covariance blocks {V, n, studyId}
 * @param {number} maxIter - Maximum iterations (default: 100)
 * @param {number} tol - Convergence tolerance (default: 1e-8)
 * @returns {Object} - {tau2, converged, iterations, logLik}
 */
function estimateTau2REML(X, y, V_blocks, maxIter = 100, tol = 1e-8) {
  const n = y.length;
  const p = X[0].length;
  const K = V_blocks.length;  // Number of studies

  // Start with DL estimator as initial value
  let tau2 = estimateTau2DL(X, y, V_blocks);

  let prevLogLik = -Infinity;
  let converged = false;
  let iter;

  // REML iteration using Fisher scoring
  for (iter = 0; iter < maxIter; iter++) {
    // Get current estimates with this tau2
    const result = solveGLSWithTau2(X, y, V_blocks, tau2);

    // Calculate log-likelihood (REML)
    const logLik = computeREMLLogLik(y, X, result.beta, V_blocks, tau2, p);

    // Check convergence
    if (Math.abs(logLik - prevLogLik) < tol) {
      converged = true;
      break;
    }

    prevLogLik = logLik;

    // Fisher scoring update
    // Calculate score (derivative of log-likelihood w.r.t. tau2)
    const score = computeREMLScore(y, X, result.beta, V_blocks, tau2);

    // Calculate Fisher information
    const fisherInfo = computeREMLFisherInfo(X, V_blocks, tau2);

    // Update tau2
    const step = score / Math.max(fisherInfo, NUMERICAL_TOLERANCE);
    tau2 += step;

    // Ensure tau2 stays non-negative
    if (tau2 < 0) tau2 = 0;
  }

  return {
    tau2: tau2,
    converged: converged,
    iterations: iter + 1,
    logLik: prevLogLik
  };
}

/**
 * Compute REML log-likelihood
 *
 * Reference: Harville, D. A. (1974). Bayesian inference for variance components
 * using only error contrasts. Biometrika, 61(2), 383-385.
 */
function computeREMLLogLik(y, X, beta, V_blocks, tau2, p) {
  const n = y.length;
  let logLik = 0;

  // Log determinant term
  let detV = 0;
  for (const block of V_blocks) {
    const blockSize = block.n;
    for (let i = 0; i < blockSize; i++) {
      detV += Math.log(Math.max(block.V[i * blockSize + i] + tau2, NUMERICAL_TOLERANCE));
    }
  }
  logLik -= 0.5 * detV;

  // Quadratic form term: (y - Xb)' V^(-1) (y - Xb)
  const Q = computeQStat(y, X, beta, V_blocks, tau2);
  logLik -= 0.5 * Q;

  // REML correction term: log |X'V^(-1)X|
  const result = solveGLSWithTau2(X, y, V_blocks, tau2);
  const XtVinvX = computeXtVinvX(X, V_blocks, tau2);
  const logDetXtVinvX = Math.log(Math.max(matrixDeterminant(XtVinvX, p), NUMERICAL_TOLERANCE));
  logLik -= 0.5 * logDetXtVinvX;

  // Constant term
  logLik -= 0.5 * (n - p) * Math.log(2 * Math.PI);

  return logLik;
}

/**
 * Compute REML score (derivative of log-likelihood w.r.t. tau2)
 */
function computeREMLScore(y, X, beta, V_blocks, tau2) {
  const K = V_blocks.length;
  const p = X[0].length;

  // Get residuals and V_inv
  const result = solveGLSWithTau2(X, y, V_blocks, tau2);

  let score = 0;
  let row = 0;

  for (const block of V_blocks) {
    const blockSize = block.n;
    const n = blockSize;

    // Build V and V_inv with current tau2
    const V = new Array(n * n);
    const V_with_tau2 = new Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        V[i * n + j] = block.V[i * n + j];
        V_with_tau2[i * n + j] = block.V[i * n + j] + (i === j ? tau2 : 0);
      }
    }

    const V_inv = invertMatrix(V_with_tau2, n);

    // Pre-compute all residuals for this block
    const residuals = [];
    for (let i = 0; i < n; i++) {
      let pred_i = 0;
      for (let k = 0; k < p; k++) {
        pred_i += beta[k] * X[row + i][k];
      }
      residuals.push(y[row + i] - pred_i);
    }

    // Compute contribution to score
    for (let i = 0; i < n; i++) {
      const resid_i = residuals[i];

      for (let j = 0; j < n; j++) {
        const resid_j = residuals[j];

        // Trace term: tr(V^(-1) * dV/dtau2)
        score += 0.5 * V_inv[i * n + j] * V_inv[i * n + j] * resid_i * resid_j;

        // Derivative of log determinant term
        if (i === j) {
          score -= 0.5 / Math.max(V[i * n + i] + tau2, NUMERICAL_TOLERANCE);
        }
      }
    }

    row += blockSize;
  }

  return score;
}

/**
 * Compute Fisher information for REML
 */
function computeREMLFisherInfo(X, V_blocks, tau2) {
  let fisherInfo = 0;

  for (const block of V_blocks) {
    const blockSize = block.n;
    const n = blockSize;

    // Build V and V_inv with current tau2
    const V_with_tau2 = new Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        V_with_tau2[i * n + j] = block.V[i * n + j] + (i === j ? tau2 : 0);
      }
    }

    const V_inv = invertMatrix(V_with_tau2, n);

    // Fisher information: tr(V^(-1) * dV/dtau2 * V^(-1) * dV/dtau2)
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        fisherInfo += 0.5 * V_inv[i * n + j] * V_inv[i * n + j];
      }
    }
  }

  return fisherInfo;
}

/**
 * Compute X'V^(-1)X matrix
 */
function computeXtVinvX(X, V_blocks, tau2) {
  const p = X[0].length;
  const XtVinvX = new Array(p * p).fill(0);

  let row = 0;
  for (const block of V_blocks) {
    const blockSize = block.n;
    const V_with_tau2 = new Array(blockSize * blockSize);

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        V_with_tau2[i * blockSize + j] = block.V[i * blockSize + j] + (i === j ? tau2 : 0);
      }
    }

    const V_inv = invertMatrix(V_with_tau2, blockSize);

    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        const w_ij = V_inv[i * blockSize + j];
        for (let k = 0; k < p; k++) {
          for (let l = 0; l < p; l++) {
            XtVinvX[k * p + l] += w_ij * X[row + i][k] * X[row + j][l];
          }
        }
      }
    }
    row += blockSize;
  }

  return XtVinvX;
}

/**
 * Compute determinant of 3x3 matrix
 */
function matrixDet3x3(A) {
  return A[0] * (A[4]*A[8] - A[5]*A[7]) -
         A[1] * (A[3]*A[8] - A[5]*A[6]) +
         A[2] * (A[3]*A[7] - A[4]*A[6]);
}

function matrixDet2x2(A) {
  return A[0] * A[3] - A[1] * A[2];
}

function matrixDeterminant(A, n) {
  if (n === 1) return A[0];
  if (n === 2) return matrixDet2x2(A);
  if (n === 3) return matrixDet3x3(A);

  const M = A.slice();
  let det = 1;

  for (let i = 0; i < n; i++) {
    let pivotRow = i;
    let pivotVal = Math.abs(M[i * n + i]);

    for (let r = i + 1; r < n; r++) {
      const candidate = Math.abs(M[r * n + i]);
      if (candidate > pivotVal) {
        pivotVal = candidate;
        pivotRow = r;
      }
    }

    if (pivotVal < NUMERICAL_TOLERANCE) return 0;

    if (pivotRow !== i) {
      for (let c = 0; c < n; c++) {
        const idx1 = i * n + c;
        const idx2 = pivotRow * n + c;
        const tmp = M[idx1];
        M[idx1] = M[idx2];
        M[idx2] = tmp;
      }
      det *= -1;
    }

    const pivot = M[i * n + i];
    det *= pivot;

    for (let r = i + 1; r < n; r++) {
      const factor = M[r * n + i] / pivot;
      if (factor === 0) continue;
      for (let c = i; c < n; c++) {
        M[r * n + c] -= factor * M[i * n + c];
      }
    }
  }

  return det;
}

/**
 * Alias for backward compatibility - now calls TRUE REML
 * Set useREML = true to use iterative REML instead of DL
 */
let useREML = false;  // Use DL for tau2 estimation (simpler, more stable)

/**
 * Estimate between-study variance (tau²) using DerSimonian-Laird method
 * CORRECTED for multivariate meta-analysis with proper degrees of freedom
 *
 * @param {Array} X - Design matrix
 * @param {Array} y - Outcome vector
 * @param {Array} V_blocks - Array of covariance blocks {V, n, studyId}
 * @returns {number} - Estimated tau²
 */
function estimateTau2DL(X, y, V_blocks) {
  const p = X[0].length;
  const K = V_blocks.length;

  if (K <= 1) return 0;

  // Fixed-effect estimate (tau2 = 0)
  const feResult = solveGLSWithTau2(X, y, V_blocks, 0);
  const Q = computeQStat(y, X, feResult.beta, V_blocks, 0);

  // Multivariate meta-analysis degrees of freedom
  const df = Math.max((K - 1) * p, 0);

  // DL denominator based on sum of within-study traces
  let sumTrV = 0;
  for (const block of V_blocks) {
    const n = block.n;
    for (let i = 0; i < n; i++) {
      sumTrV += block.V[i * n + i];
    }
  }

  const denominator = sumTrV - df;
  if (!Number.isFinite(Q) || !Number.isFinite(denominator) || denominator <= 0) {
    return 0;
  }

  return Math.max(0, (Q - df) / denominator);
}

/**
 * Wrapper function that uses the globally selected tau2 estimation method
 * @returns {Object} - tau2 estimate with method info
 */
function estimateTau2(X, y, V_blocks) {
  if (useREML) {
    const remlResult = estimateTau2REML(X, y, V_blocks);
    return {
      tau2: remlResult.tau2,
      method: 'REML',
      converged: remlResult.converged,
      iterations: remlResult.iterations,
      logLik: remlResult.logLik
    };
  } else {
    return {
      tau2: estimateTau2DL(X, y, V_blocks),
      method: 'DL',
      converged: true,
      iterations: 1,
      logLik: null
    };
  }
}

/**
 * Compute Cochran's Q statistic for heterogeneity
 *
 * Reference: Cochran, W. G. (1954). The combination of estimates from
 * different experiments. Biometrics, 10(1), 101-129.
 *
 * @param {Array} y - Outcome vector
 * @param {Array} X - Design matrix
 * @param {Array} beta - Coefficient estimates
 * @param {Array} V_blocks - Array of covariance blocks
 * @param {number} tau2 - Between-study variance
 * @returns {number} - Q statistic
 */
function computeQStat(y, X, beta, V_blocks, tau2) {
  const n = y.length;
  const p = X[0].length;
  let Q = 0;
  let row = 0;

  for (const block of V_blocks) {
    const blockSize = block.n;
    const V = new Array(blockSize * blockSize);

    // Build covariance matrix with tau2
    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        V[i * blockSize + j] = block.V[i * blockSize + j] + (i === j ? tau2 : 0);
      }
    }

    // Invert covariance matrix
    const V_inv = invertMatrix(V, blockSize);

    // Pre-compute all residuals for this block
    const residuals = [];
    for (let i = 0; i < blockSize; i++) {
      let pred_i = 0;
      for (let k = 0; k < p; k++) {
        pred_i += beta[k] * X[row + i][k];
      }
      residuals.push(y[row + i] - pred_i);
    }

    // Compute quadratic form: (y - Xb)' V^(-1) (y - Xb)
    for (let i = 0; i < blockSize; i++) {
      for (let j = 0; j < blockSize; j++) {
        Q += residuals[i] * V_inv[i * blockSize + j] * residuals[j];
      }
    }
    row += blockSize;
  }

  // Q statistic should be non-negative; negative values indicate numerical issues
  return Math.max(0, Q);
}

function solve3x3(A, b) {
  // Solve Ax = b for 3x3 matrix using Cramer's rule
  const det = A[0] * (A[4]*A[8] - A[5]*A[7]) -
              A[1] * (A[3]*A[8] - A[5]*A[6]) +
              A[2] * (A[3]*A[7] - A[4]*A[6]);

  if (Math.abs(det) < DETERMINANT_THRESHOLD) return [0, 0, 0];

  const x0 = (b[0] * (A[4]*A[8] - A[5]*A[7]) -
             b[1] * (A[1]*A[8] - A[2]*A[7]) +
             b[2] * (A[1]*A[5] - A[2]*A[4])) / det;
  const x1 = (A[0] * (b[1]*A[8] - b[2]*A[7]) -
             A[1] * (b[0]*A[8] - b[2]*A[6]) +
             A[2] * (b[0]*A[7] - b[1]*A[6])) / det;
  const x2 = (A[0] * (A[4]*b[2] - A[5]*b[1]) -
             A[1] * (A[3]*b[2] - A[5]*b[0]) +
             A[2] * (A[3]*b[1] - A[4]*b[0])) / det;

  return [x0, x1, x2];
}

function invert3x3(A) {
  const det = A[0] * (A[4]*A[8] - A[5]*A[7]) -
              A[1] * (A[3]*A[8] - A[5]*A[6]) +
              A[2] * (A[3]*A[7] - A[4]*A[6]);

  if (Math.abs(det) < DETERMINANT_THRESHOLD) return A.map(() => 0);

  const inv = new Array(9);
  inv[0] = (A[4]*A[8] - A[5]*A[7]) / det;
  inv[1] = (A[2]*A[7] - A[1]*A[8]) / det;
  inv[2] = (A[1]*A[5] - A[2]*A[4]) / det;
  inv[3] = (A[5]*A[6] - A[3]*A[8]) / det;
  inv[4] = (A[0]*A[8] - A[2]*A[6]) / det;
  inv[5] = (A[2]*A[4] - A[0]*A[5]) / det;
  inv[6] = (A[3]*A[7] - A[4]*A[6]) / det;
  inv[7] = (A[1]*A[6] - A[0]*A[7]) / det;
  inv[8] = (A[0]*A[4] - A[1]*A[3]) / det;

  return inv;
}

function chiSqCDF(x, df) {
  if (x <= 0) return 0;
  if (df === 1) return 2 * (1 - normCDF(Math.sqrt(x)));
  // Wilson-Hilferty approximation
  const z = Math.pow(x / df, 1/3) - (1 - 2 / (9 * df));
  return 1 - normCDF(z * Math.sqrt(df / 2));
}

function normCDF(x) {
  const a1 = .254829592, a2 = -.284496736, a3 = 1.421413741, a4 = -1.453152027, a5 = 1.061405429, p = .3275911;
  const s = x < 0 ? -1 : 1;
  const xAbs = Math.abs(x) / Math.sqrt(2);
  const t = 1 / (1 + p * xAbs);
  const erf = s * (1 - (((((a5 * t + a4) * t + a3) * t + a2) * t + a1) * t) * Math.exp(-xAbs * xAbs));
  return .5 * (1 + erf);
}

function normCDFinv(p) {
  // Beasley-Springer-Moro approximation
  const a = [-3.969683, -1.963510, 0.5742055, 0.6938985, -1.220525, 1.380257, -0.6726213, -0.4548294, 0.2459106];
  const b = [0.0132309, 0.0416717, 0.0097009, -0.0013922, -0.0320185, -0.0098922, 0.0019939, 0.0013942, -0.0015587];

  const q = Math.min(Math.max(p, 1e-10), 1 - 1e-10);
  if (q <= 0.02425) {
    return (a[0] + a[1] / (a[2] + q)) / ((a[3] + q) / (a[4] + q) - (a[5] + q) / (a[6] + q) + (a[7] + q) / (a[8] + q));
  }

  const r = Math.sqrt(-2 * Math.log(1 - q));
  const num = (((((b[5]*r + b[4])*r + b[3])*r + b[2])*r + b[1])*r + b[0]);
  const den = (((((b[8]*r + b[7])*r + b[6])*r + 1))*r + 1);
  return r - num / den;
}

/**
 * Compute t-distribution critical value using approximation
 * Reference: Abramowitz & Stegun, Handbook of Mathematical Functions
 *
 * @param {number} df - Degrees of freedom
 * @param {number} p - Cumulative probability (e.g., 0.975 for 95% two-sided)
 * @returns {number} - t critical value
 */
function tCriticalValue(df, p) {
  if (df <= 0) return normCDFinv(p);
  if (df > 100) return normCDFinv(p);  // Approximate as normal for large df

  // Use approximation based on normal quantile
  const z = normCDFinv(p);

  // Cornish-Fisher expansion for t-distribution
  const g1 = (z * z * z + z) / 4;
  const g2 = (5 * z * z * z * z * z + 16 * z * z * z + 3 * z) / 96;
  const g3 = (3 * z * z * z * z * z * z * z + 19 * z * z * z * z * z + 17 * z * z * z - 15 * z) / 384;

  const t = z + g1 / df + g2 / (df * df) + g3 / (df * df * df);

  return t;
}

// ================================================================
// ANALYSIS & DISPLAY
// ================================================================

function normalizeModelType(modelType, fallback = 'quadratic') {
  const key = String(modelType || '').toLowerCase().trim();
  if (key === 'gls') return 'quadratic';
  if (['linear', 'quadratic', 'spline', 'exponential'].includes(key)) return key;
  return fallback;
}

function normalizeMainModel(modelType) {
  const key = normalizeModelType(modelType);
  if (key === 'linear') return 'linear';
  return 'quadratic';
}

function getMainModelLabel(modelType) {
  const key = normalizeModelType(modelType);
  if (key === 'linear') return 'Linear';
  if (key === 'spline') return 'Restricted Cubic Spline';
  if (key === 'exponential') return 'Exponential';
  return 'Quadratic';
}

function fitModelByType(points, tau2Override = null, modelType = currentMainModel) {
  const resolvedType = normalizeModelType(modelType);
  switch (resolvedType) {
    case 'linear':
      return fitLinearModel(points, tau2Override);
    case 'spline':
      return fitSplineModel(points, tau2Override, 4);
    case 'exponential':
      return fitExponentialModel(points, tau2Override);
    default:
      return solveGLS(points, tau2Override);
  }
}

function setMainModel(modelType) {
  currentMainModel = normalizeMainModel(modelType);
  localStorage.setItem('dose_response_main_model', currentMainModel);
  const select = document.getElementById('mainModelSelect');
  if (select && select.value !== currentMainModel) {
    select.value = currentMainModel;
  }

  if (analysisResults) {
    updateAllDisplays();
  }
  logAudit('main_model_changed', { model: currentMainModel });
}


function setUseREML(v){ useREML = v; }
export { fitLinearModel, fitSplineModel, tCriticalValue, normCDFinv, chiSqCDF, estimateTau2DL, estimateTau2REML, setUseREML };
