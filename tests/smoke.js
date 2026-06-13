/*
 * Minimal smoke test for Dose Response Pro (single-file HTML app).
 *
 * Loads index.html, extracts the main inline <script>, evaluates it in a
 * sandboxed vm with stub DOM globals, and exercises the pure numerical
 * helpers of the stat engine (invertMatrix, normCDF, simpleOLS) against
 * known-correct values.
 *
 * Run: node tests/smoke.js
 * Exit 0 = all assertions passed; exit 1 = failure.
 */

const fs = require('fs');
const path = require('path');
const vm = require('vm');
const assert = require('assert');

const htmlPath = path.join(__dirname, '..', 'index.html');
const html = fs.readFileSync(htmlPath, 'utf8').replace(/\r\n/g, '\n');

// --- Structural checks on the shipped file --------------------------------
const openDivs = (html.match(/<div[\s>]/g) || []).length;
const closeDivs = (html.match(/<\/div>/g) || []).length;
assert.strictEqual(openDivs, closeDivs, `div balance mismatch: ${openDivs} open vs ${closeDivs} close`);
assert.ok(!/﻿/.test(html), 'index.html must not contain a BOM');

// Extract the main inline <script> ... </script> (the engine, before report template).
// The engine block starts at the first top-level <script> after the body markup.
const scriptOpenIdx = html.indexOf('\n<script>\n');
assert.ok(scriptOpenIdx !== -1, 'could not locate main <script> block');
const bodyStart = scriptOpenIdx + '\n<script>\n'.length;
// The engine block is closed by the next literal </script> in the file.
const scriptCloseIdx = html.indexOf('</script>', bodyStart);
assert.ok(scriptCloseIdx !== -1, 'could not locate closing </script>');
const engineSource = html.slice(bodyStart, scriptCloseIdx);

assert.ok(engineSource.includes('function invertMatrix'), 'engine must define invertMatrix');
assert.ok(engineSource.includes('function solveGLSLinear'), 'engine must define solveGLSLinear');
assert.ok(engineSource.includes('function estimateTau2Linear'), 'engine must define estimateTau2Linear');

// --- Evaluate pure helpers in a sandbox -----------------------------------
// Provide enough stub globals that top-level declarations evaluate without a DOM.
const noop = () => {};
// A deeply chainable stub: any property access returns the same callable proxy,
// so arbitrary top-level DOM init code (document.x.y(...).z = ...) is a no-op.
function makeChainable() {
  const fn = function () { return proxy; };
  const proxy = new Proxy(fn, {
    get(_t, prop) {
      if (prop === Symbol.toPrimitive || prop === 'toString' || prop === 'valueOf') return () => '';
      if (prop === 'length') return 0;
      if (prop === 'style' || prop === 'classList' || prop === 'dataset') return proxy;
      return proxy;
    },
    set: () => true,
    apply: () => proxy,
  });
  return proxy;
}
const chainable = makeChainable();
const sandbox = {
  console: { log: noop, warn: noop, error: noop },
  Math, Number, Array, Object, JSON, isFinite, isNaN, parseFloat, parseInt,
  Date, Proxy, Map, Set, String, Boolean, RegExp, Symbol,
  document: chainable,
  navigator: { userAgent: 'node-smoke' },
  location: { href: '', search: '' },
  setTimeout: noop, clearTimeout: noop, requestAnimationFrame: noop,
  addEventListener: noop, removeEventListener: noop,
  localStorage: { getItem: () => null, setItem: noop, removeItem: noop },
};
sandbox.window = sandbox;
sandbox.globalThis = sandbox;

const context = vm.createContext(sandbox);
// Function declarations are hoisted, so the export line at the top of the IIFE
// captures them even if later top-level DOM-init statements throw. We swallow
// init-time errors (no DOM here) but still surface them if exports are missing.
const wrapped =
  '(function () {\n' +
  '  try {\n' +
  '    globalThis.__exports = function () {\n' +
  '      return { invertMatrix, normCDF, simpleOLS };\n' +
  '    };\n' +
  engineSource +
  '\n  } catch (e) { globalThis.__initError = e; }\n' +
  '})();\n';
vm.runInContext(wrapped, context, { timeout: 10000 });
sandbox.__exports = sandbox.__exports ? sandbox.__exports() : {};

const { invertMatrix, normCDF, simpleOLS } = sandbox.__exports;
assert.strictEqual(typeof invertMatrix, 'function', 'invertMatrix not exported');
assert.strictEqual(typeof normCDF, 'function', 'normCDF not exported');
assert.strictEqual(typeof simpleOLS, 'function', 'simpleOLS not exported');

// invertMatrix: inverse of [[4,3],[6,3]] is [[-0.5,0.5],[1,-0.6667]].
const inv = invertMatrix([4, 3, 6, 3], 2);
assert.ok(Math.abs(inv[0] - (-0.5)) < 1e-9, `inv[0] wrong: ${inv[0]}`);
assert.ok(Math.abs(inv[1] - 0.5) < 1e-9, `inv[1] wrong: ${inv[1]}`);
assert.ok(Math.abs(inv[2] - 1.0) < 1e-9, `inv[2] wrong: ${inv[2]}`);
assert.ok(Math.abs(inv[3] - (-2 / 3)) < 1e-9, `inv[3] wrong: ${inv[3]}`);

// normCDF: known values.
assert.ok(Math.abs(normCDF(0) - 0.5) < 1e-9, `normCDF(0) wrong: ${normCDF(0)}`);
assert.ok(Math.abs(normCDF(1.959964) - 0.975) < 1e-4, `normCDF(1.96) wrong: ${normCDF(1.959964)}`);
assert.ok(normCDF(-3) < 0.01 && normCDF(3) > 0.99, 'normCDF tails wrong');

// simpleOLS: perfect line y = 2x + 1 -> slope 2, intercept 1.
const ols = simpleOLS([1, 2, 3, 4], [3, 5, 7, 9]);
assert.ok(ols && Math.abs(ols.slope - 2) < 1e-9, `OLS slope wrong: ${ols && ols.slope}`);
assert.ok(Math.abs(ols.intercept - 1) < 1e-9, `OLS intercept wrong: ${ols && ols.intercept}`);
// k<3 guard returns null.
assert.strictEqual(simpleOLS([1, 2], [1, 2]), null, 'simpleOLS must reject n<3');

console.log('smoke: 11 assertions passed');
