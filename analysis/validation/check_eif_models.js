"use strict";

const assert = require("node:assert/strict");
const { evaluate } = require("../../docs/eif-model-calculator.js");
const { inspect } = require("../../docs/eif-lab.js");

function close(actual, expected, tolerance = 1e-10) {
  assert.ok(Number.isFinite(actual) && Math.abs(actual - expected) < tolerance,
    `Expected ${expected}; found ${actual}`);
}

function atoms(model, settings) {
  const result = [];
  const { groups, q } = evaluate(model, settings);
  for (const x of [0, 1]) {
    const trialMass = settings.overlap === "uneven" ? (1 + x) / 3 : 0.5;
    const fx = q * trialMass + (1 - q) * (1 - trialMass);
    for (const group of groups) {
      for (const z of [-1, 0, 1]) {
        const a = evaluate(model, settings, { x, group, z });
        result.push({ ...a, x, fx, group, z, S: Number(group !== "external") });
      }
    }
  }
  return result;
}

// Construct smooth submodels on the fixed outcome atoms. The linear tilts
// preserve each conditional total exactly. In the common-mean model the two
// control means move by exactly the same amount, even at unequal variances.
function scoreFor(a, settings, path) {
  if (path === "case mix") return a.x - (a.x ? a.fx : 1 - a.fx);
  if (path === "membership") return (1 + 0.2 * a.x) * (a.S - a.p);
  const variance = a.group === "external" ? settings.v0 : 1;
  if (path === "control mean" && a.group !== "treated") return a.residual / variance;
  if (path === a.group + " mean") return a.residual / variance;
  if (path === a.group + " variance") return a.residual ** 2 / variance - 1;
  return 0;
}

// Recover the target only from probabilities and first two outcome moments.
function functional(atoms, model, epsilon, settings, path) {
  let numerator = 0, q = 0;
  for (const x of [0, 1]) {
    const cells = {};
    for (const a of atoms.filter(a => a.x === x)) {
      const weight = a.probability * (1 + epsilon * scoreFor(a, settings, path));
      const cell = cells[a.group] ||= { mass: 0, mean: 0, moment2: 0 };
      cell.mass += weight;
      cell.mean += weight * a.y;
      cell.moment2 += weight * a.y ** 2;
    }
    for (const c of Object.values(cells)) {
      c.mean /= c.mass;
      c.variance = c.moment2 / c.mass - c.mean ** 2;
    }
    const { treated, internal, external } = cells;
    const trialMass = treated.mass + (internal?.mass ?? 0);
    let comparator;
    if (model === "ect") comparator = external.mean;
    else if (model === "rct" || model === "transportable") comparator = internal.mean;
    else {
      const w = internal.mass / internal.variance;
      const v = external.mass / external.variance;
      comparator = (w * internal.mean + v * external.mean) / (w + v);
    }
    numerator += trialMass * (treated.mean - comparator);
    q += trialMass;
  }
  return numerator / q;
}

let paths = 0;
for (const overlap of ["balanced", "uneven"]) {
  for (const v0 of [0.25, 1, 4]) {
    for (const pattern of ["constant", "concentrated"]) {
      const settings = { overlap, v0, pattern, delta: -0.3 };
      for (const model of ["unrestricted", "transportable", "rct", "ect"]) {
        const all = atoms(model, settings);
        const r = all[0];
        close(all.reduce((s, a) => s + a.probability, 0), 1);
        close(all.reduce((s, a) => s + a.probability * a.phi, 0), 0);
        close(r.target, functional(all, model, 0, settings, ""));
        for (const a of all) {
          close(a.terms.reduce((s, t) => s + t.value, 0), a.phi);
          if (model === "ect") {
            assert.equal(a.trial, null);
            assert.equal(a.gapPhi, null);
            assert.equal(a.gap, null);
            assert.equal(a.gapBound, null);
            const m0 = 0.4 * a.x - (pattern === "concentrated" ? 2 * settings.delta * a.x : settings.delta);
            close(a.phi, (a.S * (a.y - m0 - a.target) -
              (1 - a.S) * a.p / (1 - a.p) * (a.y - m0)) / a.q);
          }
        }
        const directions = ["case mix", "membership", "control mean", "treated mean",
          "internal variance", "external variance"];
        // Changing just one control mean is an admissible path in the
        // unrestricted HCT, but would leave the common-mean model.
        if (model === "unrestricted") directions.push("internal mean", "external mean");
        for (const path of directions) {
          const h = 1e-5;
          const fd = (functional(all, model, h, settings, path) -
            functional(all, model, -h, settings, path)) / (2 * h);
          const product = all.reduce((s, a) => s + a.probability * a.phi * scoreFor(a, settings, path), 0);
          close(product, fd, 2e-8);
          if (path === "internal mean" || path === "external mean") {
            const share = path === "internal mean" ? a => a.lambda : a => 1 - a.lambda;
            // Both outcome tilts preserve variances to first order by symmetry.
            const meanShare = all.reduce((s, a) => s + a.probability * a.S * share(a) / a.q, 0);
            close(product, -meanShare);
          }
          paths++;
        }
        if (model === "unrestricted") {
          const previous = inspect({ ...settings, target: "tce", c: 0, path: "source" });
          close(r.bound, previous.totals.tceBound);
          close(r.gapBound, previous.totals.gapBound);
          close(r.gap, previous.population.gap);
        } else if (model === "transportable" || model === "rct") {
          close(r.gapBound, 0);
          close(r.gap, 0);
        }
      }
    }
  }
}
const hct = evaluate("unrestricted");
close(hct.lambda, 0.25);
close(hct.target, 1.15);
close(hct.bound, 4.528125);
close(evaluate("rct").bound, 4.5);
close(evaluate("ect").bound, 1 / 0.4 + 1 / 0.6);
close(evaluate("transportable").bound, 4.5);
const unrestrictedAtEquality = evaluate("unrestricted", { delta: 0 });
close(unrestrictedAtEquality.gap, 0);
close(unrestrictedAtEquality.gapBound, 4.5);
assert.throws(() => evaluate("ect", {}, { group: "internal" }), RangeError);
assert.throws(() => evaluate("rct", {}, { group: "external" }), RangeError);
assert.equal(paths, 312);
console.log(`Verified four EIFs and ${paths} pathwise derivatives.`);
