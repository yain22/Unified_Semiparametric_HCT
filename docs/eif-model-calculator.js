"use strict";

// Finite iid laws with a binary baseline covariate. At each conditional mean
// the standardized outcome takes -sqrt(3), 0, sqrt(3) with masses 1/6, 2/3, 1/6.
const EifModels = (() => {
  const designs = {
    unrestricted: ["treated", "internal", "external"],
    transportable: ["treated", "internal", "external"],
    rct: ["treated", "internal"],
    ect: ["treated", "external"]
  };
  const allocation = 2 / 3;

  function evaluate(model, settings = {}, record = {}) {
    if (!Object.hasOwn(designs, model)) throw new RangeError("Unknown model.");
    const groups = designs[model];
    const { delta = 0.2, v0 = 1, overlap = "balanced", pattern = "constant" } = settings;
    const { x = 0, group = groups[0], z = 0 } = record;
    if (![0, 1].includes(x) || !groups.includes(group) || ![-1, 0, 1].includes(z)) {
      throw new RangeError("The selected record must belong to the design.");
    }
    if (!Number.isFinite(delta) || !Number.isFinite(v0) || v0 <= 0 ||
        !["balanced", "even", "uneven"].includes(overlap) ||
        !["constant", "concentrated"].includes(pattern)) {
      throw new RangeError("Invalid parameters for the finite law.");
    }

    const isEct = model === "ect";
    const isRct = model === "rct";
    const isUnrestricted = model === "unrestricted";
    const q = isEct ? 0.4 : isRct ? 1 : 0.5;
    const e = isEct ? 1 : allocation;
    // Under the reference HCT, X is even in the mixture and q_H = 1/2.
    // Retain F(X|trial) and F(X|external) across designs, then use each
    // design's own q to obtain its mixture and source propensity by Bayes.
    const strata = [0, 1].map(xValue => {
      const trialMass = overlap === "uneven" ? (1 + xValue) / 3 : 0.5;
      const externalMass = 1 - trialMass;
      const fx = q * trialMass + (1 - q) * externalMass;
      const p = q * trialMass / fx;
      const baseline = 0.4 * xValue;
      const drift = model === "transportable" ? 0 :
        pattern === "concentrated" ? 2 * delta * xValue : delta;
      const means = { treated: 1 + baseline };
      if (!isEct) means.internal = baseline;
      if (!isRct) means.external = baseline - drift;
      let lambda;
      if (isEct) lambda = 0;
      else if (isRct) lambda = 1;
      else {
        const internalPrecision = p * (1 - e);
        const externalPrecision = (1 - p) / v0;
        lambda = internalPrecision / (internalPrecision + externalPrecision);
      }
      let control = 0;
      if (!isEct) control += lambda * means.internal;
      if (!isRct) control += (1 - lambda) * means.external;
      return { trialMass, fx, p, lambda, means, t: means.treated - control };
    });
    const target = strata.reduce((sum, s) => sum + s.trialMass * s.t, 0);
    const trial = isEct ? null : strata.reduce((sum, s) =>
      sum + s.trialMass * (s.means.treated - s.means.internal), 0);
    const gap = isEct ? null : isUnrestricted ? target - trial : 0;

    function contribution(stratum, arm, atom) {
      const { p, lambda, means, t } = stratum;
      const S = Number(arm !== "external");
      const variance = arm === "external" ? v0 : 1;
      const residual = Math.sqrt(3 * variance) * atom;
      const y = means[arm] + residual;
      const parts = {
        standardization: S * (t - target) / q,
        treated: arm === "treated" ? residual / (q * e) : 0,
        internal: arm === "internal" ? -lambda * residual / (q * (1 - e)) : 0,
        external: arm === "external" ? -(1 - lambda) * p * residual / (q * (1 - p)) : 0,
        source: 0, internalVariance: 0, externalVariance: 0
      };
      if (isUnrestricted) {
        const coefficient = (means.internal - means.external) * lambda * (1 - lambda) / q;
        parts.source = -coefficient * (S - p) / (1 - p);
        const squared = residual * residual / variance - 1;
        if (arm === "internal") parts.internalVariance = coefficient * squared / (1 - e);
        if (arm === "external") parts.externalVariance = -coefficient * p * squared / (1 - p);
      }
      const phi = Object.values(parts).reduce((a, b) => a + b, 0);
      let gapPhi = isEct ? null : 0;
      if (isUnrestricted) {
        let phiTrial = S * (means.treated - means.internal - trial) / q;
        if (arm === "treated") phiTrial += residual / (q * e);
        if (arm === "internal") phiTrial -= residual / (q * (1 - e));
        gapPhi = phi - phiTrial;
      }
      const sourceMass = arm === "external" ? 1 - p : arm === "internal" ? p * (1 - e) : p * e;
      const probability = stratum.fx * sourceMass * (atom === 0 ? 2 / 3 : 1 / 6);
      return { y, residual, probability, phi, gapPhi, parts };
    }

    let meanPhi = 0, bound = 0, gapBound = isEct ? null : 0;
    for (const s of strata) {
      for (const arm of groups) {
        for (const atom of [-1, 0, 1]) {
          const a = contribution(s, arm, atom);
          meanPhi += a.probability * a.phi;
          bound += a.probability * a.phi ** 2;
          if (!isEct) gapBound += a.probability * a.gapPhi ** 2;
        }
      }
    }
    const s = strata[x];
    const a = contribution(s, group, z);
    const terms = Object.entries(a.parts).map(([key, value]) => ({ key, value }));
    return { q, p: s.p, e, lambda: s.lambda,
      delta: isEct || isRct ? null : s.means.internal - s.means.external,
      target, trial, gap, y: a.y, residual: a.residual, probability: a.probability,
      phi: a.phi, gapPhi: a.gapPhi, bound, gapBound, meanPhi, terms, groups: groups.slice() };
  }

  return { evaluate };
})();
if (typeof module !== "undefined" && module.exports) module.exports = EifModels;
if (typeof window !== "undefined") window.EifModels = EifModels;
