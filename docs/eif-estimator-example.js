/* Evaluate the paper's estimators on the six illustrative records in Figure S1. */
(() => {
  "use strict";
  const trialRecords = typeof module !== "undefined" && module.exports
    ? require("./eif-record-table.js") : window.EifRecordTable;
  const sum = values => values.reduce((total, value) => total + value, 0);
  const mean = values => sum(values) / values.length;
  const e = 0.5;
  const pH = x => 0.55 + 0.15 * x;
  const treatedMean = x => 8 + 2 * x;
  const internalMean = x => 7 + 1.5 * x;
  const externalMean = x => 6.5 + 1.5 * x;
  const vInternal = 1;
  const vExternal = 1.5;
  const rho = vInternal / vExternal;

  // Simple OLS with intercept. Used both for the pooled control regression
  // and for residualizing treatment on X in the ANCOVA coefficient.
  function fitLine(records, field) {
    const xbar = mean(records.map(row => row.x));
    const ybar = mean(records.map(row => row[field]));
    const slope = sum(records.map(row => (row.x - xbar) * (row[field] - ybar))) /
      sum(records.map(row => (row.x - xbar) ** 2));
    return x => ybar + slope * (x - xbar);
  }

  function unrestricted(row) {
    const {s, a, x, y} = row;
    const p = pH(x);
    const odds = p / (1 - p);
    const lambda = p * (1 - e) / (p * (1 - e) + (1 - p) * rho);
    const m1 = treatedMean(x), m10 = internalMean(x), m0 = externalMean(x);
    const c = (m10 - m0) * lambda * (1 - lambda);
    return [
      s * (m1 - (lambda * m10 + (1 - lambda) * m0)),
      a * (y - m1) / e,
      -s * (1 - a) * lambda * (y - m10) / (1 - e),
      -(1 - s) * odds * (1 - lambda) * (y - m0),
      -c * (s - p) / (1 - p),
      s * (1 - a) * c * ((y - m10) ** 2 / vInternal - 1) / (1 - e),
      -(1 - s) * c * odds * ((y - m0) ** 2 / vExternal - 1)
    ];
  }

  function restricted(row, pooledMean) {
    const {s, a, x, y} = row;
    const p = pH(x);
    const d = p * (1 - e) + (1 - p) * rho;
    const mp = pooledMean(x);
    return [
      s * (treatedMean(x) - mp),
      a * (y - treatedMean(x)) / e,
      -p * (s * (1 - a) + (1 - s) * rho) * (y - mp) / d
    ];
  }

  function ect(row) {
    const {s, x, y} = row;
    // Keeping trial-treated and external records changes the source odds.
    const p = pH(x);
    const oddsE = e * p / (1 - p);
    const residual = y - externalMean(x);
    return [s * residual, -(1 - s) * oddsE * residual];
  }

  function trialWeights(rows, adjustment) {
    const abar = mean(rows.map(row => row.a));
    const fittedA = adjustment === "ancova" ? fitLine(rows, "a") : () => abar;
    const residual = rows.map(row => row.a - fittedA(row.x));
    const ss = sum(residual.map(value => value ** 2));
    return residual.map(value => value / ss);
  }

  const notes = {
    unrestricted: "Evaluate all seven terms at X = X₁. For this illustration, hold the nuisance predictions fixed as you vary the two external outcomes: pH(x) = 0.55 + 0.15x, e = 1/2, m₁(x) = 8 + 2x, m₁,₀(x) = 7 + 1.5x, m₀(x) = 6.5 + 1.5x, v₁,₀ = 1, and v₀ = 1.5. The two variance terms use the observed squared control residuals.",
    transportable: "Fit one pooled control regression of Y on an intercept and X₁ using the four control records. This regression is refitted when the external outcomes change. Hold the supplied predictions pH(x) = 0.55 + 0.15x and m₁(x) = 8 + 2x fixed, with e = 1/2 and working variance ratio ρ = 1/1.5.",
    ect: "Keep the two treated and two external records. Holding the supplied external mean m₀,E(x) = 6.5 + 1.5x fixed, evaluate their residuals at X = X₁. For this record restriction the ECT source odds are e pH(x)/{1 − pH(x)}, with e = 1/2 and pH(x) = 0.55 + 0.15x.",
    anova: "There are two participants in each randomized arm. The treated outcomes receive weight +½ each and the control outcomes receive weight −½ each.",
    ancova: "Here the adjustment is for X₁. Fit a least-squares regression of treatment on an intercept and X₁ across the four trial records. The Frisch–Waugh–Lovell formula uses its treatment residuals to give the weights below."
  };
  const labels = {
    unrestricted: ["Standardization", "Treated residual", "Internal residual", "External residual", "Source probability", "Internal variance", "External variance"],
    transportable: ["Standardization", "Treated residual", "Pooled control residual"],
    ect: ["Trial residual", "External odds correction"],
    rct: ["OLS weight × outcome"]
  };

  function evaluate(model, {externalShift = 0, adjustment = "ancova"} = {}) {
    if (!Object.hasOwn(labels, model) || !Number.isFinite(externalShift) ||
        !["anova", "ancova"].includes(adjustment)) throw new RangeError("Invalid estimator or settings.");
    const rows = trialRecords.getRecords().filter(row =>
      (model !== "rct" || row.s === 1) && (model !== "ect" || row.a === 1 || row.s === 0));
    for (const row of rows) {
      if (!row.s) row.y += externalShift;
      row.group = row.a ? "Trial treated" : row.s ? "Internal control" : "External control";
    }
    const pooledMean = model === "transportable" ? fitLine(rows.filter(row => !row.a), "y") : null;
    const weights = model === "rct" ? trialWeights(rows, adjustment) : null;
    for (const [i, row] of rows.entries()) {
      row.weight = weights ? weights[i] : null;
      if (model === "unrestricted") row.parts = unrestricted(row);
      if (model === "transportable") row.parts = restricted(row, pooledMean);
      if (model === "ect") row.parts = ect(row);
      if (model === "rct") row.parts = [row.weight * row.y];
      row.value = sum(row.parts);
    }
    const numerator = sum(rows.map(row => row.value));
    const denominator = weights ? 1 : sum(rows.map(row => row.s));
    const denominatorLabel = weights ? "Already scaled OLS weights" :
      model === "ect" ? "Treated trial participants" : "Trial participants";
    return {estimate: numerator / denominator, numerator, denominator, denominatorLabel, rows,
      termLabels: labels[model].slice(), note: notes[model === "rct" ? adjustment : model]};
  }

  if (typeof window !== "undefined") window.EifEstimatorExample = {evaluate};
  if (typeof module !== "undefined" && module.exports) module.exports = {evaluate};
})();
