"use strict";

// Enumerate an iid hybrid law. Each outcome distribution uses three fixed
// support points along a path. Exponential tilting changes their probabilities.
const stratumValues = [0, 1];
const groups = ["treated", "internal", "external"];
const zValues = [-Math.sqrt(3), 0, Math.sqrt(3)];
const zMass = [1 / 6, 2 / 3, 1 / 6];
const allocation = 2 / 3;

function shiftedOdds(p, offset) {
  const odds = p / (1 - p) * Math.exp(offset);
  return odds / (1 + odds);
}

function direction(path, x, group, z, settings) {
  if (path === "external-variance" && group === "external") return (z * z - 1) / 2;
  if (path === "external-mean" && group === "external") return z / Math.sqrt(settings.v0);
  if (path === "internal-variance" && group === "internal") return (z * z - 1) / 2;
  return 0;
}

function model(settings, path = "", epsilon = 0) {
  const px1 = shiftedOdds(0.5, path === "case-mix" ? epsilon : 0);
  return stratumValues.map(x => {
    const pBase = settings.overlap === "uneven" ? (1 + x) / 3 : 0.5;
    const p = shiftedOdds(pBase, path === "source" ? epsilon : 0);
    const drift = settings.pattern === "concentrated" ? 2 * settings.delta * x : settings.delta;
    const trialBaseline = 0.4 * x;
    const laws = {};
    for (const group of groups) {
      const center = group === "treated" ? trialBaseline + 1 :
        group === "external" ? trialBaseline - drift : trialBaseline;
      const spread = group === "external" ? Math.sqrt(settings.v0) : 1;
      const atoms = zValues.map((z, i) => ({
        y: center + spread * z,
        score: direction(path, x, group, z, settings),
        mass: zMass[i] * Math.exp(epsilon * direction(path, x, group, z, settings))
      }));
      const sum = atoms.reduce((acc, a) => acc + a.mass, 0);
      atoms.forEach(a => { a.mass /= sum; });
      const mean = atoms.reduce((acc, a) => acc + a.mass * a.y, 0);
      const variance = atoms.reduce((acc, a) => acc + a.mass * (a.y - mean) ** 2, 0);
      laws[group] = { atoms, mean, variance };
    }
    const internalPrecision = p * (1 - allocation) / laws.internal.variance;
    const externalPrecision = (1 - p) / laws.external.variance;
    const lambda = internalPrecision / (internalPrecision + externalPrecision);
    const tTrial = laws.treated.mean - laws.internal.mean;
    const delta = laws.internal.mean - laws.external.mean;
    return { x, fx: x === 0 ? 1 - px1 : px1, p, laws, lambda, delta, tTrial,
      tBorrow: tTrial + (1 - lambda) * delta };
  });
}

function targets(strata) {
  const q = strata.reduce((a, s) => a + s.fx * s.p, 0);
  const trial = strata.reduce((a, s) => a + s.fx * s.p * s.tTrial / q, 0);
  const tce = strata.reduce((a, s) => a + s.fx * s.p * s.tBorrow / q, 0);
  return { q, trial, tce, gap: tce - trial };
}

function inspect(settings) {
  const { path, target, c } = settings;
  const strata = model(settings, path);
  const population = targets(strata);
  const { q, tce, trial } = population;
  const totals = { fixed: 0, source: 0, internal: 0, external: 0, derivative: 0,
    bound: 0, augmented: 0, centered: 0, orthogonal: 0, trialBound: 0, gapBound: 0, tceBound: 0 };
  for (const { x, fx, p, laws, lambda, delta, tTrial, tBorrow } of strata) {
    for (const group of groups) {
      const S = Number(group !== "external"), A = Number(group === "treated");
      const IC = Number(group === "internal"), EC = 1 - S;
      const groupMass = A * p * allocation + IC * p * (1 - allocation) + EC * (1 - p);
      const law = laws[group];
      for (const atom of law.atoms) {
        const w = fx * groupMass * atom.mass;
        const resid = atom.y - law.mean;
        const score = atom.score + (path === "source" ? S - p : 0) +
          (path === "case-mix" ? x - 0.5 : 0);
        const trialIF = (S * (tTrial - trial) + A * resid / allocation -
          IC * resid / (1 - allocation)) / q;
        const fixedIF = (S * (tBorrow - tce) + A * resid / allocation -
          IC * lambda * resid / (1 - allocation) -
          EC * p / (1 - p) * (1 - lambda) * resid) / q;
        const b = delta * lambda * (1 - lambda) / q;
        const sourceIF = -b * (S - p) / (1 - p);
        const squared = resid * resid / law.variance - 1;
        const internalIF = b * IC * squared / (1 - allocation);
        const externalIF = -b * EC * p / (1 - p) * squared;
        const tceIF = fixedIF + sourceIF + internalIF + externalIF;
        const gapIF = tceIF - trialIF;
        const phi = target === "trial" ? trialIF : target === "gap" ? gapIF : tceIF;
        const assignment = S * (A - allocation);
        const fixedTarget = target === "trial" ? trialIF : target === "gap" ? fixedIF - trialIF : fixedIF;
        totals.fixed += w * fixedTarget * score;
        if (target !== "trial") {
          totals.source += w * sourceIF * score;
          totals.internal += w * internalIF * score;
          totals.external += w * externalIF * score;
        }
        totals.derivative += w * phi * score;
        totals.bound += w * phi * phi;
        totals.augmented += w * (phi + c * assignment) ** 2;
        totals.centered += w * phi;
        totals.orthogonal += w * assignment * score;
        totals.trialBound += w * trialIF ** 2;
        totals.gapBound += w * gapIF ** 2;
        totals.tceBound += w * tceIF ** 2;
      }
    }
  }
  const h = 0.0001;
  const forward = targets(model(settings, path, h))[target];
  const backward = targets(model(settings, path, -h))[target];
  return { strata, population, totals, finiteDifference: (forward - backward) / (2 * h) };
}

const pathNotes = {
  "external-variance": "Tilt external outcomes toward both tails. The external mean stays fixed. Its log variance grows at unit speed at ε = 0, so the precision rule gives more weight to internal controls.",
  "external-mean": "Tilt external outcomes toward larger values. At ε = 0 the external mean increases by one outcome unit per unit ε. Symmetry makes its first-order variance change zero.",
  "internal-variance": "Tilt internal control outcomes toward both tails. Their mean stays fixed and their log variance grows at unit speed at ε = 0. The precision rule shifts weight toward external controls.",
  source: "Increase the log odds of trial membership by ε in both strata. This can change both the borrowing share and the distribution of X among trial participants.",
  "case-mix": "Increase the population log odds of X = 1 by ε. The conditional outcome laws and borrowing shares stay fixed. Choose discrepancy concentrated at X = 1 to see the effect of changing the trial mix."
};
const directions = {
  "external-variance": "I₀ (Z² − 1) / 2",
  "external-mean": "I₀ Z / √v₀",
  "internal-variance": "I₁,₀ (Z² − 1) / 2",
  source: "S − p(X)",
  "case-mix": "X − ½"
};
const targetNames = { tce: "Borrowing target", trial: "Trial effect", gap: "Causal-gap" };

function showNumber(value, digits = 4) {
  if (Math.abs(value) < 0.5 * 10 ** (-digits)) value = 0;
  return value.toFixed(digits).replace("-", "−");
}

function drawPath(settings, baseline, gradient) {
  const container = document.getElementById("chart-drawing");
  container.replaceChildren();
  function add(name, attributes, label) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", name);
    for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
    if (label !== undefined) node.textContent = label;
    container.append(node);
  }
  const half = 0.5;
  const points = Array.from({ length: 49 }, (_, i) => {
    const eps = -half + i / 48;
    return { eps, change: targets(model(settings, settings.path, eps))[settings.target] - baseline };
  });
  const max = Math.max(0.01, Math.abs(half * gradient), ...points.map(p => Math.abs(p.change))) * 1.08;
  const left = 100, right = 578, top = 35, bottom = 162;
  const x = value => left + (value + half) * (right - left);
  const y = value => top + (1 - value / max) * (bottom - top) / 2;
  add("text", { x: left, y: 19 }, "Change in target (outcome units)");
  for (const v of [-max, 0, max]) {
    add("line", {x1: left, x2: right, y1: y(v), y2: y(v), stroke: "#ccd8db"});
    add("text", {x: left - 8, y: y(v) + 4, "text-anchor": "end"}, showNumber(v, 3));
  }
  for (const v of [-half, 0, half]) {
    add("line", {x1: x(v), x2: x(v), y1: top, y2: bottom, stroke: "#e4ebec"});
    add("text", {x: x(v), y: bottom + 18, "text-anchor": "middle"}, showNumber(v, 1));
  }
  add("text", {x: (left + right) / 2, y: 202, "text-anchor": "middle"}, "Perturbation ε");
  add("polyline", {points: points.map(p => x(p.eps) + "," + y(p.change)).join(" "),
    stroke: "#115b70", "stroke-width": 3.5, fill: "none"});
  add("line", {x1: left, x2: right, y1: y(-half * gradient), y2: y(half * gradient),
    stroke: "#88531c", "stroke-width": 2, "stroke-dasharray": "7 5"});
  add("circle", {cx: x(0), cy: y(0), r: 3.3, fill: "#263b47"});
}

function startLab() {
  const get = id => document.getElementById(id);
  const write = (id, value) => { get(id).textContent = value; };
  function render() {
    const settings = {
      delta: Number(get("delta").value), v0: Number(get("v0").value),
      overlap: get("overlap").value, pattern: get("pattern").value,
      path: get("path").value, target: get("target").value,
      c: Number(get("penalty").value), n: Number(get("size").value)
    };
    const { strata, population, totals, finiteDifference } = inspect(settings);
    write("delta-value", showNumber(settings.delta, 2));
    write("v0-value", showNumber(settings.v0, 2));
    write("c-value", showNumber(settings.c, 2));
    const cells = get("strata-rows");
    cells.replaceChildren();
    for (const s of strata) {
      const row = cells.insertRow();
      const th = document.createElement("th");
      th.scope = "row";
      th.textContent = s.x;
      row.append(th);
      for (const value of [s.p, s.fx * s.p / population.q, s.delta, s.lambda, (1 - s.lambda) * s.delta]) {
        row.insertCell().textContent = showNumber(value, 3);
      }
    }
    for (const key of ["trial", "tce", "gap"]) write(key + "-value", showNumber(population[key], 3));
    drawPath(settings, population[settings.target], totals.derivative);
    write("path-note", pathNotes[settings.path]);
    write("path-score", directions[settings.path]);
    write("selected-target", targetNames[settings.target]);
    write("gradient", showNumber(totals.derivative, 5));
    write("difference", showNumber(finiteDifference, 5));
    for (const key of ["fixed", "source", "internal", "external"]) write(key + "-part", showNumber(totals[key], 5));
    write("orthogonal-pair", showNumber(totals.orthogonal, 5));
    write("efficient-bound", showNumber(totals.bound));
    write("alternative-var", showNumber(totals.augmented));
    write("efficient-se", showNumber(Math.sqrt(totals.bound / settings.n), 4));
    write("alternative-se", showNumber(Math.sqrt(totals.augmented / settings.n), 4));
    write("extra-var", showNumber(totals.augmented - totals.bound));
    get("efficient-bar").style.width = (100 * totals.bound / totals.augmented) + "%";
    get("restriction").textContent = settings.delta === 0
      ? "At this law the two control means agree in both strata. For the trial effect θ, the per-observation bound is " +
        showNumber(totals.trialBound) + " if the means can vary separately, and " +
        showNumber(totals.tceBound) + " if equality is imposed along every path. The gap is zero at this law, yet its unrestricted bound is " +
        showNumber(totals.gapBound) + ". Under the common-mean restriction the gap stays zero on every path."
      : "At the selected law, the control means differ in at least one stratum. Set δ to zero to compare the unrestricted and common-mean models at the same law.";
  }
  for (const el of document.querySelectorAll("input, select")) el.addEventListener("input", render);
  get("equal-means").addEventListener("click", () => { get("delta").value = 0; render(); });
  get("reset").addEventListener("click", () => { get("settings").reset(); get("directions").reset(); get("efficiency").reset(); render(); });
  render();
}
if (typeof module !== "undefined" && module.exports) module.exports = { inspect, model, targets };
if (typeof document !== "undefined") startLab();
