"use strict";

const EfficiencyBounds = (() => {
  const defaults = { nT: 100, nI: 50, nE: 150, vT: 1, vI: 1, vE: 1, delta: 0, plan: "retain" };
  const models = {
    unrestricted: { label: "Unrestricted HCT", color: "#176cb0", shape: "circle" },
    transportable: { label: "Common-mean HCT", color: "#076b60", shape: "diamond" },
    ect: { label: "ECT", color: "#c75a18", shape: "open-circle" },
    rct: { label: "RCT", color: "#50575d", shape: "square" }
  };

  // Bounds are on the estimator variance scale (E[phi²]/N for each design).
  // Counts specify probabilities in an iid mixture. Randomization is known,
  // source membership is random, and the outcome variances are unknown.
  function evaluate(settings = {}) {
    const { nT, nI, nE, vT, vI, vE, delta, plan } = { ...defaults, ...settings };
    if (![nT, nI, nE, vT, vI, vE].every(x => Number.isFinite(x) && x > 0))
      throw new RangeError("Enter a positive expected count and variance for each group.");
    if (!Number.isFinite(delta)) throw new RangeError("Enter a finite control-mean difference.");
    if (!["retain", "equal-total"].includes(plan)) throw new RangeError("Select a sampling plan.");

    const precision = nI / vI + nE / vE;
    if (!Number.isFinite(precision) || precision <= 0)
      throw new RangeError("These settings exceed the numerical range of the calculator.");
    const share = (nI / vI) / precision;
    const externalShare = (nE / vE) / precision;
    const base = vT / nT + 1 / precision;
    const sensitivitySquared = (delta * share * externalShare) ** 2;
    // The source score estimates trial vs external membership. Conditional on
    // trial membership, treatment allocation is fixed by the known probability.
    const extraSource = sensitivitySquared * (1 / (nT + nI) + 1 / nE);
    // Centered squared Gaussian residuals have variance 2 on the standardized scale.
    const extraInternal = 2 * sensitivitySquared / nI;
    const extraExternal = 2 * sensitivitySquared / nE;
    const rctControls = plan === "retain" ? nI : nI + nE;
    const ectControls = plan === "retain" ? nE : nI + nE;
    const common = delta === 0;
    const rows = [
      { model: "unrestricted", nT, nI, nE, estimand: "Precision-weighted contrast",
        target: 1 + externalShare * delta, bound: base + extraSource + extraInternal + extraExternal },
      { model: "transportable", nT, nI, nE, estimand: "Trial treatment effect",
        target: common ? 1 : null, bound: common ? base : null },
      { model: "ect", nT, nI: 0, nE: ectControls, estimand: "External comparison",
        target: 1 + delta, bound: vT / nT + vE / ectControls },
      { model: "rct", nT, nI: rctControls, nE: 0, estimand: "Trial treatment effect",
        target: 1, bound: vT / nT + vI / rctControls }
    ];
    for (const row of rows) {
      row.n = row.nT + row.nI + row.nE;
      row.se = row.bound === null ? null : Math.sqrt(row.bound);
      if (!Number.isFinite(row.n) || (row.bound !== null &&
          (!Number.isFinite(row.target) || !Number.isFinite(row.bound) || row.bound <= 0))) {
        throw new RangeError("These settings exceed the numerical range of the calculator.");
      }
    }
    return { share, externalShare, base, extraSource, extraInternal, extraExternal, rows };
  }

  // Show enough precision to distinguish nearby bounds. Raw numeric values
  // are also attached to each table cell for machine-readable comparisons.
  function format(value, digits = 9) {
    const rounded = value.toPrecision(digits);
    const number = Number(rounded);
    return (Number.isFinite(number) ? String(number) : rounded).replace(/^-/, "−");
  }

  function svgElement(tag, attributes = {}, content) {
    const node = document.createElementNS("http://www.w3.org/2000/svg", tag);
    for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, value);
    if (content !== undefined) node.textContent = content;
    return node;
  }

  function makeChart(host) {
    const svg = svgElement("svg", { role: "img", width: "100%",
      "aria-labelledby": "bounds-chart-title bounds-chart-description",
      "font-family": "system-ui, sans-serif", fill: "#253f43" });
    svg.append(svgElement("title", { id: "bounds-chart-title" }, "Efficiency bounds on one variance scale"));
    const description = svgElement("desc", { id: "bounds-chart-description" });
    const grid = svgElement("g", { "aria-hidden": "true" });
    svg.append(description, grid);
    const lanes = {};
    for (const [index, [model, style]] of Object.entries(models).entries()) {
      const y = 23 + 37 * index;
      const lane = svgElement("g", { "data-model": model });
      const rail = svgElement("line", { y1: y, y2: y, stroke: "#c9d6d3" });
      const label = svgElement("text", { x: 0, y: y + 4.5, "font-size": 13,
        "font-weight": 600 }, style.label);
      const value = svgElement("text", { y: y + 4, "font-size": 12, "text-anchor": "end" });
      const missing = svgElement("text", { y: y + 4, "font-size": 12, fill: "#526567" }, "Inapplicable");
      let marker;
      const open = style.shape === "open-circle";
      if (style.shape === "circle" || open) marker = svgElement("circle", { r: 6 });
      if (style.shape === "square") marker = svgElement("rect", { x: -5.5, y: -5.5, width: 11, height: 11 });
      if (style.shape === "diamond") marker = svgElement("polygon", { points: "0,-7 7,0 0,7 -7,0" });
      marker.setAttribute("class", "bound-marker");
      marker.setAttribute("fill", open ? "#fff" : style.color);
      marker.setAttribute("stroke", open ? style.color : "#fff");
      marker.setAttribute("stroke-width", open ? "2" : "1.2");
      lane.append(rail, label, value, marker);
      svg.append(lane);
      lanes[model] = { lane, rail, value, marker, missing, y };
    }
    host.replaceChildren(svg);
    return function draw(rows) {
      const width = Math.max(240, host.clientWidth || 600);
      const compact = width < 520;
      const left = compact ? 130 : 153;
      const right = width - (compact ? 16 : 104);
      const axisY = 157;
      const intervalCount = compact ? 2 : 4;
      const maxBound = Math.max(...rows.map(row => row.bound ?? 0));
      const rough = maxBound / intervalCount;
      const unit = 10 ** Math.floor(Math.log10(rough));
      const multiplier = [1, 2, 2.5, 5, 10].find(m => m * unit >= rough);
      const roundedEnd = multiplier * unit * intervalCount;
      // Rounding upward can overflow for a bound near the largest float.
      // Very small subnormal values may not admit distinct intermediate ticks.
      const end = Number.isFinite(roundedEnd) && roundedEnd >= maxBound ? roundedEnd : maxBound;
      const x = bound => left + (bound / end) * (right - left);
      svg.setAttribute("viewBox", "0 0 " + width + " 214");
      description.textContent = rows.map(row => models[row.model].label + ", " +
        (row.bound === null ? "inapplicable" : format(row.bound))).join(". ") +
        ". Smaller values are farther left. The table gives the identified target of each design.";
      grid.replaceChildren();
      let previousTick;
      for (let i = 0; i <= intervalCount; ++i) {
        const tick = end * (i / intervalCount);
        if (tick === previousTick) continue;
        previousTick = tick;
        const at = x(tick);
        grid.append(svgElement("line", { x1: at, x2: at, y1: 10, y2: axisY + 4, stroke: "#c9d6d3",
          "stroke-dasharray": "2 3" }));
        grid.append(svgElement("text", { x: at, y: axisY + 20, "font-size": 12,
          "text-anchor": tick === end ? "end" : "middle" }, format(tick, 4)));
      }
      grid.append(svgElement("line", { x1: left, x2: right, y1: axisY, y2: axisY, stroke: "#42595a" }));
      grid.append(svgElement("path", { d: "M " + right + " " + (axisY - 3.5) + " l 6 3.5 l -6 3.5 Z",
        fill: "#42595a" }));
      grid.append(svgElement("text", { x: width / 2, y: axisY + 46, "text-anchor": "middle",
        "font-size": 13 }, "Variance bound  B/N  ·  smaller to the left"));
      for (const row of rows) {
        const { lane, rail, value, marker, missing, y } = lanes[row.model];
        rail.setAttribute("x1", left);
        rail.setAttribute("x2", right);
        value.setAttribute("x", width - 2);
        value.textContent = compact || row.bound === null ? "" : format(row.bound, 7);
        if (row.bound === null) {
          marker.remove();
          rail.setAttribute("visibility", "hidden");
          missing.setAttribute("x", left + 4);
          lane.append(missing);
        } else {
          missing.remove();
          rail.setAttribute("visibility", "visible");
          if (!marker.parentNode) lane.append(marker);
          marker.setAttribute("transform", "translate(" + x(row.bound) + " " + y + ")");
          marker.setAttribute("data-number", row.bound);
        }
      }
    };
  }

  function renderTable(host, rows) {
    const table = document.createElement("table");
    table.createCaption().textContent = "Expected counts, targets and bounds";
    const headers = table.createTHead().insertRow();
    for (const label of ["Model", "Treated / internal / external", "Target", "Variance bound B/N", "Efficient SE"]) {
      const th = document.createElement("th");
      th.scope = "col";
      th.textContent = label;
      headers.append(th);
    }
    const body = table.createTBody();
    for (const row of rows) {
      const tr = body.insertRow();
      tr.dataset.model = row.model;
      const heading = document.createElement("th");
      heading.scope = "row";
      heading.textContent = models[row.model].label;
      tr.append(heading);
      tr.insertCell().textContent = [row.nT, row.nI, row.nE].map(n => format(n)).join(" / ");
      for (const [key, number] of [["target", row.target], ["variance", row.bound], ["se", row.se]]) {
        const cell = tr.insertCell();
        cell.dataset.value = key;
        if (number === null) {
          cell.textContent = "Inapplicable";
        } else {
          cell.dataset.number = String(number);
          cell.textContent = key === "target" ? row.estimand + " = " + format(number) : format(number);
        }
      }
    }
    host.replaceChildren(table);
  }

  function mount(root) {
    const host = root.querySelector("#bound-line");
    const results = root.querySelector("#bound-results");
    const status = root.querySelector("#bound-status");
    const draw = makeChart(host);
    let current;
    function refresh() {
      const settings = {};
      for (const key of Object.keys(defaults)) {
        const field = root.querySelector('[name="' + key + '"]:checked') ||
          root.querySelector('[name="' + key + '"]');
        settings[key] = key === "plan" ? field.value : (field.value.trim() === "" ? NaN : Number(field.value));
        if (field.id) {
          const output = root.querySelector('output[for="' + field.id + '"]');
          if (output) output.value = field.value;
        }
      }
      try {
        current = evaluate(settings);
      } catch (error) {
        current = null;
        host.hidden = true;
        results.replaceChildren();
        status.textContent = error.message;
        return;
      }
      host.hidden = false;
      draw(current.rows);
      renderTable(results, current.rows);
      const share = format(100 * current.share, 5) + "%";
      if (settings.delta === 0) {
        let ordering;
        if (settings.plan === "retain")
          ordering = "Using both control groups gives the hybrid designs a smaller bound than either single-control design.";
        else if (settings.vI === settings.vE)
          ordering = "All four bounds coincide because each design has the same total control count and control variance.";
        else
          ordering = "At the same total control count, the hybrid bound lies between the RCT and ECT bounds.";
        status.textContent = "Both hybrid bounds coincide at " + format(current.base) +
          ". All four contrasts equal 1 at these control means. The hybrid comparison assigns " + share +
          " of the control mean to internal controls. " + ordering;
      } else {
        status.textContent = "The designs now estimate different contrasts. The hybrid target is " +
          format(current.rows.find(row => row.model === "unrestricted").target) + " with " + share +
          " internal weight. The trial effect is 1 and the external comparison is " +
          format(current.rows.find(row => row.model === "ect").target) +
          ". The common-mean model is inapplicable.";
      }
    }
    root.addEventListener("input", refresh);
    root.addEventListener("change", refresh);
    const examples = {
      retain: defaults,
      "equal-total": { ...defaults, plan: "equal-total", vE: 4 },
      "unequal-means": { ...defaults, delta: 1 }
    };
    for (const button of root.querySelectorAll("button[data-bound-example]")) {
      button.addEventListener("click", () => {
        const example = examples[button.dataset.boundExample];
        for (const [key, value] of Object.entries(example)) {
          root.querySelector('[name="' + key + '"]').value = value;
        }
        refresh();
      });
    }
    refresh();
    if (typeof ResizeObserver !== "undefined") {
      let previousWidth = host.clientWidth;
      new ResizeObserver(() => {
        if (current && previousWidth !== host.clientWidth) {
          previousWidth = host.clientWidth;
          draw(current.rows);
        }
      }).observe(host);
    }
    return refresh;
  }

  return { evaluate, mount };
})();

if (typeof module !== "undefined" && module.exports) module.exports = EfficiencyBounds;
if (typeof window !== "undefined") {
  window.EfficiencyBounds = EfficiencyBounds;
  function openLinkedLaw() {
    const id = window.location.hash.slice(1);
    if (!["law-unrestricted", "law-common", "law-rct", "law-ect"].includes(id)) return;
    const details = document.getElementById(id);
    if (details && !details.open) {
      details.open = true;
      details.scrollIntoView();
    }
  }
  const start = () => {
    const root = document.getElementById("efficiency-bounds");
    if (root) EfficiencyBounds.mount(root);
    openLinkedLaw();
    window.addEventListener("hashchange", openLinkedLaw);
  };
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", start, { once: true });
  else start();
}
