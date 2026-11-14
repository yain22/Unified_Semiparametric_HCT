/* Six hypothetical records from Figure S1. Keep the guide usable offline. */
(() => {
  "use strict";
  const records = [
    ["P01", 1, 1, .20, 0, 8.4], ["P02", 1, 1, .70, 1, 10.1],
    ["P03", 1, 0, .30, 1, 7.2], ["P04", 1, 0, .80, 0, 8.5],
    ["P05", 0, 0, .25, 0, 6.9], ["P06", 0, 0, .65, 1, 8.0]
  ];
  const captions = {
    unrestricted: "Unrestricted hybrid trial · all six records",
    transportable: "Common-mean hybrid trial · all six records",
    rct: "Randomized trial · the four trial participants",
    ect: "Externally controlled trial · two treated participants and two external controls"
  };
  const interpretation = {
    unrestricted: "Both control sources are observed. At the same baseline covariate values, their mean outcomes may differ.",
    transportable: "The same six rows are observed. The model assumes equal conditional means across control sources at shared covariate values. Individual control outcomes may differ.",
    rct: "Only trial participants are observed, with both treated and control outcomes represented.",
    ect: "All trial participants receive treatment. The external cohort supplies the observed controls. Trial participants’ outcomes under control are unobserved."
  };
  function render(model, mount = document.getElementById("record-table")) {
    if (!mount || !Object.hasOwn(captions, model)) return;
    const table = document.createElement("table");
    table.createCaption().textContent = captions[model];
    const head = table.createTHead().insertRow();
    for (const title of ["ID", "S", "A", "X₁", "X₂", "Y", "Observed group"]) {
      const th = document.createElement("th");
      th.scope = "col";
      th.textContent = title;
      head.append(th);
    }
    const body = table.createTBody();
    for (const [id, s, a, x1, x2, y] of records) {
      if ((model === "rct" && !s) || (model === "ect" && s && !a)) continue;
      const tr = body.insertRow();
      const label = document.createElement("th");
      label.scope = "row";
      label.textContent = id;
      tr.append(label);
      const group = a ? "Trial treated (Y₁)" : s ? "Internal control (Y₁,₀)" : "External control (Y₀)";
      for (const value of [s, a, x1.toFixed(2), x2, y.toFixed(1), group]) {
        tr.insertCell().textContent = value;
      }
    }
    const viewport = document.createElement("div");
    viewport.className = "observed-records-scroll";
    viewport.tabIndex = 0;
    viewport.setAttribute("role", "region");
    viewport.setAttribute("aria-label", "Observed records table");
    viewport.append(table);
    const note = document.createElement("p");
    note.textContent = interpretation[model];
    const key = document.createElement("p");
    key.textContent = "S identifies trial membership and A identifies treatment. The observed group indicators are I₁ = SA (trial treated), I₁,₀ = S(1 − A) (internal control), and I₀ = 1 − S (external control). Each record has one observed outcome Y. The subscripts on Y name its group.";
    mount.replaceChildren(viewport, note, key);
  }
  window.EifRecordTable = { render };
})();
