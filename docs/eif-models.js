"use strict";

(() => {
  const by = id => document.getElementById(id);
  const designs = {
    unrestricted: { title: "Unrestricted hybrid trial", formula: "full", n: 300 },
    transportable: { title: "Hybrid trial under control-mean transportability", formula: "transport", n: 300 },
    rct: { title: "Randomized controlled trial", formula: "rct", n: 150 },
    ect: { title: "Externally controlled single-arm trial", formula: "ect", n: 250 }
  };
  const labels = { treated: "Trial treated", internal: "Internal control", external: "External control" };
  const names = { standardization: "Trial standardization", treated: "Treated outcome residual",
    internal: "Internal control residual", external: "External control residual", source: "Source probability",
    internalVariance: "Internal control variance", externalVariance: "External control variance" };
  const groups = { unrestricted: Object.keys(labels), transportable: Object.keys(labels),
    rct: ["treated", "internal"], ect: ["treated", "external"] };

  const fmt = n => n === null ? "Unidentified" :
    (Math.abs(n) < 5e-9 ? 0 : n).toLocaleString("en-US", { minimumFractionDigits: 3, maximumFractionDigits: 6 });

  // Group numerical terms as they are written in each displayed formula.
  function termsFor(model, terms) {
    const a = Object.fromEntries(terms.map(t => [t.key, t.value]));
    if (model === "ect") return [
      ["Treated trial term (including standardization)", a.standardization + a.treated],
      [names.external, a.external]
    ];
    const rows = [[names.standardization, a.standardization], [names.treated, a.treated]];
    if (model === "rct") return rows.concat([[names.internal, a.internal]]);
    if (model === "transportable") {
      return rows.concat([["Combined control residual", a.internal + a.external]]);
    }
    return terms.map(t => [names[t.key], t.value]);
  }

  const notes = {
    unrestricted: "The two control means and variances can differ. The final three terms account for changes in the precision share.",
    transportable: "The means are equal by restriction throughout this model. Both groups contribute to the common control mean and can have different variances.",
    rct: "Only trial records are sampled, with a two-thirds probability of treatment. The external mean and variance play no role.",
    ect: "The external shift specifies an observed external mean for the example. The ECT has no concurrent control group from which to estimate a trial–external difference."
  };
  const structural = {
    unrestricted: "All seven components belong to this model. Individual terms can be zero for this record or these settings.",
    transportable: "The common-mean formula combines the two control residuals. Its three share-derivative terms vanish throughout this model.",
    rct: "The RCT formula has three terms. The external control and share-derivative terms are absent throughout this model.",
    ect: "The ECT formula combines trial standardization and the treated residual in its first term. The internal control and share-derivative terms are absent throughout this model."
  };
  let lastModel;
  // Store the current free setting while the common-mean restriction fixes it at zero.
  let savedDelta;
  let traceModel = null;
  let frameIndex = 0;
  let playTimer = null;

  function stopPlayback() {
    clearTimeout(playTimer);
    playTimer = null;
    const atEnd = traceModel && frameIndex === EifReductionSteps[traceModel].length - 1;
    by("reduction-play").textContent = atEnd ? "Replay" : "Play";
  }

  function showFrame(announce = true) {
    window.EifNotationHelp.close();
    const frames = EifReductionSteps[traceModel];
    const frame = frames[frameIndex];
    const position = `Step ${frameIndex + 1} of ${frames.length}`;
    by("reduction-count").textContent = position;
    by("reduction-title").textContent = frame.title;
    by("reduction-explanation").textContent = frame.explanation;
    by("reduction-equation").innerHTML = frame.html;
    by("reduction-key").querySelector(".term-zero-key").hidden = !by("reduction-equation").querySelector(".reduction-removed");
    const reducedMotion = window.matchMedia("(prefers-reduced-motion: reduce)").matches;
    if (!reducedMotion) {
      by("reduction-equation").animate(
        [{ opacity: .25, transform: "translateY(4px)" }, { opacity: 1, transform: "translateY(0)" }],
        { duration: 280, easing: "ease-out" }
      );
    }
    by("reduction-prev").disabled = frameIndex === 0;
    by("reduction-next").disabled = frameIndex === frames.length - 1;
    if (announce) by("reduction-announcement").textContent = position + ". " + frame.title;
    if (playTimer === null || by("reduction-next").disabled) stopPlayback();
    markWideEquations();
  }

  function openTrace(model) {
    stopPlayback();
    traceModel = model;
    frameIndex = 0;
    by("reduction-player").hidden = false;
    by("formula").hidden = true;
    by("trace-toggle").setAttribute("aria-expanded", "true");
    by("trace-toggle").textContent = "Return to model equation";
    showFrame();
  }

  function closeTrace() {
    stopPlayback();
    traceModel = null;
    by("reduction-player").hidden = true;
    by("formula").hidden = false;
    by("trace-toggle").setAttribute("aria-expanded", "false");
    by("trace-toggle").textContent = "Trace this reduction";
    by("reduction-announcement").textContent = "";
  }

  function scheduleNextFrame() {
    by("reduction-play").textContent = "Pause";
    playTimer = setTimeout(() => {
      frameIndex += 1;
      showFrame();
      if (playTimer !== null) scheduleNextFrame();
    }, 6500);
  }

  function refresh() {
    const model = document.querySelector('input[name="model"]:checked').value;
    const design = designs[model];
    const isTrial = model === "rct", common = model === "transportable", isEct = model === "ect";
    if (lastModel !== model) {
      window.EifNotationHelp.close();
      const old = by("group").value;
      by("group").replaceChildren(...groups[model].map(g => new Option(labels[g], g)));
      by("group").value = groups[model].includes(old) ? old : groups[model][0];
      if (common) {
        savedDelta = by("delta").value;
        by("delta").value = "0";
      } else if (lastModel === "transportable") {
        by("delta").value = savedDelta;
      }
      by("model-name").textContent = design.title;
      by("formula").innerHTML = EifFormulaPanels[design.formula];
      EifRecordTable.render(model);
      if (fittedSelected()) EifEstimatorView.render(model);
      const choosingFromFull = model === "unrestricted";
      by("reduction-path").hidden = !choosingFromFull;
      document.querySelector('label[for="reduction-path"]').hidden = !choosingFromFull;
      if (!choosingFromFull) by("reduction-path").value = model;
      if (traceModel !== null) {
        if (choosingFromFull) closeTrace();
        else openTrace(model);
      }
      lastModel = model;
    }
    by("delta").disabled = isTrial || common;
    by("pattern").disabled = isTrial || common;
    by("v0").disabled = isTrial;
    by("delta").closest("label").hidden = isTrial;
    by("pattern").closest("label").hidden = isTrial || common;
    by("v0").closest("label").hidden = isTrial;
    by("delta-label").textContent = isEct ? "External mean shift below 0.4X" : "Control-mean difference";
    by("law-note").textContent = notes[model];
    const x = Number(by("x").value);
    const r = EifModels.evaluate(model, {
      delta: Number(by("delta").value), v0: Number(by("v0").value),
      pattern: by("pattern").value, overlap: by("overlap").value
    }, { x, group: by("group").value, z: Number(by("z").value) });

    by("target-label").textContent = isTrial || common ? "Trial treatment effect" :
      isEct ? "External comparison target" : "Hybrid borrowing target";
    by("target").textContent = fmt(r.target);
    by("gap").textContent = fmt(r.gap);
    by("score").textContent = fmt(r.phi);
    by("bound").textContent = fmt(r.bound);
    const counts = isTrial ? "100 treated and 50 internal controls" :
      isEct ? "100 treated and 150 external controls" : "100 treated, 50 internal controls, and 150 external controls";
    by("plan-detail").textContent = `For N = ${design.n} independent records (${counts} in expectation), ` +
      `the variance bound B/N is ${fmt(r.bound / design.n)} and its square root, the efficient standard error, is ${fmt(Math.sqrt(r.bound / design.n))}.`;
    by("record-detail").textContent = `At X = ${x}, the trial source probability is ${fmt(r.p)} ` +
      `and the internal control share is ${fmt(r.lambda)}. The selected outcome is ${fmt(r.y)} and its residual is ${fmt(r.residual)}.`;
    by("terms").replaceChildren(...termsFor(model, r.terms).map(([name, value]) => {
      const row = document.createElement("tr");
      const term = document.createElement("th");
      term.scope = "row";
      term.textContent = name;
      row.append(term);
      row.insertCell().textContent = fmt(value);
      return row;
    }));
    let zeroNote = structural[model];
    if (model === "unrestricted" && Number(by("delta").value) === 0) {
      zeroNote += " At these equal means the three share-derivative contributions are zero; unequal-mean paths remain available in the unrestricted model.";
    } else if (model !== "unrestricted") {
      zeroNote += " A zero shown in the table is the contribution for this record at the chosen settings.";
    }
    by("zero-note").textContent = zeroNote;
  }
  function markWideEquations() {
    for (const block of document.querySelectorAll(".formula-block")) {
      const wide = block.scrollWidth > block.clientWidth + 1;
      block.classList.toggle("wide-equation", wide);
      block.tabIndex = wide ? 0 : -1;
    }
  }
  document.querySelectorAll(".experiment select, #models input").forEach(el => el.addEventListener("change", () => {
    refresh();
    markWideEquations();
  }));
  window.addEventListener("resize", markWideEquations);
  by("formula").addEventListener("toggle", markWideEquations, true);
  by("estimator-formula").addEventListener("toggle", markWideEquations, true);
  function fittedSelected() {
    return document.querySelector('input[name="equation-view"]:checked').value === "estimator";
  }
  function switchView() {
    const fitted = fittedSelected();
    window.EifNotationHelp.close();
    closeTrace();
    for (const id of ["trace-toolbar", "formula", "population-records", "population-example", "population-r"]) {
      by(id).hidden = fitted;
    }
    by("fitted-estimator").hidden = !fitted;
    if (fitted) EifEstimatorView.render(document.querySelector('input[name="model"]:checked').value);
    markWideEquations();
  }
  by("estimators").addEventListener("change", switchView);
  by("adjustment-control").addEventListener("change", () => {
    window.EifNotationHelp.close();
    EifEstimatorView.updateTrial();
    markWideEquations();
  });
  by("trace-toggle").addEventListener("click", () => {
    if (traceModel !== null) {
      window.EifNotationHelp.close();
      closeTrace();
      markWideEquations();
    } else {
      const model = by("reduction-path").value;
      document.querySelector(`#models input[value="${model}"]`).checked = true;
      refresh();
      openTrace(model);
    }
  });
  function moveFrame(offset) {
    stopPlayback();
    frameIndex += offset;
    showFrame();
  }
  by("reduction-prev").addEventListener("click", () => moveFrame(-1));
  by("reduction-next").addEventListener("click", () => moveFrame(1));
  by("reduction-play").addEventListener("click", () => {
    if (playTimer !== null) {
      stopPlayback();
      return;
    }
    if (frameIndex === EifReductionSteps[traceModel].length - 1) {
      frameIndex = 0;
      showFrame();
    }
    scheduleNextFrame();
  });
  for (const eventName of ["focusin", "pointerover"]) {
    by("reduction-equation").addEventListener(eventName, event => {
      if (playTimer !== null && event.target.closest("[data-notation]")) stopPlayback();
    });
  }
  document.addEventListener("visibilitychange", () => {
    if (document.hidden) stopPlayback();
  });
  by("trace-toolbar").hidden = false;
  function openLinkedView() {
    if (location.hash === "#rct-estimators") {
      document.querySelector('input[name="model"][value="rct"]').checked = true;
      document.querySelector('input[name="equation-view"][value="estimator"]').checked = true;
      refresh();
    } else if (location.hash === "#estimators") {
      document.querySelector('input[name="equation-view"][value="estimator"]').checked = true;
    } else if (["#record-table", "#models"].includes(location.hash)) {
      document.querySelector('input[name="equation-view"][value="eif"]').checked = true;
    } else return;
    switchView();
    if (location.hash === "#rct-estimators") by("estimator-heading").scrollIntoView();
  }
  refresh();
  openLinkedView();
  window.addEventListener("hashchange", openLinkedView);
  markWideEquations();
})();
