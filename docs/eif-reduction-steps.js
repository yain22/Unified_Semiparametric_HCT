/* Successive equations for the three model reductions. The first and last
   equations come from the same panels as the stationary formula display. */
(() => {
  "use strict";
  const h = window.EifMath, panels = window.EifFormulaPanels;
  const { id, op, roman, sub, sup, frac, wrap, braces, I, m, res, mark,
    S, A, Y, zero, one, q, theta, eq, plus, minus, diff, complement,
    p, e, r, d, l, delta, k, t, pooled, expected, trial, common,
    trialResidual, pE, qE, mE, thetaE, Eres, fR, fRE } = h;

  function markRows(html, states) {
    let index = 0;
    return html.replace(/<mtr\b([^>]*)>/g, (tag, attributes) => {
      const state = states[index++] || "retained";
      if (/\bclass="/.test(attributes)) {
        return '<mtr' + attributes.replace(/\bclass="/, 'class="reduction-' + state + ' ') + '>';
      }
      return '<mtr class="reduction-' + state + '"' + attributes + '>';
    });
  }
  function panelEquation(panel, states) {
    const start = panel.indexOf('<div class="formula-block"');
    const endTag = '</math></div>';
    const end = panel.indexOf(endTag, start);
    if (start < 0 || end < 0) throw new Error("Missing model equation in the formula panel.");
    return markRows(panel.slice(start, end + endTag.length), states);
  }
  function display(label, lines) {
    return markRows(h.display(label, lines.map(([a, b]) => [a, b])), lines.map(line => line[2]));
  }

  const emph = "emphasized", gone = "removed", keep = "retained";
  const prediction = (a, b, center) => braces(a + minus + b + minus + center);
  const epsP = braces(diff(Y, pooled));
  const tP = diff(m(1), pooled);
  const tTrial = diff(m(1), m("10"));
  const indexH0 = zero + op(",") + roman("H");
  const thetaH0 = mark("hybridExternalTarget", sub('<mover>' + id("θ") + '<mo>~</mo></mover>', indexH0));
  const phiH0 = mark("hybridExternalEIF", sup(sub(id("φ"), indexH0), roman("eff")) + wrap(id("O")));
  const t0 = diff(m(1), m(0));
  const tE = diff(mE(1), mE(0));
  const conditionalTrial = sub(roman("P"), roman("R"));
  const meanGroupTerm = (ind, denominator, residual) => frac(ind, denominator) + residual;
  const first = (explanation, highlights) => ({
    title: "Start with the unrestricted hybrid trial",
    explanation, html: panelEquation(panels.full, highlights)
  });
  const last = (panel, title, explanation, highlights) => ({
    title, explanation, html: panelEquation(panels[panel], highlights)
  });

  const transportable = [
    first("All three groups are sampled. The last three terms describe how the target changes when the source probability or a control variance changes the precision share.",
      [keep, keep, emph, emph, emph, emph, emph]),
    {
      title: "Impose a common conditional control mean",
      explanation: "On shared support the two control means now agree and their variances may still differ. Both residuals use the same mean. The discrepancy multiplier is zero, so all three share-derivative terms vanish.",
      html: display("Common control mean and vanishing derivative terms", [
        [m("10"), eq + m(0) + eq + pooled, emph],
        [t, eq + tP, emph],
        [delta, eq + zero + op(",") + k + eq + zero, gone],
        [common, eq + frac(S, q) + prediction(m(1), pooled, theta) + plus + trialResidual, keep],
        ["", minus + frac(I("10") + l, q + braces(complement(e))) + epsP, emph],
        ["", minus + frac(I(0) + p + braces(complement(l)), q + braces(complement(p))) + epsP, emph]
      ])
    },
    {
      title: "Combine the two control residuals",
      explanation: "The internal and external control terms now multiply the same residual. Substituting the precision share gives one coefficient that still depends on the source and on the ratio of the two control variances.",
      html: display("Precision identities and common residual contribution", [
        [frac(l, complement(e)), eq + frac(p, d), emph],
        [frac(p + braces(complement(l)), complement(p)), eq + frac(p + r, d), emph],
        ["", minus + frac(I("10") + p, q + d) + epsP + minus + frac(I(0) + p + r, q + d) + epsP, keep],
        ["", eq + minus + frac(p + braces(I("10") + plus + I(0) + r), q + d) + epsP, emph]
      ])
    },
    last("transport", "EIF in the common-mean hybrid model",
      "The target is the trial effect throughout this restricted model. Its causal gap is identically zero. At an equal-mean law the unrestricted borrowing EIF takes this same value, although the unrestricted gap can change along paths that move the two means apart.",
      [keep, keep, emph])
  ];

  const rct = [
    first("Follow the trial covariate, treated outcome, and internal control terms. A conventional randomized trial obtains its comparison from these two trial arms.",
      [emph, emph, emph, keep, keep, keep, keep]),
    {
      title: "Use the trial-only contrast in the hybrid sample",
      explanation: "First choose a fixed internal share of one. The target becomes the trial effect. The external residual has coefficient zero, and a fixed share has no derivative. The factor S/q still standardizes over trial participants within the hybrid mixture.",
      html: display("Trial-only target in a hybrid mixture", [
        [l, op("≡") + one + op(",") + t + eq + tTrial, emph],
        [braces(complement(l)), eq + zero, gone],
        [trial, eq + frac(S, q) + prediction(m(1), m("10"), theta), keep],
        ["", plus + trialResidual + minus + meanGroupTerm(I("10"), q + braces(complement(e)), res("10")), emph]
      ])
    },
    {
      title: "Condition the observation law on trial membership",
      explanation: "Retain only trial records. Every record now has S = 1, the trial fraction is one, and the two arm indicators reduce to A and 1 − A. The trial’s randomization probability at X is unchanged.",
      html: display("From mixture normalization to trial randomization weights", [
        [conditionalTrial, eq + roman("P") + wrap(op("·") + h.given(S + eq + one)), emph],
        [S, eq + q + eq + one, emph],
        [frac(S,q), eq + one, keep],
        [frac(I(1),q + e), eq + frac(A,e), emph],
        [frac(I("10"),q + braces(complement(e))), eq + frac(complement(A),complement(e)), emph]
      ])
    },
    last("rct", "EIF in the randomized trial",
      "All expectations and variances now use the trial-only law. The external cohort does not enter the target or its variance, and the causal gap is zero by design.",
      [keep, emph, emph])
  ];

  const ect = [
    first("An externally controlled single-arm trial retains treated trial records and external controls. Follow both the target and the law supplying its source weights.",
      [emph, emph, keep, emph, keep, keep, keep]),
    {
      title: "Define the external comparison in the hybrid population",
      explanation: "First choose a fixed internal share of zero. The trial-standardized contrast compares the treated mean with the external control mean. Fixed-share derivative terms and the internal residual vanish. At this stage p(X), q, and the randomization probability still refer to the hybrid law.",
      html: display("External comparison while all hybrid records are sampled", [
        [l, op("≡") + zero + op(",") + thetaH0 + eq + expected + braces(t0), emph],
        [phiH0, eq + frac(S,q) + prediction(m(1),m(0),thetaH0), keep],
        ["", plus + trialResidual, emph],
        ["", minus + meanGroupTerm(I(0) + p, q + braces(complement(p)), res(0)), emph]
      ])
    },
    {
      title: "Recompute source weights and combine the trial terms",
      explanation: "Retain treated trial records and external controls. Conditioning on retention gives the source weights and trial covariate density below. All retained trial participants receive treatment, so the treated residual combines with trial standardization. A separately planned ECT uses its own sampling law.",
      html: display("Source probabilities and trial density after deleting internal controls", [
        [pE, eq + frac(p + e, complement(p) + plus + p + e), emph],
        [qE, eq + frac(q + expected + braces(e), complement(q) + plus + q + expected + braces(e)), emph],
        [fRE, eq + frac(e + fR, expected + braces(e)), emph]
      ]) + display("Combining the two trial contributions and reweighting the external residual", [
        [frac(S,qE) + braces(tE + minus + thetaE + plus + wrap(diff(Y,mE(1)))), eq + frac(S,qE) + braces(Y + minus + mE(0) + minus + thetaE), emph],
        ["", minus + frac(wrap(complement(S)) + pE, qE + braces(complement(pE))) + Eres, emph]
      ])
    },
    last("ect", "EIF in the externally controlled trial",
      "With external support throughout the trial covariate distribution, the two-source data identify this standardized external contrast. The trial effect and its causal gap still require information about the trial control outcome. When forming an ECT by deleting internal controls, covariate-dependent randomization can change the retained trial population.",
      [emph, emph, emph])
  ];

  window.EifReductionSteps = Object.freeze({ transportable, rct, ect });
})();
