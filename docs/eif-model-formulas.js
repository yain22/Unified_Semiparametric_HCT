/* Equation details for the design selector. MathML renders locally in the browser. */
(function () {
  'use strict';

  const id = s => '<mi>' + s + '</mi>';
  const num = s => '<mn>' + s + '</mn>';
  const op = s => '<mo>' + s + '</mo>';
  const roman = s => '<mi mathvariant="normal">' + s + '</mi>';
  const row = s => '<mrow>' + s + '</mrow>';
  const sub = (s, t) => '<msub>' + row(s) + row(t) + '</msub>';
  const sup = (s, t) => '<msup>' + row(s) + row(t) + '</msup>';
  const frac = (s, t) => '<mfrac>' + row(s) + row(t) + '</mfrac>';
  const wrap = (s, a = '(', b = ')') => row(op(a) + s + op(b));
  const braces = s => wrap(s, '{', '}');
  const X = id('X'), O = id('O'), S = id('S'), A = id('A'), Y = id('Y');
  const zero = num(0), one = num(1), q = id('q'), theta = id('θ');
  const eq = op('='), plus = op('+'), minus = op('−');
  const diff = (s, t) => s + minus + t;
  const complement = s => diff(one, s);
  const given = s => op('∣') + s;
  const term = (s, t) => id(s) + wrap(t);
  const atX = s => s + wrap(X);
  const idx = g => g === '10' ? one + op(',') + zero : num(g);
  const I = g => sub(id('I'), idx(g));
  const m = g => atX(sub(id('m'), idx(g)));
  const v = g => atX(sub(id('v'), idx(g)));
  const res = g => sub(id('ε'), idx(g)) + wrap(O);
  const squaredResidual = g => braces(diff(frac(sup(res(g), num(2)), v(g)), one));
  const p = term('p', X), e = term('e', X), r = term('r', X), d = term('d', X);
  const l = term('λ', X), delta = term('δ', X), k = term('c', X);
  const t = atX(sub(id('t'), id('λ')));
  const ttilde = '<mover>' + theta + '<mo>~</mo></mover>';
  const target = sub(ttilde, id('λ'));
  const gap = sub(id('ξ'), id('λ') + op(',') + id('δ'));
  const pooled = atX(sub(id('m'), roman('P')));
  const expected = sub(roman('E'), roman('R'));
  const phi = g => sup(sub(id('φ'), roman(g)), roman('eff')) + wrap(O);
  const hct = phi('U'), trial = phi('trial,H'), common = phi('HCT,trans');
  const math = s => '<math xmlns="http://www.w3.org/1998/Math/MathML">' + s + '</math>';
  const display = (label, lines) => {
    const cells = lines.map(([a,b]) => '<mtr><mtd>' + a + '</mtd><mtd>' + b + '</mtd></mtr>').join('');
    return '<div class="formula-block" role="region" tabindex="0" aria-label="' + label + '">' +
      '<math xmlns="http://www.w3.org/1998/Math/MathML" display="block"><mtable columnalign="right left" ' +
      'columnspacing=".18em" rowspacing=".35em">' + cells + '</mtable></math></div>';
  };
  const para = s => '<p>' + s + '</p>';
  const disclosure = (label, content) => '<details><summary>' + label + '</summary>' + content + '</details>';
  const subscript = (s, index) => sub(s, roman(index));
  const trialResidual = frac(I(1), q + e) + res(1);
  const internalResidual = frac(I('10') + l, q + braces(complement(e))) + res('10');
  const externalResidual = frac(I(0) + p + braces(complement(l)), q + braces(complement(p))) + res(0);
  const cdef = delta + l + braces(complement(l));
  const eg = sub(id('ε'), id('g')) + wrap(O);
  const mg = atX(sub(id('m'), id('g')));
  const vg = atX(sub(id('v'), id('g')));
  const groupCondition = given(X + op(',') + sub(id('I'), id('g')) + eq + one);
  const groupDefinition =
    display('Record and group definitions', [
      [O, eq + wrap(X + op(',') + S + op(',') + A + op(',') + Y)],
      [I(1), eq + S + A],
      [I('10'), eq + S + wrap(complement(A))],
      [I(0), eq + complement(S)],
      [p, eq + roman('Pr') + wrap(S + eq + one + given(X))],
      [q, eq + roman('Pr') + wrap(S + eq + one)],
      [e, eq + roman('Pr') + wrap(A + eq + one + given(X + op(',') + S + eq + one))],
      [mg, eq + roman('E') + wrap(Y + groupCondition)],
      [vg, eq + roman('Var') + wrap(Y + groupCondition)],
      [eg, eq + diff(Y, mg)]
    ]);
  const weightDefinition =
    display('Precision share and borrowing target', [
      [r, eq + frac(v('10'), v(0))],
      [d, eq + p + braces(complement(e)) + plus + braces(complement(p)) + r],
      [l, eq + frac(p + braces(complement(e)), d)],
      [delta, eq + diff(m('10'), m(0))],
      [k, eq + cdef],
      [t, eq + m(1) + minus + l + m('10') + minus + braces(complement(l)) + m(0)],
      [target, eq + expected + braces(t)],
      [gap, eq + diff(target, theta) + eq + expected + braces(braces(complement(l)) + delta)]
    ]);
  const hctNote = disclosure('Notation and overlap', para(
    'The trial is indicated by S = 1. Groups 1, (1, 0), and 0 denote trial treatment, trial control, and external control. ' +
    'Treatment allocation is known. The expectation E<sub>R</sub> averages over trial covariates.' ) +
    groupDefinition + weightDefinition + para(
    'The precision rule is defined on shared support with positive control variances. At trial-only profiles use the internal share of one and evaluate only the trial contributions. At external-only profiles the EIF is zero. The variance bound assumes square-integrable contributions.'));
  const full = para('The target is the treated mean minus the precision-weighted control mean, standardized to trial covariates. Write ' +
    math(k + eq + cdef) + ' and ' + math(eg + eq + diff(Y, mg)) + '.') +
    display('Unrestricted hybrid trial efficient influence function', [
      [hct, eq + frac(S, q) + braces(diff(t, target))],
      ['', plus + trialResidual],
      ['', minus + internalResidual],
      ['', minus + externalResidual],
      ['', minus + frac(k, q) + frac(diff(S, p), complement(p))],
      ['', plus + frac(k, q) + frac(I('10'), complement(e)) + squaredResidual('10')],
      ['', minus + frac(k, q) + frac(I(0) + p, complement(p)) + squaredResidual(0)]
    ]) +
    para('The first four terms reflect trial standardization and the three outcome means. The last three arise because the precision share changes with source membership and the internal and external control variances. With an entirely prespecified share, those three derivatives vanish for that fixed-share target.') +
    hctNote +
    disclosure('The causal-gap and its uncertainty',
      para('Randomization identifies the trial effect θ under the trial identification assumptions. The same hybrid records identify both targets, so the gradient of their difference is the difference of their gradients.') +
      display('Causal-gap and its gradient under the HCT law', [
        [gap, eq + diff(target, theta)],
        [phi('ξ'), eq + diff(hct, trial)],
        [trial, eq + frac(S, q) + braces(m(1) + minus + m('10') + minus + theta)],
        ['', plus + trialResidual + minus + frac(I('10'), q + braces(complement(e))) + res('10')]
      ]) +
      para('Pair the two gradients for each participant when estimating the gap variance. Their shared trial data induce covariance.'));

  const transport = para('Retain the three groups and impose equal conditional control means on shared support, ' +
    math(m('10') + eq + m(0) + eq + pooled) + '. The control variances may differ.') +
    display('Transportable hybrid trial efficient influence function', [
      [common, eq + frac(S, q) + braces(m(1) + minus + pooled + minus + theta)],
      ['', plus + trialResidual],
      ['', minus + frac(p + braces(I('10') + plus + I(0) + r), q + d) + braces(diff(Y, pooled))]
    ]) +
    para('The discrepancy and its multiplier c(X) are zero, so all three share-derivative terms vanish. Both control residuals become ' +
      math(diff(Y, pooled)) + '. The two precision identities give their combined coefficient.') +
    display('Combining control residuals under transportability', [
      [frac(l, complement(e)), eq + frac(p,d)],
      [frac(p + braces(complement(l)), complement(p)), eq + frac(p + r,d)]
    ]) +
    disclosure('Which efficiency bound applies?',
      para('Under the common-mean restriction the borrowing target equals the trial effect throughout the model and the gap is identically zero. At a law where the two control means happen to agree, the unrestricted borrowing EIF has the displayed value. If departures from equality are allowed, the trial-effect EIF still uses only trial outcomes, and the unrestricted gap gradient can be nonzero.') +
      para('For the same hybrid sample, the restricted trial-effect gradient is the projection of the unrestricted trial-effect gradient onto the common-mean tangent space. Its variance cannot exceed the original bound. This comparison holds the sampling law and the number of records fixed.')) +
    hctNote;

  const rct = para('Removing external records leaves a randomized trial. The internal share is fixed at one and the borrowing target equals the trial effect. Each record now belongs to the trial, so S = q = 1.') +
    display('Randomized controlled trial efficient influence function', [
      [phi('RCT'), eq + m(1) + minus + m('10') + minus + theta],
      ['', plus + frac(A,e) + res(1)],
      ['', minus + frac(complement(A),complement(e)) + res('10')]
    ]) +
    para('The external residual and all three share derivatives drop out. The two remaining residuals correct the treatment and control means with their known randomization probabilities.') +
    disclosure('Using only trial outcomes in a hybrid sample',
      para('If all hybrid records are retained but the internal share is fixed at one, the external records contribute zero. The trial contribution still has to be normalized for the hybrid mixture. Conditioning the sample on trial membership gives the pure RCT expression above.') +
      display('Trial-only gradient normalized for HCT mixture sampling', [
        [trial, eq + frac(S,q) + braces(m(1) + minus + m('10') + minus + theta)],
        ['', plus + trialResidual + minus + frac(I('10'), q + braces(complement(e))) + res('10')]
      ]) +
      para('Here S indicates trial membership, q is its mixture probability, and the two indicators select randomized treatment and control. The EIF requires both randomized arms on trial support.'));

  const pE = atX(subscript(id('p'), 'E')), qE = subscript(q, 'E');
  const mE = g => atX(sub(id('m'), idx(g) + op(',') + roman('E')));
  const thetaE = subscript(ttilde,'E'), exE = subscript(roman('E'),'E');
  const Eres = braces(diff(Y,mE(0)));
  const ect = para('All trial participants receive treatment, and all observed controls come from the external source. Under this two-group law, the estimand is the trial-standardized external comparison.') +
    display('Externally controlled trial target and efficient influence function', [
      [thetaE, eq + exE + braces(mE(1) + minus + mE(0) + given(S + eq + one))],
      [phi('ECT'), eq + frac(S,qE) + braces(Y + minus + mE(0) + minus + thetaE)],
      ['', minus + frac(wrap(complement(S)) + pE, qE + braces(complement(pE))) + Eres]
    ]) +
    para('The internal share is fixed at zero. Under the ECT law the treated residual and trial standardization combine because all trial participants are treated. The external residual is reweighted by the ECT source odds. With no internal control arm, the data identify the external comparison. A causal interpretation requires an assumption linking the unobserved trial control mean to the external mean.') +
    display('Combining the two treated terms', [
      [braces(mE(1) + minus + mE(0) + minus + thetaE) + plus + braces(diff(Y,mE(1))), eq + Y + minus + mE(0) + minus + thetaE]
    ]) +
    disclosure('Recomputing the observation law',
      para('Here ' + math(pE + eq + subscript(roman('Pr'),'E') + wrap(S + eq + one + given(X))) + ' and ' +
      math(qE + eq + subscript(roman('Pr'),'E') + wrap(S + eq + one)) +
      '. If the ECT is formed by deleting internal controls from an HCT, condition the original observation law on retaining a treated trial or external record. This gives') +
      display('Two-source probabilities after deleting internal controls', [
        [pE, eq + frac(p + e, complement(p) + plus + p + e)],
        [qE, eq + frac(q + expected + braces(e), complement(q) + plus + q + expected + braces(e))],
        [atX(sub(id('f'), roman('R,E'))), eq + frac(e + atX(sub(id('f'), roman('R'))), expected + braces(e))]
      ]) + para('The last line gives the retained trial covariate density. Constant allocation preserves the original trial covariate distribution. With covariate-dependent allocation, retaining the original trial target requires restandardization. In a separately planned ECT, use its actual two-source law. External controls must cover the target trial support, and the displayed terms must be square-integrable.'));

  const panels = Object.freeze({full, transport, rct, ect});
  if (typeof window !== 'undefined') window.EifFormulaPanels = panels;
  if (typeof module !== 'undefined' && module.exports) module.exports = panels;
}());
