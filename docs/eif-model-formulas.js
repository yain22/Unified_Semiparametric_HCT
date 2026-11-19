/* Equation details for the design selector. MathML renders locally in the browser. */
(function () {
  'use strict';

  const entry = (label, definition, role) => Object.freeze({label, definition, role});
  const notation = Object.freeze({
    record: entry('Observed record O', 'The baseline covariates X, trial membership R, treatment A, and observed outcome Y for one participant.', 'The influence function assigns a contribution to each observed record.'),
    covariates: entry('Baseline covariates X', 'Patient characteristics measured before treatment that enter the source and outcome models.', 'The comparison averages over the distribution of these covariates in the target trial.'),
    source: entry('Trial membership R', 'R is one for a trial participant and zero for an external control.', 'It selects the records used to standardize the contrast to the trial population.'),
    treatment: entry('Treatment indicator A', 'A is one for the investigational treatment and zero for the comparator. All external controls have A = 0.', 'In a randomized trial it identifies the assigned arm. In a single-arm trial all trial participants have A = 1.'),
    outcome: entry('Observed outcome Y', 'The measured outcome for this participant under the treatment they received.', 'Subtracting the appropriate group mean from Y gives the outcome residual.'),
    trialFraction: entry('Trial sampling fraction q', 'The marginal probability of trial membership under the hybrid sampling law, Pr(R = 1).', 'Dividing by q expresses trial-population averages in terms of all hybrid records.'),
    trialEffect: entry('Trial treatment effect θ', 'The mean potential outcome under treatment minus that under control in the target trial population.', 'Randomized arms identify this effect under the trial identification assumptions.'),
    indicator1: entry('Treated trial indicator I₁', 'I₁ equals R times A and is one only for a treated trial participant.', 'Selects the treated-outcome residual.'),
    indicator10: entry('Trial control indicator I₁,₀', 'I₁,₀ equals R times (1 − A) and is one only for a randomized trial control.', 'Selects the internal-control outcome and variance contributions.'),
    indicator0: entry('External control indicator I₀', 'I₀ equals 1 − R and is one only for an external control.', 'Selects the external-control outcome and variance contributions.'),
    indicatorG: entry('Group indicator I_g', 'An indicator that the record belongs to group g, one of trial treatment, trial control, or external control.', 'Conditioning on I_g = 1 selects the group whose outcome law is being described.'),
    sourceProbability: entry('Conditional trial probability p(X)', 'The probability of trial membership given the baseline covariates under the hybrid sampling law.', 'It determines both the external source odds and the relative availability of internal controls.'),
    randomization: entry('Randomization probability e(X)', 'The known probability of treatment among randomized trial participants at covariates X.', 'Its inverse and the inverse of 1 − e(X) weight the two randomized arms.'),
    meanG: entry('Group outcome mean m_g(X)', 'The conditional mean of Y given X and membership in group g.', 'It defines the fitted or population mean subtracted from the observed outcome.'),
    mean1: entry('Treated trial mean m₁(X)', 'The mean outcome at X among treated trial participants.', 'It supplies the treated side of the treatment contrast.'),
    mean10: entry('Trial control mean m₁,₀(X)', 'The mean outcome at X among randomized trial controls.', 'It supplies the internal comparator and, on shared support, the internal side of the source discrepancy.'),
    mean0: entry('External control mean m₀(X)', 'The mean outcome at X among external controls.', 'Its contribution to the hybrid comparator is weighted by the external-control share.'),
    varianceG: entry('Group outcome variance v_g(X)', 'The conditional variance of Y given X and membership in group g.', 'The relative variances of the two control groups enter the precision share.'),
    variance1: entry('Treated trial variance v₁(X)', 'The conditional variance of the outcome among treated trial participants at X.', 'It contributes to the sampling variance of the treated-outcome term.'),
    variance10: entry('Trial control variance v₁,₀(X)', 'The conditional variance of the outcome among randomized trial controls at X.', 'It enters the variance ratio and the correction for estimating internal-control precision.'),
    variance0: entry('External control variance v₀(X)', 'The conditional variance of the outcome among external controls at X.', 'It enters the variance ratio and the correction for estimating external-control precision.'),
    residualG: entry('Group outcome residual ε_g(O)', 'The observed outcome minus m_g(X) for a record in group g.', 'Its conditional mean is zero within that group.'),
    residual1: entry('Treated trial residual ε₁(O)', 'The observed outcome Y minus the treated trial mean m₁(X).', 'The treated-outcome term weights this residual by the inverse randomization probability.'),
    residual10: entry('Trial control residual ε₁,₀(O)', 'The observed outcome Y minus the trial control mean m₁,₀(X).', 'This residual enters the internal-control mean term. Its square also enters the precision-share correction.'),
    residual0: entry('External control residual ε₀(O)', 'The observed outcome Y minus the external control mean m₀(X).', 'This residual enters the external-control mean term. Its square also enters the precision-share correction.'),
    varianceRatio: entry('Control variance ratio r(X)', 'The trial control variance v₁,₀(X) divided by the external control variance v₀(X) on shared support.', 'At fixed source availability, a larger ratio gives external controls more weight.'),
    weightDenominator: entry('Precision-share denominator d(X)', 'The sum p(X){1 − e(X)} + {1 − p(X)}r(X).', 'It combines the relative availability and precision of the two control groups.'),
    internalShare: entry('Internal-control share λ(X)', 'The weight assigned to the trial control mean at X. The external control mean receives the remaining weight, 1 − λ(X).', 'The precision rule sets this share to p(X){1 − e(X)}/d(X) on shared support.'),
    discrepancy: entry('Control-mean discrepancy δ(X)', 'The internal control mean minus the external control mean at X, on their shared support.', 'The trial average of this discrepancy weighted by the external share is the causal-gap.'),
    derivativeMultiplier: entry('Share-derivative multiplier c(X)', 'The product δ(X)λ(X){1 − λ(X)} on shared support.', 'It scales all three corrections for estimating the precision share and vanishes when the two control means agree.'),
    conditionalContrast: entry('Conditional borrowing contrast tλ(X)', 'The treated trial mean minus the hybrid control mean at X, where λ(X) is the internal-control share.', 'Averaging this contrast over trial covariates gives the borrowing target.'),
    borrowingTarget: entry('Hybrid borrowing target θ̃λ', 'The average of the conditional borrowing contrast tλ(X) over the trial covariate distribution.', 'This is the statistical target of AIPW-U. It can differ from the trial treatment effect when the control means differ.'),
    hybridExternalTarget: entry('Hybrid target with a fixed external-only comparator', 'The trial mean under treatment minus the external control mean standardized to all trial covariates, with the internal share fixed at zero.', 'This target requires external control support throughout the trial covariate distribution.'),
    hybridExternalEIF: entry('Hybrid EIF with a fixed external-only comparator', 'The efficient influence function for the external-only comparator target under the original hybrid sampling law.', 'All trial records contribute to standardization. The treated outcome residual still uses the trial randomization probability.'),
    causalGap: entry('Causal-gap ξλ,δ', 'The hybrid borrowing target minus the trial treatment effect.', 'It measures how far borrowing shifts the treatment contrast on the outcome scale.'),
    pooledMean: entry('Common control mean mP(X)', 'The mean of Y among all controls given X. Under control-mean transportability both source-specific control means equal it on shared support.', 'The restricted hybrid EIF can use both groups to estimate this one mean.'),
    trialExpectation: entry('Trial expectation E_R', 'An average over the covariate distribution among trial participants, R = 1.', 'It standardizes the treatment contrast to the trial population.'),
    eifU: entry('Unrestricted hybrid EIF φᵉᶠᶠ_U(O)', 'The efficient influence function for the hybrid borrowing target when the two conditional control means may differ.', 'It includes the effects of estimating source availability and both control variances.'),
    eifTrial: entry('Trial-effect EIF under hybrid sampling', 'The efficient influence function for the trial effect in the unrestricted hybrid model, using only randomized trial outcomes.', 'Subtracting it from the borrowing-target EIF gives the unrestricted causal-gap EIF.'),
    eifGap: entry('Hybrid causal-gap EIF', 'The borrowing-target EIF minus the trial-effect EIF for the same hybrid record.', 'Its variance accounts for the covariance between the two paired contributions.'),
    eifTransport: entry('Hybrid EIF under control-mean transportability', 'The efficient influence function for the trial effect when the conditional control means agree across the two sources.', 'The source-specific control residuals combine and the share-derivative terms vanish.'),
    eifRCT: entry('Randomized trial EIF', 'The efficient influence function for the trial treatment effect under the randomized trial observation law.', 'The randomized treatment and control residuals are corrected using the known allocation probability.'),
    eifECT: entry('Externally controlled trial EIF', 'The efficient influence function for the trial-standardized external comparison under the two-source ECT observation law.', 'The trial records supply treated outcomes and the external records supply weighted control residuals.'),
    probability: entry('Probability Pr', 'The probability of the indicated event under the current observation law.', 'A vertical bar introduces the covariates or group on which the probability is conditioned.'),
    meanOperator: entry('Conditional expectation E', 'The average outcome among records satisfying the stated conditions.', 'Here it defines the mean outcome at given covariates within an observed group.'),
    varianceOperator: entry('Conditional variance Var', 'The variance of outcomes among records satisfying the stated conditions.', 'Here it gives the within-group outcome variability at the specified covariates.'),
    ectProbability: entry('ECT probability operator Pr_E', 'A probability under the two-source law with only treated trial records and external control records.', 'Conditional and marginal probabilities in the ECT use this sampling law.'),
    ectExpectation: entry('ECT expectation E_E', 'An average under the law of treated trial participants and external controls in the ECT.', 'Conditioning this average on R = 1 standardizes to the ECT trial population.'),
    ectSource: entry('ECT conditional trial probability p_E(X)', 'The probability that a record belongs to the trial given X under the two-source ECT law.', 'Its odds reweight external-control outcome residuals toward the ECT trial population.'),
    ectFraction: entry('ECT trial sampling fraction q_E', 'The marginal probability that a record in the ECT belongs to its treated trial group.', 'It normalizes the EIF for the two-source ECT sample.'),
    ectMean1: entry('ECT treated mean m₁,E(X)', 'The mean outcome at X among treated participants in the ECT trial.', 'It is the treated component of the ECT target.'),
    ectMean0: entry('ECT external control mean m₀,E(X)', 'The mean outcome at X among external controls in the ECT.', 'Standardizing this mean to ECT trial covariates supplies the external comparator.'),
    ectTarget: entry('ECT borrowing target θ̃_E', 'The mean treated-minus-external-control outcome contrast standardized to the ECT trial covariate distribution.', 'A causal interpretation as the trial effect requires an assumption linking the unobserved trial control mean to the external mean.'),
    trialDensity: entry('Original trial covariate density f_R(X)', 'The covariate density among trial participants before the hybrid design is reduced.', 'With covariate-dependent allocation, selecting only treated trial records reweights this density by e(X).'),
    ectTrialDensity: entry('Retained trial covariate density f_R,E(X)', 'The covariate density among the retained treated trial participants after internal controls are removed.', 'It equals e(X)f_R(X)/E_R{e(X)} under this reduction.')
  });

  // A compound token has one tab stop, including its subscripts and arguments.
  function mark(key, expression, vocabulary = notation) {
    const item = vocabulary[key];
    if (!item) throw new Error('Unknown EIF notation: ' + key);
    const flat = expression.replace(/<mrow data-notation="[^"]+" role="button" tabindex="0" aria-label="[^"]+">/g, '<mrow>');
    const label = item.label.replace(/&/g, '&amp;').replace(/"/g, '&quot;').replace(/</g, '&lt;');
    return '<mrow data-notation="' + key + '" role="button" tabindex="0" aria-label="' + label + '">' + flat + '</mrow>';
  }

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
  const rawX = id('X'), rawO = id('O'), rawTheta = id('θ'), rawQ = id('q');
  const X = mark('covariates', rawX), O = mark('record', rawO);
  const R = mark('source', id('R')), A = mark('treatment', id('A')), Y = mark('outcome', id('Y'));
  const zero = num(0), one = num(1), q = mark('trialFraction', rawQ), theta = mark('trialEffect', rawTheta);
  const eq = op('='), plus = op('+'), minus = op('−');
  const diff = (s, t) => s + minus + t;
  const complement = s => diff(one, s);
  const given = s => op('∣') + s;
  const term = (s, t) => id(s) + wrap(t);
  const atX = s => s + wrap(rawX);
  const idx = g => String(g) === '10' ? one + op(',') + zero : num(g);
  const I = g => mark('indicator' + g, sub(id('I'), idx(g)));
  const m = g => mark('mean' + g, atX(sub(id('m'), idx(g))));
  const v = g => mark('variance' + g, atX(sub(id('v'), idx(g))));
  const res = g => mark('residual' + g, sub(id('ε'), idx(g)) + wrap(rawO));
  const squaredResidual = g => braces(diff(frac(sup(res(g), num(2)), v(g)), one));
  const p = mark('sourceProbability', term('p', rawX));
  const e = mark('randomization', term('e', rawX));
  const r = mark('varianceRatio', term('r', rawX));
  const d = mark('weightDenominator', term('d', rawX));
  const l = mark('internalShare', term('λ', rawX));
  const delta = mark('discrepancy', term('δ', rawX));
  const k = mark('derivativeMultiplier', term('c', rawX));
  const t = mark('conditionalContrast', atX(sub(id('t'), id('λ'))));
  const ttilde = '<mover>' + rawTheta + '<mo>~</mo></mover>';
  const target = mark('borrowingTarget', sub(ttilde, id('λ')));
  const gap = mark('causalGap', sub(id('ξ'), id('λ') + op(',') + id('δ')));
  const pooled = mark('pooledMean', atX(sub(id('m'), roman('P'))));
  const expected = mark('trialExpectation', sub(roman('E'), roman('R')));
  const phiKeys = {U: 'eifU', 'trial,H': 'eifTrial', 'HCT,trans': 'eifTransport', ξ: 'eifGap', RCT: 'eifRCT', ECT: 'eifECT'};
  const phi = g => mark(phiKeys[g], sup(sub(id('φ'), roman(g)), roman('eff')) + wrap(rawO));
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
  const eg = mark('residualG', sub(id('ε'), id('g')) + wrap(rawO));
  const mg = mark('meanG', atX(sub(id('m'), id('g'))));
  const vg = mark('varianceG', atX(sub(id('v'), id('g'))));
  const groupCondition = given(X + op(',') + mark('indicatorG', sub(id('I'), id('g'))) + eq + one);
  const groupDefinition =
    display('Record and group definitions', [
      [O, eq + wrap(X + op(',') + R + op(',') + A + op(',') + Y)],
      [I(1), eq + R + A],
      [I('10'), eq + R + wrap(complement(A))],
      [I(0), eq + complement(R)],
      [p, eq + mark('probability', roman('Pr')) + wrap(R + eq + one + given(X))],
      [q, eq + mark('probability', roman('Pr')) + wrap(R + eq + one)],
      [e, eq + mark('probability', roman('Pr')) + wrap(A + eq + one + given(X + op(',') + R + eq + one))],
      [mg, eq + mark('meanOperator', roman('E')) + wrap(Y + groupCondition)],
      [vg, eq + mark('varianceOperator', roman('Var')) + wrap(Y + groupCondition)],
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
    'The trial is indicated by R = 1. Groups 1, (1, 0), and 0 denote trial treatment, trial control, and external control. ' +
    'Treatment allocation is known. The expectation E<sub>R</sub> averages over trial covariates.' ) +
    groupDefinition + weightDefinition + para(
    'The precision rule is defined on shared support with positive control variances. At trial-only profiles use the internal share of one and evaluate only the trial contributions. At external-only profiles the EIF is zero. The variance bound assumes square-integrable contributions.'));
  const full = para('The target is the treated mean minus the precision-weighted control mean, standardized to trial covariates. Write ' +
    math(k + eq + cdef) + ' and ' + math(eg + eq + diff(Y, mg)) + '.') +
    display('Unrestricted hybrid trial efficient influence function', [
      [hct, eq + frac(R, q) + braces(diff(t, target))],
      ['', plus + trialResidual],
      ['', minus + internalResidual],
      ['', minus + externalResidual],
      ['', minus + frac(k, q) + frac(diff(R, p), complement(p))],
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
        [trial, eq + frac(R, q) + braces(m(1) + minus + m('10') + minus + theta)],
        ['', plus + trialResidual + minus + frac(I('10'), q + braces(complement(e))) + res('10')]
      ]) +
      para('Pair the two gradients for each participant when estimating the gap variance. Their shared trial data induce covariance.'));

  const transport = para('Retain the three groups and impose equal conditional control means on shared support, ' +
    math(m('10') + eq + m(0) + eq + pooled) + '. The control variances may differ.') +
    display('Transportable hybrid trial efficient influence function', [
      [common, eq + frac(R, q) + braces(m(1) + minus + pooled + minus + theta)],
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

  const rct = para('Removing external records leaves a randomized trial. The internal share is fixed at one and the borrowing target equals the trial effect. Each record now belongs to the trial, so R = q = 1.') +
    display('Randomized controlled trial efficient influence function', [
      [phi('RCT'), eq + m(1) + minus + m('10') + minus + theta],
      ['', plus + frac(A,e) + res(1)],
      ['', minus + frac(complement(A),complement(e)) + res('10')]
    ]) +
    para('The external residual and all three share derivatives drop out. The two remaining residuals correct the treatment and control means with their known randomization probabilities.') +
    disclosure('Using only trial outcomes in a hybrid sample',
      para('If all hybrid records are retained but the internal share is fixed at one, the external records contribute zero. The trial contribution still has to be normalized for the hybrid mixture. Conditioning the sample on trial membership gives the pure RCT expression above.') +
      display('Trial-only gradient normalized for HCT mixture sampling', [
        [trial, eq + frac(R,q) + braces(m(1) + minus + m('10') + minus + theta)],
        ['', plus + trialResidual + minus + frac(I('10'), q + braces(complement(e))) + res('10')]
      ]) +
      para('Here R indicates trial membership, q is its mixture probability, and the two indicators select randomized treatment and control. The EIF requires both randomized arms on trial support.'));

  const pE = mark('ectSource', atX(subscript(id('p'), 'E')));
  const qE = mark('ectFraction', subscript(rawQ, 'E'));
  const mE = g => mark('ectMean' + g, atX(sub(id('m'), idx(g) + op(',') + roman('E'))));
  const thetaE = mark('ectTarget', subscript(ttilde,'E'));
  const exE = mark('ectExpectation', subscript(roman('E'),'E'));
  const prE = mark('ectProbability', subscript(roman('Pr'),'E'));
  const fR = mark('trialDensity', atX(sub(id('f'), roman('R'))));
  const fRE = mark('ectTrialDensity', atX(sub(id('f'), roman('R,E'))));
  const Eres = braces(diff(Y,mE(0)));
  const ect = para('All trial participants receive treatment, and all observed controls come from the external source. Under this two-group law, the estimand is the trial-standardized external comparison.') +
    display('Externally controlled trial target and efficient influence function', [
      [thetaE, eq + exE + braces(mE(1) + minus + mE(0) + given(R + eq + one))],
      [phi('ECT'), eq + frac(R,qE) + braces(Y + minus + mE(0) + minus + thetaE)],
      ['', minus + frac(wrap(complement(R)) + pE, qE + braces(complement(pE))) + Eres]
    ]) +
    para('The internal share is fixed at zero. Under the ECT law the treated residual and trial standardization combine because all trial participants are treated. The external residual is reweighted by the ECT source odds. With no internal control arm, the data identify the external comparison. A causal interpretation requires an assumption linking the unobserved trial control mean to the external mean.') +
    display('Combining the two treated terms', [
      [braces(mE(1) + minus + mE(0) + minus + thetaE) + plus + braces(diff(Y,mE(1))), eq + Y + minus + mE(0) + minus + thetaE]
    ]) +
    disclosure('Recomputing the observation law',
      para('Here ' + math(pE + eq + prE + wrap(R + eq + one + given(X))) + ' and ' +
      math(qE + eq + prE + wrap(R + eq + one)) +
      '. If the ECT is formed by deleting internal controls from an HCT, condition the original observation law on retaining a treated trial or external record. This gives') +
      display('Two-source probabilities after deleting internal controls', [
        [pE, eq + frac(p + e, complement(p) + plus + p + e)],
        [qE, eq + frac(q + expected + braces(e), complement(q) + plus + q + expected + braces(e))],
        [fRE, eq + frac(e + fR, expected + braces(e))]
      ]) + para('The last line gives the retained trial covariate density. Constant allocation preserves the original trial covariate distribution. With covariate-dependent allocation, retaining the original trial target requires restandardization. In a separately planned ECT, use its actual two-source law. External controls must cover the target trial support, and the displayed terms must be square-integrable.'));

  const panels = Object.freeze({full, transport, rct, ect});
  const helpers = Object.freeze({
    id, num, op, roman, row, sub, sup, frac, wrap, braces, math, display, term, atX, idx, I, m, v, res,
    squaredResidual, phi, X, O, R, A, Y, zero, one, q, theta, eq, plus, minus, diff, complement,
    given, p, e, r, d, l, delta, k, t, target, gap, pooled, expected, hct, trial, common,
    trialResidual, internalResidual, externalResidual, cdef, eg, mg, vg, pE, qE, mE, thetaE, exE,
    Eres, fR, fRE, mark
  });
  if (typeof window !== 'undefined') {
    window.EifFormulaPanels = panels;
    window.EifNotation = notation;
    window.EifMath = helpers;
  }
  if (typeof module !== 'undefined' && module.exports) module.exports = panels;
}());
