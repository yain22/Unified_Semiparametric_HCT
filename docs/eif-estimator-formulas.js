/* Sample estimators in Section 4. The formulas use the same MathML and glossary
   as the population influence functions on this page. */
(function () {
  'use strict';

  const {id, roman, num, op, sub, sup, frac, wrap, braces, display, math, mark,
    eq, minus, plus, one, diff, complement} = window.EifMath;
  const notation = {...window.EifNotation};
  const i = id('i'), comma = op(','), zero = num(0);
  const hat = a => '<mover accent="true">' + a + '<mo>^</mo></mover>';
  const bar = a => '<mover accent="true">' + a + '<mo>¯</mo></mover>';
  const sq = a => sup(a, num(2));
  const transpose = a => sup(a, op('⊤'));
  const at = (f, a) => f + wrap(a);
  const indexed = (name, index) => sub(id(name), index);
  const bracesSq = a => sq(braces(a));
  const P = text => '<p>' + text + '</p>';
  const details = (label, body) => '<details><summary>' + label + '</summary>' + body + '</details>';
  const sum = (upper, index = i) => '<munderover><mo>∑</mo><mrow>' + index + eq + one +
    '</mrow><mrow>' + upper + '</mrow></munderover>';
  const sumIn = condition => '<munder><mo>∑</mo><mrow>' + i + op(':') + condition + '</mrow></munder>';

  function symbol(key, expression, label, definition, role) {
    notation[key] = {label, definition, role};
    return mark(key, expression, notation);
  }
  const X = mark('covariates', indexed('X', i));
  const O = mark('record', indexed('O', i));
  const S = mark('source', indexed('S', i));
  const A = mark('treatment', indexed('A', i));
  const Y = mark('outcome', indexed('Y', i));
  const e = mark('randomization', at(id('e'), X));
  const g10 = one + comma + zero;
  const G = (key, index) => mark(key, indexed('I', index + comma + i));
  const I1 = G('indicator1', one), I10 = G('indicator10', g10), I0 = G('indicator0', zero);
  const N = symbol('fitN', id('N'), 'Hybrid sample size N',
    'The total number of treated trial participants, randomized controls, and external controls.',
    'All N records enter the hybrid sum.');
  const nR = symbol('fitNR', indexed('n', roman('RCT')), 'Randomized trial size n_RCT',
    'The total number in the two randomized arms.', 'This count normalizes the hybrid sum. It is also the sample size of the trial-only regression.');
  const nE = symbol('fitNE', indexed('n', roman('ECT')), 'Externally controlled trial size n_ECT',
    'The number of treated trial records plus the external controls in this analysis.', 'Both groups enter the ECT sum.');
  const n1 = symbol('fitN1', indexed('n', one), 'Treated trial sample size n₁',
    'The number of treated trial participants retained in this analysis.', 'In the ECT they define the population used for standardization.');
  const n10 = symbol('fitN10', indexed('n', g10), 'Internal control sample size n₁,₀',
    'The number of randomized trial controls.', 'Each control outcome enters the unadjusted trial contrast with weight −1/n₁,₀.');

  function fittedAt(key, name, index, argument, label, definition, role) {
    const f = index ? sub(hat(id(name)), index) : hat(id(name));
    return symbol(key, at(f, argument), label, definition, role);
  }
  const groupNames = ['treated trial participants', 'internal controls', 'external controls'];
  const groupIndices = [one, g10, zero];
  const groupLabels = ['₁', '₁,₀', '₀'];
  const m = groupIndices.map((g, j) => fittedAt('fittedMean' + j, 'm', g, X,
    'Fitted mean m̂' + groupLabels[j] + '(Xᵢ)', 'The predicted outcome for ' + groupNames[j] + ' at participant i’s baseline covariates.',
    'For cross-fitting the prediction comes from the training folds that exclude participant i.'));
  const residual = groupIndices.map((g, j) => fittedAt('fittedResidual' + j, 'ε', g, O,
    'Fitted residual ε̂' + groupLabels[j] + '(Oᵢ)', 'The observed outcome minus its prediction from the outcome model for ' + groupNames[j] + '.',
    'The corresponding group indicator selects the records that supply this residual.'));
  const v = (j) => fittedAt('fittedVariance' + j, 'v', groupIndices[j], X,
    'Fitted variance v̂' + groupLabels[j] + '(Xᵢ)', 'The predicted conditional outcome variance for ' + groupNames[j] + ' at these covariates.',
    'It determines the precision share and scales the raw squared residual in the share correction.');
  const v10 = v(1), v0 = v(2);
  const p = fittedAt('fittedSource', 'p', '', X, 'Fitted trial probability p̂(Xᵢ)',
    'The estimated probability of trial membership at Xᵢ under the three-group hybrid design.',
    'It determines the source odds and the relative availability of internal and external controls.');
  const omegaH = fittedAt('fittedOddsH', 'ω', roman('H'), X, 'Hybrid source odds ω̂_H(Xᵢ)',
    'The fitted hybrid probability p̂(Xᵢ) divided by 1 − p̂(Xᵢ) on shared support.',
    'It reweights an external record toward the distribution of trial covariates.');
  const lambda = fittedAt('fittedShare', 'λ', '', X, 'Fitted internal-control share λ̂(Xᵢ)',
    'The internal-control weight from the fitted source probability and ratio of the two fitted control variances.',
    'The complementary weight goes to the external mean. Use a share of one at trial-only profiles.');
  const delta = fittedAt('fittedDrift', 'δ', '', X, 'Fitted control-mean difference δ̂(Xᵢ)',
    'The predicted internal control mean minus the predicted external control mean at the same covariates.',
    'It multiplies the three corrections for fitting the precision share.');
  const t = fittedAt('fittedContrast', 't', '', X, 'Fitted conditional contrast t̂(Xᵢ)',
    'The fitted treated mean minus the precision-weighted average of the two fitted control means.',
    'The trial participants supply the covariates at which this contrast is averaged.');
  const mP = fittedAt('fittedPooled', 'm', roman('P'), X, 'Fitted pooled control mean m̂_P(Xᵢ)',
    'The predicted outcome from one control regression using both sources on shared support.',
    'Under conditional control-mean transportability it estimates the common control mean.');
  const rho = symbol('workingRatio', at(id('ρ'), X), 'Working variance ratio ρ(Xᵢ)',
    'A positive working ratio of internal to external control variances.',
    'The true ratio gives the efficient weights under equal control means. The clinical comparator uses one.');
  const dRho = fittedAt('fittedDenom', 'd', id('ρ'), X, 'Fitted weight denominator d̂ρ(Xᵢ)',
    'The sum p̂(Xᵢ){1 − e(Xᵢ)} + {1 − p̂(Xᵢ)}ρ(Xᵢ).',
    'It scales the two groups’ residual corrections for the pooled control mean.');
  const mE = fittedAt('fittedExternalECT', 'm', zero + comma + roman('E'), X,
    'Fitted external mean m̂₀,E(Xᵢ)', 'The external control outcome predicted at Xᵢ from the ECT analysis cohort.',
    'Predictions at the treated patients’ covariates standardize the comparator to the treated group.');
  const oddsE = fittedAt('fittedOddsECT', 'ω', roman('E'), X, 'ECT source odds ω̂_E(Xᵢ)',
    'The fitted trial-membership odds under the treated-plus-external sampling law.',
    'External residuals receive these weights. The clinical fit clips the source probability to [0.01, 0.99] before forming the odds.');
  const resE = fittedAt('fittedResidualECT', 'ε', zero + comma + roman('E'), O, 'ECT residual ε̂₀,E(Oᵢ)',
    'The observed outcome Yᵢ minus m̂₀,E(Xᵢ), evaluated for either source.',
    'The treated records supply the standardized contrast and the external records correct its comparator.');

  function pointEstimate(key, label, name, definition) {
    return symbol('estimate' + key, sub(hat(id('θ')), roman(label)), name, definition,
      'The result is measured in the same units as the outcome.');
  }
  const thetaU = pointEstimate('U','AIPW-U','Unrestricted hybrid estimate',
    'The fitted precision-share treatment contrast standardized to the trial covariate distribution.');
  const thetaR = pointEstimate('R','AIPW-R','Restricted hybrid estimate',
    'The hybrid estimator using a common control regression. Under conditional control-mean transportability it estimates the trial effect.');
  const thetaE = pointEstimate('E','E','Externally controlled trial estimate',
    'The treated-versus-external outcome contrast standardized to the treated trial participants.');
  const thetaTrial = pointEstimate('Trial','RCT','Randomized trial estimate',
    'The fitted treatment coefficient from the selected regression using only randomized participants.');
  const H = symbol('fittedH', at(sub(id('H'), sub(hat(id('η')), i)), O), 'Fitted hybrid numerator Hη̂ᵢ(Oᵢ)',
    'The sum of the four mean and standardization contributions and the three share-derivative contributions.',
    'Each observation supplies one numerator. The sum is divided by the number of trial participants.');
  const HM = symbol('fittedHM', sub(hat(id('H')), roman('M') + comma + i), 'Mean contribution Ĥ_M,i',
    'Trial standardization plus the three group-specific residual corrections.',
    'This component persists when the entire borrowing share is prespecified.');
  const HL = symbol('fittedHL', sub(hat(id('H')), id('λ') + comma + i), 'Share contribution Ĥ_λ,i',
    'The corrections for fitting the source probability and the two control variances.',
    'These terms enter because the precision share determines the borrowing target.');
  const phi = symbol('fittedInfluence', sub(hat(id('φ')), i), 'Fitted participant influence value φ̂ᵢ',
    'The fitted numerator minus the estimated effect times Sᵢ, divided by the observed trial fraction.',
    'Its sample variance divided by N estimates the estimator variance under iid sampling and the stated convergence conditions.');
  const phiMean = symbol('fittedInfluenceMean', bar(hat(id('φ'))), 'Mean fitted influence value',
    'The arithmetic mean of the N fitted participant influence values.', 'Subtract it to form the sample variance in equation (17).');
  const qhat = symbol('fittedFraction', hat(id('q')), 'Observed trial fraction q̂',
    'The number of randomized trial participants divided by the total hybrid sample size.', 'It normalizes the fitted influence values under mixture sampling.');
  const varianceLabel = (key, index, label, definition) => symbol(key, sq(sub(hat(id('s')), roman(index))),
    label, definition, 'Taking the square root gives the estimated standard error.');
  const sU = varianceLabel('fitVarianceIID','iid','Estimated iid sampling variance',
    'Sample variance of the fitted influence values divided by the total number of observations.');
  const sE = varianceLabel('fitVarianceECT','E','ECT plug-in variance',
    'The sum of the squared fitted EIF contributions on the scale of the ECT estimator.');
  const sRCT = varianceLabel('fitVarianceRCT','RCT','HC2 variance for the treatment coefficient',
    'A heteroskedasticity-consistent variance using squared regression weights and leverage-adjusted residual squares.');

  const ratioCorrection = (j, variance) => braces(diff(frac(sq(residual[j]), variance), one));
  const hybridSum = frac(one, nR) + sum(N);
  const treatedTerm = frac(I1, e) + residual[0];
  const internalTerm = frac(I10 + lambda, complement(e)) + residual[1];
  const externalTerm = I0 + omegaH + braces(complement(lambda)) + residual[2];
  const multiplier = delta + lambda + braces(complement(lambda));
  const expandedU = display('Equation 16 expanded over all hybrid records', [
    [thetaU, eq + hybridSum + op('[') + S + t + plus + treatedTerm],
    ['', minus + internalTerm + minus + externalTerm],
    ['', plus + multiplier + op('(') + minus + frac(diff(S, p), complement(p))],
    ['', plus + frac(I10, complement(e)) + ratioCorrection(1, v10)],
    ['', minus + I0 + omegaH + ratioCorrection(2, v0) + op(')') + op(']')]
  ]);
  const fitDefinitions = display('Quantities evaluated from the hybrid fits', [
    [t, eq + m[0] + minus + lambda + m[1] + minus + braces(complement(lambda)) + m[2]],
    [delta, eq + diff(m[1], m[2])],
    [omegaH, eq + frac(p, complement(p))],
    [lambda, eq + frac(p + braces(complement(e)),
      p + braces(complement(e)) + plus + braces(complement(p)) + frac(v10, v0))]
  ]);
  const unrestricted = P('Fit separate means for trial treatment, internal control, and external control. Each hybrid record contributes to the sum in equation (16). The denominator counts the trial participants whose covariates define the target.') +
    display('Unrestricted hybrid estimator, equation 16', [
      [thetaU, eq + hybridSum + H], [H, eq + HM + plus + HL]
    ]) +
    details('Expand the participant sum', expandedU + P('The first four terms use the fitted contrast and the three outcome residuals. The remaining three correct for fitting the source probability and each control variance. The variance corrections use raw squared residuals.') + fitDefinitions) +
    details('Uncertainty and support', display('Plug-in variance under independent mixture sampling, equation 17', [
      [phi, eq + frac(diff(H, S + thetaU), qhat)],
      [qhat, eq + frac(nR, N)],
      [sU, eq + frac(one, N + wrap(diff(N, one))) + sum(N) + bracesSq(diff(phi, phiMean))]
    ]) + P('With cross-fitting, the fitted quantities for participant i come from folds that exclude that participant. Equation (17) uses independent mixture sampling and requires score convergence and negligible first-order bias from nuisance fitting. The interval measures uncertainty around the hybrid statistical target.') +
    P('At trial-only covariate profiles, set the internal share to one, use only trial terms, and set the share correction to zero. External-only profiles contribute zero. The paper uses a within-group sandwich for its fixed-count simulations and a refitted, stratified bootstrap for the clinical intervals.'));

  const pooledCorrection = frac(p + braces(I10 + plus + I0 + rho), dRho) + braces(diff(Y, mP));
  const HR = symbol('fittedHR', indexed('H', roman('R') + comma + i), 'Restricted numerator H_R,i',
    'The terms inside the square brackets in equation (18), evaluated at participant i.',
    'It combines trial standardization with a treated residual and a pooled control residual.');
  const transportable = P('Fit one control regression on shared support and specify a positive working variance ratio. Write H<sub>R,i</sub> for the fitted participant contribution in equation (18).') +
    display('Restricted hybrid estimator', [[thetaR, eq + hybridSum + HR]]) +
    details('Expand the participant sum', display('Restricted hybrid estimator, equation 18', [
      [thetaR, eq + hybridSum + op('[') + S + braces(diff(m[0], mP)) + plus + treatedTerm],
      ['', minus + pooledCorrection + op(']')],
      [dRho, eq + p + braces(complement(e)) + plus + braces(complement(p)) + rho]
    ]) + P('The two control groups share one fitted mean. Their residual weights depend on source membership and the working variance ratio. In the clinical comparison the working ratio equals one.')) +
    details('Interpretation and support', P('Under equal conditional control means, a consistent estimate of the true control variance ratio gives efficient weighting when the nuisance conditions hold. Trial-only profiles use the trial control regression. External-only profiles contribute zero.') +
    P('If the conditional means differ, the conventional score with correct source and allocation probabilities targets the contrast induced by the working ratio. The refitted clinical bootstrap propagates sampling variation in this procedure.'));

  const weightE = braces(diff(S, wrap(complement(S)) + oddsE));
  const ect = P('Fit the external outcome regression and source odds using treated trial participants and external controls. Divide by the number of treated participants on the selected target support.') +
    display('Externally controlled trial estimator, equation 19', [
      [thetaE, eq + frac(one, n1) + sum(nE) + weightE + resE]
    ]) + details('Expand the participant sum', display('ECT treated and external sums', [
      [thetaE, eq + frac(one, n1) + sumIn(S + eq + one) + braces(diff(Y, mE))],
      ['', minus + frac(one, n1) + sumIn(S + eq + zero) + oddsE + braces(diff(Y, mE))]
    ]) + P('The first group compares its observed outcomes with predicted external control outcomes at its covariates. The weighted external residuals correct the standardized predictions.')) +
    details('Uncertainty and support', display('ECT plug-in variance, equation 19', [
      [sE, eq + frac(one, sq(n1)) + sum(nE) + op('[') + S + bracesSq(diff(resE, thetaE))],
      ['', plus + wrap(complement(S)) + sq(oddsE) + sq(resE) + op(']')]
    ]) + P('External controls must cover the trial target’s covariate profiles. Interpreting the statistical contrast as a trial effect also requires an assumption about the missing trial control mean.') +
    P('The implementation clips fitted ECT source probabilities to [0.01, 0.99]. If clipping changes the limiting odds, the usual protection against outcome-regression misspecification can fail. The plug-in variance requires convergence of the fitted score and a negligible remainder.'));

  const betaA = symbol('fittedBetaA', sub(hat(id('β')), id('A')), 'Treatment coefficient β̂_A',
    'The coefficient of randomized treatment in the fitted trial-only regression.', 'The selected design matrix determines whether baseline covariates adjust the estimate.');
  const b = symbol('trialBasis', at(id('b'), X), 'Baseline adjustment terms b(Xᵢ)',
    'Prespecified baseline covariates or transformations, excluding the treatment indicator and intercept.', 'The ANCOVA regression includes these terms with a common slope.');
  const L = symbol('trialWeights', id('L'), 'Least-squares matrix L',
    'For full-column-rank trial design matrix D, L = (DᵀD)⁻¹Dᵀ.', 'Multiplying L by the outcome vector gives the fitted regression coefficients.');
  const LA = symbol('trialOutcomeWeight', indexed('L', id('A') + i), 'Outcome weight L_Ai',
    'The entry for participant i in the treatment-coefficient row of L.', 'Multiplying each trial outcome by its entry and summing gives the treatment estimate.');
  const leverage = symbol('trialLeverage', indexed('h', i), 'Leverage hᵢ',
    'The ith diagonal entry of DL, the projection matrix of the fitted regression.', 'HC2 divides the squared residual by 1 − hᵢ.');
  const rTrial = symbol('trialResidualFitted', sub(hat(id('r')), i), 'Trial regression residual r̂ᵢ',
    'The observed trial outcome minus the outcome predicted by the fitted regression.', 'Its squared value enters the HC2 variance.');
  const betaB = symbol('trialBaselineSlope', sub(hat(id('β')), id('b')), 'Fitted common slope β̂_b',
    'The coefficients of the prespecified baseline terms in the ANCOVA regression.', 'They translate the arm difference in baseline covariates into the adjustment to the raw contrast.');
  function armMean(name, treated) {
    const group = treated ? 'treated arm' : 'randomized control arm';
    const what = name === 'Y' ? 'observed outcomes' : 'baseline adjustment vectors';
    return symbol('armMean' + name + (treated ? 'T' : 'C'), sub(bar(id(name)), treated ? one : g10),
      'Mean of ' + what + ' in the ' + group, 'The arithmetic mean of the ' + what + ' across records in this trial arm.',
      name === 'Y' ? 'The difference between the two arm means gives the unadjusted estimate.' :
        'The difference between the arms’ mean covariates determines the size of the ANCOVA correction.');
  }
  const differenceY = diff(armMean('Y', true), armMean('Y', false));
  const baselineDifference = diff(armMean('b', true), armMean('b', false));
  const D = symbol('trialDesignMatrix', id('D'), 'Trial regression design matrix D',
    'Each randomized participant supplies a row containing the intercept and treatment assignment. ANCOVA also includes prespecified baseline terms.',
    'The design matrix determines the least-squares weights and each record’s leverage.');
  const fittedAssignment = symbol('fittedTrialAssignment', sub(hat(id('A')), i), 'Fitted treatment indicator Âᵢ',
    'The fitted value from regressing treatment on an intercept and the baseline adjustment terms across the randomized participants.',
    'Subtracting this value from observed treatment gives the residual used to calculate the ANCOVA outcome weight.');
  const tildeA = '<mover accent="true">' + id('A') + '<mo>~</mo></mover>';
  const residualAssignment = symbol('residualTrialAssignment', sub(tildeA, i), 'Residualized treatment',
    'The observed treatment indicator minus its fitted value from the trial regression of treatment on the baseline terms and an intercept.',
    'Its value for participant i supplies the numerator of that participant’s weight. The denominator sums its square over all trial participants.');
  const trialSum = sub(op('∑'), i);
  const sumOverArm = assignment => sumIn(S + eq + one + comma + A + eq + assignment);
  function trialCase(label, terms) {
    // A phone can wrap at a mathematical operator without separating a sum from its outcome.
    return '<div class="trial-case"><h4>' + label + '</h4><div class="trial-case-equation">' +
      terms.map(term => math('<mstyle displaystyle="true">' + term + '</mstyle>')).join('') +
      '</div></div>';
  }
  const trial = P('Both analyses use the randomized participants. Each estimates the treatment coefficient in its own regression. Adding baseline covariates to the design matrix changes the weight assigned to each outcome.') +
    display('Shared trial regression form', [
      [thetaTrial, eq + betaA + eq + sum(nR) + LA + Y],
      [L, eq + sup(wrap(transpose(D) + D), minus + one) + transpose(D)]
    ]) +
    P('The ANOVA design row is ' + math(wrap(one + comma + A)) + '. Common-slope ANCOVA uses ' +
      math(wrap(one + comma + A + comma + transpose(b))) + '. Here ' + math(residualAssignment) +
      ' is the residual from regressing treatment on the baseline terms and an intercept in the trial.') +
    '<div class="trial-case-display" role="group" aria-label="Two forms of the randomized trial estimator">' +
    '<div class="trial-cases-lhs">' + math(thetaTrial + eq) + '</div><div class="trial-cases">' +
    '<svg class="trial-brace" aria-hidden="true" viewBox="0 0 12 100" preserveAspectRatio="none"><path d="M11 0 Q6 0 6 6 V43 Q6 49 1 50 Q6 51 6 57 V94 Q6 100 11 100" vector-effect="non-scaling-stroke"/></svg>' +
    trialCase('ANOVA (unadjusted)', [
      frac(one, n1) + sumOverArm(one) + Y, minus + frac(one, n10) + sumOverArm(zero) + Y, eq + differenceY
    ]) +
    trialCase('ANCOVA (baseline-adjusted)', [
      frac(trialSum + residualAssignment + Y, trialSum + sq(residualAssignment)),
      eq + wrap(differenceY), minus + transpose(wrap(baselineDifference)) + betaB
    ]) +
    '</div></div>' + details('HC2 variance for either regression', display('HC2 variance of the treatment coefficient', [
      [sRCT, eq + sum(nR) + sq(LA) + frac(sq(rTrial), complement(leverage))]
    ]) + P('Use the weights, fitted outcome residuals and leverages from the selected regression. The design matrix must have full column rank, with leverage below one for each contributing record. The paper assumes constant treatment allocation for these coefficient comparisons. With covariate-dependent allocation, recovering the trial effect generally requires trial standardization or appropriate weighting.'));

  const j = id('j');
  function trialWeights(adjustment) {
    if (adjustment === 'anova') {
      return display('ANOVA outcome weight', [[LA, eq + frac(A, n1) + minus + frac(complement(A), n10)]]);
    }
    return display('ANCOVA outcome weights by treatment residualization', [
      [residualAssignment, eq + diff(A, fittedAssignment)],
      [LA, eq + frac(residualAssignment, sum(nR, j) + sq(mark('residualTrialAssignment', sub(tildeA, j), notation)))]
    ]);
  }

  const views = {unrestricted, transportable, ect, rct: trial};
  window.EifNotation = Object.freeze(notation);
  window.EifEstimatorFormulas = Object.freeze({
    render: model => views[model],
    trialWeights
  });
}());
