'use strict';

// Run from the repository root: node analysis/validation/check_interactive_labs.js
const assert = require('node:assert/strict');
const { inspect } = require('../../docs/eif-lab.js');
const { designResult, evaluateDesigns } = require('../../docs/design-lab.js');

function close(actual, expected, tolerance, label) {
  const difference = Math.abs(actual - expected);
  assert.ok(difference <= tolerance,
    label + ': expected ' + expected + ', found ' + actual + ', difference ' + difference);
  return difference;
}

// Section 4.2 gives the bound and derivative at this law. The lab's discrete
// conditional outcomes share the Gaussian example's first four moments.
const reference = {
  delta: 0.2, v0: 1, pattern: 'constant', overlap: 'even',
  c: 1, n: 300, path: 'external-variance', target: 'tce'
};
const base = inspect(reference);
close(base.strata[0].lambda, 0.25, 1e-12, 'internal precision share');
close(base.population.gap, 0.15, 1e-12, 'causal gap');
close(base.totals.derivative, -0.0375, 1e-12, 'external log-variance derivative');
close(base.totals.bound, 4.528125, 1e-12, 'per-observation TCE bound');

const common = inspect({ ...reference, delta: 0, target: 'gap' });
close(common.population.gap, 0, 1e-12, 'gap at equal control means');
close(common.totals.trialBound, 9, 1e-12, 'unrestricted trial effect bound');
close(common.totals.tceBound, 4.5, 1e-12, 'TCE bound at the common-mean law');
close(common.totals.gapBound, 4.5, 1e-12, 'gap bound under unrestricted perturbations');

// The finite difference tilts cell probabilities and recomputes the target.
// The other derivative averages the analytic EIF times the observed-data score.
const paths = ['external-variance', 'internal-variance', 'external-mean', 'source', 'case-mix'];
const targets = ['tce', 'gap', 'trial'];
let cases = 0;
let largestDerivativeError = 0;
let largestVarianceError = 0;
for (const path of paths) {
  for (const delta of [-1, 0, 0.4, 1]) {
    for (const overlap of ['even', 'uneven']) {
      for (const pattern of ['constant', 'concentrated']) {
        for (const v0 of [0.25, 1, 4]) {
          for (const target of targets) {
            const settings = { ...reference, path, delta, overlap, pattern, v0, target };
            const { totals, finiteDifference } = inspect(settings);
            const label = [path, delta, overlap, pattern, v0, target].join('/');
            close(totals.centered, 0, 1e-12, 'EIF mean ' + label);
            close(totals.orthogonal, 0, 1e-12, 'randomization term times score ' + label);
            const derivativeError = close(totals.derivative, finiteDifference, 1e-7, 'derivative ' + label);
            // With q=1/2, e=2/3 and c=1, adding c*R*(A-e) adds q*e*(1-e)=1/9 to the variance.
            const varianceError = close(totals.augmented - totals.bound, 1 / 9, 1e-10, 'variance ' + label);
            largestDerivativeError = Math.max(largestDerivativeError, derivativeError);
            largestVarianceError = Math.max(largestVarianceError, varianceError);
            cases++;
          }
        }
      }
    }
  }
}
assert.equal(cases, 720);

// Scenario I, with 120 treated participants. Expected probabilities were
// independently calculated in base R with pnorm() and qnorm(0.975).
const planning = { rct: 120, hct: 60, external: 300 };
const designNames = ['rct', 'hct', 'ect'];
for (const row of [
  { delta: 0, gap: [0, 0, 0],
    nullRate: [0.025, 0.025, 0.025],
    benefitRate: [0.664645516115, 0.831779970078, 0.813487824514] },
  { delta: 0.5, gap: [0, 5 / 12, 0.5],
    nullRate: [0.025, 0.088286329858, 0.106149744417],
    benefitRate: [0.664645516115, 0.941769596057, 0.945589967694] }
]) {
  const designs = evaluateDesigns({ ...planning, delta: row.delta });
  designNames.forEach((key, i) => {
    const d = designs[key];
    const label = key + ', delta=' + row.delta;
    assert.equal(d.enrolled, [240, 180, 120][i], 'newly recruited ' + label);
    close(d.gap, row.gap[i], 1e-12, 'gap ' + label);
    close(d.falsePositive, row.nullRate[i], 8e-8, 'null rejection ' + label);
    close(d.efficacy, row.benefitRate[i], 8e-8, 'efficacy declaration ' + label);
  });
}

const n266 = designResult(266, 0, 0).efficacy;
const n267 = designResult(267, 0, 0).efficacy;
close(n266, 0.799588179177, 8e-8, 'RCT power with 266 controls');
close(n267, 0.800045210575, 8e-8, 'RCT power with 267 controls');
assert.ok(n266 < 0.8 && n267 >= 0.8, 'minimum control count for 80% power');
console.log(cases + ' EIF derivative checks passed. Largest absolute error: ' + largestDerivativeError.toExponential(2));
console.log('Variance identities passed. Largest absolute error: ' + largestVarianceError.toExponential(2));
console.log('Scenario I design probabilities and the 267-control threshold passed.');
