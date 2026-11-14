/* The participant table shows the same numerator as the selected estimator. */
(() => {
  'use strict';

  const by = id => document.getElementById(id);
  const trialAdjustment = () => document.querySelector('input[name="trial-adjustment"]:checked').value;
  const titles = {
    unrestricted: 'AIPW-U with separate control means',
    transportable: 'AIPW-R with a pooled control mean',
    ect: 'AIPW with external controls',
    rct: 'Randomized trial estimators'
  };
  const format = value => (Math.abs(value) < .00005 ? 0 : value).toFixed(4).replace('-', '−');
  let participant = 'P01';
  let result;

  function showParticipant() {
    const record = result.rows.find(row => row.id === participant);
    by('estimator-parts-caption').textContent = 'Terms for ' + record.id + ' · ' + record.group;
    const parts = result.termLabels.map((label, i) => [label, record.parts[i]]);
    if (record.weight !== null) parts.unshift(['OLS weight', record.weight]);
    const nodes = [];
    for (const [label, value] of parts) {
      const name = document.createElement('dt');
      const amount = document.createElement('dd');
      name.textContent = label;
      amount.textContent = format(value);
      nodes.push(name, amount);
    }
    by('estimator-parts').replaceChildren(...nodes);
    for (const button of by('estimator-rows').querySelectorAll('button')) {
      button.setAttribute('aria-pressed', String(button.dataset.record === participant));
    }
  }

  function updateNumbers(model, adjustment) {
    result = EifEstimatorExample.evaluate(model, {
      adjustment, externalShift: Number(by('external-outcomes').value)
    });
    if (!result.rows.some(row => row.id === participant)) participant = result.rows[0].id;
    by('estimator-fit-note').textContent = result.note;
    const isTrial = model === 'rct';
    by('trial-weight-formula').innerHTML = isTrial ? EifEstimatorFormulas.trialWeights(adjustment) : '';
    by('trial-weight-heading').hidden = !isTrial;
    by('estimator-result-label').textContent = isTrial ? adjustment.toUpperCase() + ' estimate' : 'Estimated treatment contrast';
    by('estimator-rows-caption').textContent = isTrial ? adjustment.toUpperCase() + ' participant weights and contributions' : 'Participant contributions to the sum';
    by('estimator-result').textContent = format(result.estimate);
    by('estimator-arithmetic').textContent = model === 'rct'
      ? 'Sum of the weighted trial outcomes'
      : format(result.numerator) + ' ÷ ' + result.denominator + ' ' + result.denominatorLabel.toLowerCase();
    const rows = result.rows.map(record => {
      const tr = document.createElement('tr');
      const label = document.createElement('th');
      label.scope = 'row';
      const button = document.createElement('button');
      button.type = 'button';
      button.textContent = record.id;
      button.dataset.record = record.id;
      button.setAttribute('aria-label', 'Inspect ' + record.id + ', ' + record.group);
      button.setAttribute('aria-controls', 'estimator-parts');
      label.append(button);
      tr.append(label);
      const cells = [record.group, record.x.toFixed(2), record.y.toFixed(1)];
      if (isTrial) cells.push(format(record.weight));
      cells.push(format(record.value));
      for (const value of cells) {
        tr.insertCell().textContent = value;
      }
      return tr;
    });
    by('estimator-rows').replaceChildren(...rows);
    showParticipant();
  }

  function render(model) {
    by('estimator-heading').textContent = titles[model];
    by('adjustment-control').hidden = model !== 'rct';
    by('external-outcomes-control').hidden = model === 'rct';
    by('estimator-formula').innerHTML = EifEstimatorFormulas.render(model);
    updateNumbers(model, trialAdjustment());
  }

  by('estimator-rows').addEventListener('click', event => {
    const button = event.target.closest('button[data-record]');
    if (!button) return;
    participant = button.dataset.record;
    by('estimator-parts').parentElement.open = true;
    showParticipant();
  });
  by('external-outcomes').addEventListener('change', () => {
    const model = document.querySelector('input[name="model"]:checked').value;
    updateNumbers(model, trialAdjustment());
  });
  window.EifEstimatorView = { render, updateTrial: () => updateNumbers('rct', trialAdjustment()) };
})();
