'use strict';

// The one-sided 0.025 test and prognostic R-squared are fixed by Scenario I.
const z975 = 1.959963984540054;
const residualVariance = 0.75 * 7.5 ** 2;
const designNames = {hct: 'HCT', ect: 'ECT', rct: 'RCT (ANCOVA)'};

function normalCdf(z) {
  // Abramowitz & Stegun 26.2.17. Absolute error below 8e-8.
  const a = Math.abs(z);
  const t = 1 / (1 + 0.2316419 * a);
  const poly = t * (0.319381530 + t * (-0.356563782 + t * (1.781477937 +
    t * (-1.821255978 + t * 1.330274429))));
  const tail = Math.exp(-0.5 * a * a) / Math.sqrt(2 * Math.PI) * poly;
  return z >= 0 ? 1 - tail : tail;
}

function designResult(internal, external, delta) {
  const controls = internal + external;
  const externalShare = external / controls;
  const gap = externalShare * delta;
  const se = Math.sqrt(residualVariance * (1 / 120 + 1 / controls));
  return {
    enrolled: 120 + internal, externalShare, gap, se,
    falsePositive: normalCdf(gap / se - z975),
    efficacy: normalCdf((2 + gap) / se - z975)
  };
}

function controlCounts(settings) {
  return {
    hct: [settings.hct, settings.external],
    ect: [0, settings.external],
    rct: [settings.rct, 0]
  };
}

function evaluateDesigns(settings) {
  return Object.fromEntries(Object.entries(controlCounts(settings)).map(([key, [internal, external]]) =>
    [key, designResult(internal, external, settings.delta)]));
}

function rCode(s) {
  const groups = Object.entries(controlCounts(s));
  return [
    '# Scenario I: each design enrolls 120 treated participants.',
    'design <- c(' + groups.map(([key]) => '"' + designNames[key] + '"').join(', ') + ')',
    'internal <- c(' + groups.map(([, [internal]]) => internal).join(', ') + ')',
    'external <- c(' + groups.map(([, [, external]]) => external).join(', ') + ')',
    'delta <- ' + s.delta.toFixed(2) + '  # internal minus external control mean',
    'v <- 0.75 * 7.5^2',
    'external_share <- external / (internal + external)',
    'gap <- external_share * delta',
    'se <- sqrt(v * (1/120 + 1/(internal + external)))',
    'data.frame(design, enrolled = 120 + internal, external_share, gap, se,',
    '           null_rejection = pnorm(gap/se - qnorm(0.975)),',
    '           declare_benefit = pnorm((2 + gap)/se - qnorm(0.975)))'
  ].join('\n');
}

function plotCurves(canvas, settings, focus) {
  const ctx = canvas.getContext('2d');
  const width = Math.max(1, canvas.clientWidth);
  const height = canvas.clientHeight;
  const pixels = window.devicePixelRatio || 1;
  canvas.width = Math.round(width * pixels);
  canvas.height = Math.round(height * pixels);
  ctx.scale(pixels, pixels);
  const left = 48, right = width - 16, top = 16, bottom = height - 50;
  const trend = evaluateDesigns({...settings, delta: 1});
  const maxRate = Math.ceil(10 * Math.max(0.10, trend.hct.falsePositive, trend.ect.falsePositive)) / 10;
  const x = delta => left + (delta + 1) * (right - left) / 2;
  const y = probability => bottom - probability * (bottom - top) / maxRate;
  ctx.font = '12px "Segoe UI", sans-serif';
  ctx.fillStyle = '#425660';
  ctx.textAlign = 'right';
  const tickStep = maxRate <= 0.30 ? 0.05 : 0.1;
  for (let r = 0; r <= maxRate + 1e-8; r += tickStep) {
    ctx.beginPath();
    ctx.moveTo(left, y(r)); ctx.lineTo(right, y(r));
    ctx.strokeStyle = '#d7e0e3'; ctx.stroke();
    ctx.fillText(Math.round(100 * r) + '%', left - 6, y(r) + 4);
  }
  ctx.textAlign = 'center';
  [-1, -0.5, 0, 0.5, 1].forEach(delta => {
    ctx.fillText(String(delta), x(delta), bottom + 19);
  });
  ctx.fillText('Internal minus external control mean (points)', (left + right) / 2, height - 6);
  ctx.save();
  ctx.beginPath(); ctx.moveTo(x(settings.delta), top); ctx.lineTo(x(settings.delta), bottom);
  ctx.strokeStyle = '#81939c'; ctx.setLineDash([2, 3]); ctx.stroke();
  ctx.restore();
  const order = Object.keys(trend).filter(key => key !== focus).concat(focus);
  const colors = Object.fromEntries(order.map(key =>
    [key, getComputedStyle(document.querySelector('.key .' + key)).color]));
  const selected = evaluateDesigns(settings);
  order.forEach(key => {
    const curve = trend[key];
    ctx.save();
    ctx.lineWidth = key === focus ? 2.7 : 1.7;
    ctx.strokeStyle = colors[key]; ctx.beginPath();
    for (let i = 0; i <= 160; i++) {
      const delta = -1 + i / 80;
      const rate = normalCdf(curve.externalShare * delta / curve.se - z975);
      if (i === 0) ctx.moveTo(x(delta), y(rate)); else ctx.lineTo(x(delta), y(rate));
    }
    ctx.stroke(); ctx.restore();
  });
  const positions = {hct: [-.75, .40], ect: [-.55, .62], rct: [-.32, .83]};
  function mark(key, delta, rate, size) {
    const cx = x(delta), cy = y(rate);
    ctx.beginPath();
    if (key === 'rct') {
      ctx.moveTo(cx, cy - size);
      ctx.lineTo(cx + size, cy + size);
      ctx.lineTo(cx - size, cy + size);
      ctx.closePath();
    } else {
      ctx.arc(cx, cy, size, 0, 2 * Math.PI);
    }
    ctx.fillStyle = key === 'ect' ? '#fff' : colors[key];
    ctx.fill();
    if (key === 'ect') {
      ctx.strokeStyle = colors[key]; ctx.lineWidth = 1.5; ctx.stroke();
    }
  }
  order.forEach(key => {
    const curve = trend[key];
    positions[key].forEach(delta =>
      mark(key, delta, normalCdf(curve.externalShare * delta / curve.se - z975), 3));
  });
  mark(focus, settings.delta, selected[focus].falsePositive, 4.5);
}

function startLab() {
  const el = id => document.getElementById(id);
  const inputs = {rct: el('rct-ic'), hct: el('hct-ic'), external: el('ec'), delta: el('delta')};
  const settings = () => Object.fromEntries(Object.entries(inputs).map(([k, v]) => [k, Number(v.value)]));
  const percent = x => (x * 100).toFixed(1) + '%';
  const text = (id, value) => { el(id).textContent = value; };
  function refresh() {
    const s = settings();
    Object.values(inputs).forEach(input => {
      text(input.id + '-out', input === inputs.delta ? s.delta.toFixed(2) + ' points' : input.value);
    });
    const results = evaluateDesigns(s);
    const tbody = el('results');
    tbody.replaceChildren();
    for (const [key, row] of Object.entries(results)) {
      const tr = document.createElement('tr');
      const th = document.createElement('th');
      th.scope = 'row';
      th.textContent = designNames[key];
      tr.appendChild(th);
      [row.enrolled, percent(row.externalShare), row.gap.toFixed(3), row.se.toFixed(3),
        percent(row.falsePositive), percent(row.efficacy)].forEach(value => {
        const td = document.createElement('td');
        td.textContent = value; tr.appendChild(td);
      });
      tbody.appendChild(tr);
    }
    const h = results.hct;
    let insight = 'With ' + s.hct + ' concurrent and ' + s.external +
      ' external controls, the HCT assigns ' + percent(h.externalShare) + ' of the control mean to external patients. ';
    if (s.delta === 0) {
      insight += 'At equal control means, the HCT has ' + percent(h.efficacy) + ' power and a 2.5% false-positive rate.';
    } else {
      insight += 'At this discrepancy it declares benefit in ' + percent(h.falsePositive) +
        ' of trials with no treatment effect and ' + percent(h.efficacy) + ' of trials with a two-point effect.';
    }
    if (s.rct === 266 || s.rct === 267) {
      const status = results.rct.efficacy >= 0.80 ? 'above' : 'below';
      insight += ' With ' + s.rct + ' randomized controls, ANCOVA power is ' +
        (100 * results.rct.efficacy).toFixed(3) + '%, just ' + status + ' the 80% target.';
    }
    text('interpretation', insight);
    text('rcode', rCode(s));
    const focus = el('focus').value;
    const atMinus = evaluateDesigns({...s, delta: -1})[focus].falsePositive;
    const atPlus = evaluateDesigns({...s, delta: 1})[focus].falsePositive;
    text('chart-note', 'The vertical line marks the selected discrepancy. For ' + designNames[focus] +
      ', null rejection is ' + percent(atMinus) + ' at −1 point, 2.5% at zero, and ' + percent(atPlus) + ' at +1 point.');
    plotCurves(el('plot'), s, focus);
  }
  for (const node of Object.values(inputs)) node.addEventListener('input', refresh);
  el('focus').addEventListener('change', refresh);
  el('reset').addEventListener('click', () => {
    for (const node of Object.values(inputs)) node.value = node.defaultValue;
    el('focus').value = 'hct'; refresh();
  });
  el('case-drift').addEventListener('click', () => { inputs.delta.value = 0.5; refresh(); });
  el('case-rct').addEventListener('click', () => { inputs.rct.value = 267; refresh(); });
  el('download-r').addEventListener('click', () => {
    const url = URL.createObjectURL(new Blob([rCode(settings()) + '\n'], {type: 'text/plain'}));
    const a = document.createElement('a');
    a.href = url; a.download = 'trial-design.R'; a.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
  window.addEventListener('resize', refresh);
  refresh();
}

if (typeof document !== 'undefined') startLab();
if (typeof module !== 'undefined' && module.exports) module.exports = {designResult, evaluateDesigns, rCode};
