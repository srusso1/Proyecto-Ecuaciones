/* =============================================
   SIMULACIÓN EPIDEMIOLÓGICA SIR — simulation.js
   ============================================= */

// ==================== 1. CONSTANTES Y ESTADO GLOBAL ====================

const COLORS = {
  susceptible: '#58a6ff',
  infectado: '#f85149',
  recuperado: '#3fb950',
  euler: '#f0883e',
  heun: '#d2a8ff',
  rk4: '#58a6ff'
};

const appState = {
  N: 10000,
  I0: 5,
  beta: 0.30,
  gamma: 0.10,
  vac: 0,
  days: 160,
  dt: 1.0
};

const scenarioPresets = {
  parciales: { N: 10000, I0: 5, beta: 0.45, gamma: 0.10, vac: 0, days: 160, dt: 1.0 },
  vacunacion: { N: 10000, I0: 5, beta: 0.30, gamma: 0.10, vac: 40, days: 160, dt: 1.0 },
  virtual: { N: 10000, I0: 5, beta: 0.08, gamma: 0.10, vac: 0, days: 160, dt: 1.0 },
  brote: { N: 10000, I0: 20, beta: 0.80, gamma: 0.05, vac: 0, days: 160, dt: 1.0 }
};

// ==================== 2. SOLVER RK4 ====================

function solveRK4(beta, gamma, N, S0, I0, R0, days, dt) {
  const n = Math.ceil(days / dt) + 1;
  const t = new Array(n);
  const S = new Array(n);
  const I = new Array(n);
  const R = new Array(n);

  t[0] = 0;
  S[0] = S0;
  I[0] = I0;
  R[0] = R0;

  for (let i = 0; i < n - 1; i++) {
    const s = S[i], inf = I[i], r = R[i];

    const k1s = -beta * s * inf / N;
    const k1i = beta * s * inf / N - gamma * inf;
    const k1r = gamma * inf;

    const s2 = s + 0.5 * dt * k1s;
    const i2 = inf + 0.5 * dt * k1i;
    const r2 = r + 0.5 * dt * k1r;
    const k2s = -beta * s2 * i2 / N;
    const k2i = beta * s2 * i2 / N - gamma * i2;
    const k2r = gamma * i2;

    const s3 = s + 0.5 * dt * k2s;
    const i3 = inf + 0.5 * dt * k2i;
    const r3 = r + 0.5 * dt * k2r;
    const k3s = -beta * s3 * i3 / N;
    const k3i = beta * s3 * i3 / N - gamma * i3;
    const k3r = gamma * i3;

    const s4 = s + dt * k3s;
    const i4 = inf + dt * k3i;
    const r4 = r + dt * k3r;
    const k4s = -beta * s4 * i4 / N;
    const k4i = beta * s4 * i4 / N - gamma * i4;
    const k4r = gamma * i4;

    S[i + 1] = s + (dt / 6) * (k1s + 2 * k2s + 2 * k3s + k4s);
    I[i + 1] = inf + (dt / 6) * (k1i + 2 * k2i + 2 * k3i + k4i);
    R[i + 1] = r + (dt / 6) * (k1r + 2 * k2r + 2 * k3r + k4r);

    if (S[i + 1] < 0) S[i + 1] = 0;
    if (I[i + 1] < 0) I[i + 1] = 0;
    if (R[i + 1] < 0) R[i + 1] = 0;

    t[i + 1] = (i + 1) * dt;
  }

  return { t, S, I, R };
}

// ==================== 3. SOLVER EULER ====================

function solveEuler(beta, gamma, N, S0, I0, R0, days, dt) {
  const n = Math.ceil(days / dt) + 1;
  const t = new Array(n);
  const S = new Array(n);
  const I = new Array(n);
  const R = new Array(n);

  t[0] = 0;
  S[0] = S0;
  I[0] = I0;
  R[0] = R0;

  for (let i = 0; i < n - 1; i++) {
    const s = S[i], inf = I[i], r = R[i];
    S[i + 1] = s + dt * (-beta * s * inf / N);
    I[i + 1] = inf + dt * (beta * s * inf / N - gamma * inf);
    R[i + 1] = r + dt * (gamma * inf);

    if (S[i + 1] < 0) S[i + 1] = 0;
    if (I[i + 1] < 0) I[i + 1] = 0;
    if (R[i + 1] < 0) R[i + 1] = 0;

    t[i + 1] = (i + 1) * dt;
  }

  return { t, S, I, R };
}

// ==================== 4. SOLVER HEUN (EULER MEJORADO) ====================

function solveHeun(beta, gamma, N, S0, I0, R0, days, dt) {
  const n = Math.ceil(days / dt) + 1;
  const t = new Array(n);
  const S = new Array(n);
  const I = new Array(n);
  const R = new Array(n);

  t[0] = 0;
  S[0] = S0;
  I[0] = I0;
  R[0] = R0;

  for (let i = 0; i < n - 1; i++) {
    const s = S[i], inf = I[i], r = R[i];

    const f1s = -beta * s * inf / N;
    const f1i = beta * s * inf / N - gamma * inf;
    const f1r = gamma * inf;

    const sPred = s + dt * f1s;
    const iPred = inf + dt * f1i;
    const rPred = r + dt * f1r;

    const f2s = -beta * sPred * iPred / N;
    const f2i = beta * sPred * iPred / N - gamma * iPred;
    const f2r = gamma * iPred;

    S[i + 1] = s + (dt / 2) * (f1s + f2s);
    I[i + 1] = inf + (dt / 2) * (f1i + f2i);
    R[i + 1] = r + (dt / 2) * (f1r + f2r);

    if (S[i + 1] < 0) S[i + 1] = 0;
    if (I[i + 1] < 0) I[i + 1] = 0;
    if (R[i + 1] < 0) R[i + 1] = 0;

    t[i + 1] = (i + 1) * dt;
  }

  return { t, S, I, R };
}

// ==================== 5. FUNCIÓN computeStats ====================

function computeStats(solution) {
  const { t, S, I, R } = solution;
  let peakIdx = 0;
  let peakVal = I[0];
  for (let i = 1; i < I.length; i++) {
    if (I[i] > peakVal) {
      peakVal = I[i];
      peakIdx = i;
    }
  }
  const peakDay = t[peakIdx];
  const totalRecovered = R[R.length - 1];
  const totalAffected = R[R.length - 1];
  const population = S[0] + I[0] + R[0];
  const pctAffected = ((totalAffected / population) * 100);

  return { peakDay, peakVal, totalRecovered, pctAffected };
}

// ==================== 6. ESTADO DEL MODO PASO A PASO (FEATURE 1) ====================

const stepState = {
  active: false,
  currentIdx: 0,
  totalSteps: 0,
  cursorShapeIdx: null
};

function computeRK4Slopes(S, I, R, N, beta, gamma, dt) {
  const k1s = -beta * S * I / N;
  const k1i = beta * S * I / N - gamma * I;
  const k1r = gamma * I;

  const s2 = S + 0.5 * dt * k1s;
  const i2 = I + 0.5 * dt * k1i;
  const r2 = R + 0.5 * dt * k1r;
  const k2s = -beta * s2 * i2 / N;
  const k2i = beta * s2 * i2 / N - gamma * i2;
  const k2r = gamma * i2;

  const s3 = S + 0.5 * dt * k2s;
  const i3 = I + 0.5 * dt * k2i;
  const r3 = R + 0.5 * dt * k2r;
  const k3s = -beta * s3 * i3 / N;
  const k3i = beta * s3 * i3 / N - gamma * i3;
  const k3r = gamma * i3;

  const s4 = S + dt * k3s;
  const i4 = I + dt * k3i;
  const r4 = R + dt * k3r;
  const k4s = -beta * s4 * i4 / N;
  const k4i = beta * s4 * i4 / N - gamma * i4;
  const k4r = gamma * i4;

  return {
    S: [k1s, k2s, k3s, k4s],
    I: [k1i, k2i, k3i, k4i],
    R: [k1r, k2r, k3r, k4r]
  };
}

function updateStepDisplay() {
  const sol = window._lastSolution;
  if (!sol) return;
  const idx = stepState.currentIdx;
  const S = sol.S[idx];
  const I = sol.I[idx];
  const R = sol.R[idx];
  const beta = appState.beta;
  const gamma = appState.gamma;
  const dt = appState.dt;

  document.getElementById('step-val-S').textContent = Math.round(S);
  document.getElementById('step-val-I').textContent = Math.round(I);
  document.getElementById('step-val-R').textContent = Math.round(R);

  const dS = -beta * S * I / appState.N;
  const dI = beta * S * I / appState.N - gamma * I;
  const dR = gamma * I;

  document.getElementById('step-dS').textContent = dS.toFixed(2);
  document.getElementById('step-dI').textContent = dI.toFixed(2);
  document.getElementById('step-dR').textContent = dR.toFixed(2);

  document.getElementById('step-dS').style.color = dS < 0 ? 'var(--color-infectado)' : 'var(--color-recuperado)';
  document.getElementById('step-dI').style.color = dI < 0 ? 'var(--color-infectado)' : (dI > 0 ? 'var(--color-infectado)' : 'var(--text-secondary)');
  document.getElementById('step-dR').style.color = dR > 0 ? 'var(--color-recuperado)' : 'var(--text-secondary)';

  const slopes = computeRK4Slopes(S, I, R, appState.N, beta, gamma, dt);
  const labels = ['S', 'I', 'R'];
  labels.forEach(function(lbl) {
    for (let j = 0; j < 4; j++) {
      document.getElementById('k-' + lbl + (j + 1)).textContent = slopes[lbl][j].toFixed(4);
    }
  });

  const day = (idx * dt).toFixed(1);
  document.getElementById('step-indicator').textContent = 'Paso ' + idx + ' / ' + stepState.totalSteps + ' (día ' + day + ')';
}

function renderStepCursor() {
  const sol = window._lastSolution;
  if (!sol) return;
  const day = stepState.currentIdx * appState.dt;
  const chartDiv = document.getElementById('mainChart');
  if (!chartDiv.data) return;

  const existingLayout = chartDiv.layout || {};
  const existingShapes = (existingLayout.shapes || []).filter(function(s) {
    return s;
  });

  let cursorExists = false;
  const cleanShapes = existingShapes.map(function(s) {
    if (s.x0 !== undefined && s.x0 === stepState._lastDay &&
        s.line && s.line.color === '#58a6ff' && s.line.dash === 'dash') {
      cursorExists = true;
      return { ...s, x0: day, x1: day };
    }
    return s;
  });

  if (!cursorExists) {
    cleanShapes.push({
      type: 'line',
      x0: day,
      y0: 0,
      x1: day,
      y1: 1,
      yref: 'paper',
      line: { color: '#58a6ff', width: 2, dash: 'dash' }
    });
  }

  stepState._lastDay = day;

  Plotly.relayout('mainChart', { shapes: cleanShapes });
}

function removeStepCursor() {
  const chartDiv = document.getElementById('mainChart');
  if (!chartDiv || !chartDiv.layout) return;
  const existingShapes = (chartDiv.layout.shapes || []).filter(function(s) {
    return s && !(s.line && s.line.color === '#58a6ff' && s.line.dash === 'dash');
  });
  Plotly.relayout('mainChart', { shapes: existingShapes });
}

function toggleStepMode() {
  stepState.active = !stepState.active;
  const sol = window._lastSolution;
  const stepNav = document.getElementById('step-nav');
  const stepDetails = document.getElementById('step-details');

  if (stepState.active) {
    if (sol) {
      stepState.currentIdx = 0;
      stepState.totalSteps = sol.t.length - 1;
    }
    stepNav.style.display = 'flex';
    stepDetails.style.display = 'grid';
    updateStepDisplay();
    renderStepCursor();
  } else {
    stepNav.style.display = 'none';
    stepDetails.style.display = 'none';
    removeStepCursor();
  }
}

function stepNext() {
  if (!stepState.active) return;
  const sol = window._lastSolution;
  if (!sol || stepState.currentIdx >= sol.t.length - 1) return;
  stepState.currentIdx++;
  updateStepDisplay();
  renderStepCursor();
}

function stepPrev() {
  if (!stepState.active) return;
  if (stepState.currentIdx <= 0) return;
  stepState.currentIdx--;
  updateStepDisplay();
  renderStepCursor();
}

// ==================== 7. FUNCIÓN updateR0Display ====================

function updateR0Display(beta, gamma) {
  const r0 = beta / gamma;
  const el = document.getElementById('stat-r0');
  el.textContent = r0.toFixed(2);
  el.className = 'stat-value';
  if (r0 < 1) {
    el.classList.add('r0-green');
  } else if (r0 <= 2) {
    el.classList.add('r0-yellow');
  } else {
    el.classList.add('r0-red');
  }

  const r0Row = document.getElementById('r0-row');
  if (r0Row) {
    r0Row.className = '';
    if (r0 < 1) {
      r0Row.style.background = 'rgba(63, 185, 80, 0.08)';
    } else if (r0 <= 2) {
      r0Row.style.background = 'rgba(210, 153, 34, 0.08)';
    } else {
      r0Row.style.background = 'rgba(248, 81, 73, 0.08)';
    }
  }

  return r0;
}

// ==================== 8. FUNCIÓN updateDerivativesLive (FEATURE 3) ====================

function updateDerivativesLive(S, I, R, N, beta, gamma) {
  const dS = -beta * S * I / N;
  const dI = beta * S * I / N - gamma * I;
  const dR = gamma * I;

  const substS = '−' + beta.toFixed(2) + '·' + Math.round(S) + '·' + Math.round(I) + '/' + N;
  const substI = beta.toFixed(2) + '·' + Math.round(S) + '·' + Math.round(I) + '/' + N + ' − ' + gamma.toFixed(2) + '·' + Math.round(I);
  const substR = gamma.toFixed(2) + '·' + Math.round(I);

  const setDeriv = function(id, val, subst) {
    document.getElementById('deriv-' + id + '-val').textContent = val.toFixed(2);
    document.getElementById('deriv-' + id + '-subst').textContent = subst;
    const el = document.getElementById('deriv-' + id + '-val');
    el.className = 'deriv-value';
    if (val < 0) {
      el.classList.add('negative');
    } else if (val > 0) {
      el.classList.add('positive');
    }
  };

  setDeriv('dS', dS, substS);
  setDeriv('dI', dI, substI);
  setDeriv('dR', dR, substR);
}

// ==================== 9. FUNCIÓN renderMainChart ====================

function renderMainChart(solution, peakDay) {
  const { t, S, I, R } = solution;

  const shapes = [];
  const annotations = [];

  if (peakDay !== undefined && peakDay !== null && t.length > 0 && peakDay >= t[0] && peakDay <= t[t.length - 1]) {
    shapes.push({
      type: 'line',
      x0: peakDay,
      y0: 0,
      x1: peakDay,
      y1: 1,
      yref: 'paper',
      line: { color: '#f85149', width: 2, dash: 'dot' }
    });

    annotations.push({
      x: peakDay,
      y: 1,
      yref: 'paper',
      xref: 'x',
      text: 'Pico: día ' + Math.round(peakDay),
      showarrow: true,
      arrowhead: 2,
      ax: 40,
      ay: -30,
      font: { color: '#f85149', size: 12 },
      bgcolor: 'rgba(13, 17, 23, 0.8)',
      bordercolor: '#f85149',
      borderwidth: 1
    });
  }

  const traceS = {
    x: t, y: S,
    type: 'scatter', mode: 'lines',
    name: 'Susceptibles S(t)',
    line: { color: COLORS.susceptible, width: 2.5 },
    fill: 'tozeroy',
    fillcolor: 'rgba(88, 166, 255, 0.08)',
    hovertemplate: 'Día %{x}<br>S = %{y:,.0f}<extra></extra>'
  };

  const traceI = {
    x: t, y: I,
    type: 'scatter', mode: 'lines',
    name: 'Infectados I(t)',
    line: { color: COLORS.infectado, width: 2.5 },
    fill: 'tozeroy',
    fillcolor: 'rgba(248, 81, 73, 0.15)',
    hovertemplate: 'Día %{x}<br>I = %{y:,.0f}<extra></extra>'
  };

  const traceR = {
    x: t, y: R,
    type: 'scatter', mode: 'lines',
    name: 'Recuperados R(t)',
    line: { color: COLORS.recuperado, width: 2.5 },
    fill: 'tozeroy',
    fillcolor: 'rgba(63, 185, 80, 0.08)',
    hovertemplate: 'Día %{x}<br>R = %{y:,.0f}<extra></extra>'
  };

  const layout = {
    title: { text: 'Evolución de la epidemia SIR', font: { color: '#e6edf3', size: 16 } },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b949e', family: 'Segoe UI, system-ui, sans-serif' },
    xaxis: {
      title: { text: 'Días', font: { color: '#8b949e' } },
      gridcolor: '#21262d',
      zerolinecolor: '#30363d',
      tickfont: { color: '#8b949e' }
    },
    yaxis: {
      title: { text: 'Número de personas', font: { color: '#8b949e' } },
      gridcolor: '#21262d',
      zerolinecolor: '#30363d',
      tickfont: { color: '#8b949e' }
    },
    legend: {
      font: { color: '#e6edf3', size: 12 },
      bgcolor: 'rgba(22, 27, 34, 0.8)',
      bordercolor: '#30363d',
      borderwidth: 1
    },
    hovermode: 'x unified',
    shapes: shapes,
    annotations: annotations,
    margin: { t: 50, r: 20, b: 50, l: 60 }
  };

  const config = {
    responsive: true,
    displayModeBar: false
  };

  Plotly.newPlot('mainChart', [traceS, traceI, traceR], layout, config);
}

// ==================== 10. FUNCIÓN renderNumericalComparisonChart ====================

function renderNumericalComparisonChart(beta, gamma, N, S0, I0, R0_init, days, dt) {
  const solEuler = solveEuler(beta, gamma, N, S0, I0, R0_init, days, dt);
  const solHeun = solveHeun(beta, gamma, N, S0, I0, R0_init, days, dt);
  const solRK4 = solveRK4(beta, gamma, N, S0, I0, R0_init, days, dt);

  const traceEuler = {
    x: solEuler.t, y: solEuler.I,
    type: 'scatter', mode: 'lines',
    name: 'Euler',
    line: { color: COLORS.euler, width: 2, dash: 'dot' },
    fill: 'tozeroy',
    fillcolor: 'rgba(240, 136, 62, 0.06)',
    hovertemplate: 'Día %{x}<br>I_Euler = %{y:,.2f}<extra></extra>'
  };

  const traceHeun = {
    x: solHeun.t, y: solHeun.I,
    type: 'scatter', mode: 'lines',
    name: 'Heun',
    line: { color: COLORS.heun, width: 2, dash: 'dash' },
    fill: 'tozeroy',
    fillcolor: 'rgba(210, 168, 255, 0.06)',
    hovertemplate: 'Día %{x}<br>I_Heun = %{y:,.2f}<extra></extra>'
  };

  const traceRK4 = {
    x: solRK4.t, y: solRK4.I,
    type: 'scatter', mode: 'lines',
    name: 'RK4',
    line: { color: COLORS.rk4, width: 2.5 },
    fill: 'tozeroy',
    fillcolor: 'rgba(88, 166, 255, 0.08)',
    hovertemplate: 'Día %{x}<br>I_RK4 = %{y:,.2f}<extra></extra>'
  };

  const layout = {
    title: { text: 'Comparación de métodos numéricos — I(t)', font: { color: '#e6edf3', size: 16 } },
    paper_bgcolor: 'rgba(0,0,0,0)',
    plot_bgcolor: 'rgba(0,0,0,0)',
    font: { color: '#8b949e', family: 'Segoe UI, system-ui, sans-serif' },
    xaxis: {
      title: { text: 'Días', font: { color: '#8b949e' } },
      gridcolor: '#21262d',
      zerolinecolor: '#30363d',
      tickfont: { color: '#8b949e' }
    },
    yaxis: {
      title: { text: 'Infectados I(t)', font: { color: '#8b949e' } },
      gridcolor: '#21262d',
      zerolinecolor: '#30363d',
      tickfont: { color: '#8b949e' }
    },
    legend: {
      font: { color: '#e6edf3', size: 12 },
      bgcolor: 'rgba(22, 27, 34, 0.8)',
      bordercolor: '#30363d',
      borderwidth: 1
    },
    hovermode: 'x unified',
    margin: { t: 50, r: 20, b: 50, l: 60 }
  };

  const config = {
    responsive: true,
    displayModeBar: false
  };

  Plotly.newPlot('numericalChart', [traceEuler, traceHeun, traceRK4], layout, config);

  // Compute relative errors
  let maxErrEuler = 0;
  let maxErrHeun = 0;
  for (let i = 0; i < solRK4.I.length; i++) {
    const ref = solRK4.I[i];
    if (ref > 0) {
      const errE = Math.abs(solEuler.I[i] - ref) / ref;
      const errH = Math.abs(solHeun.I[i] - ref) / ref;
      if (errE > maxErrEuler) maxErrEuler = errE;
      if (errH > maxErrHeun) maxErrHeun = errH;
    }
  }

  document.getElementById('err-euler').textContent = (maxErrEuler * 100).toFixed(4) + '%';
  document.getElementById('err-heun').textContent = (maxErrHeun * 100).toFixed(4) + '%';
}

// ==================== 11. FUNCIÓN PRINCIPAL runSimulation ====================

function runSimulation() {
  const N = appState.N;
  const I0 = appState.I0;
  const beta = appState.beta;
  const gamma = appState.gamma;
  const pctVac = appState.vac;
  const days = appState.days;
  const dt = appState.dt;

  const vacunados = Math.floor(N * pctVac / 100);
  const S0 = N - I0 - vacunados;
  const R0_init = vacunados;

  const solution = solveRK4(beta, gamma, N, S0, I0, R0_init, days, dt);
  const stats = computeStats(solution);
  const r0 = updateR0Display(beta, gamma);

  window._lastSolution = solution;

  renderMainChart(solution, stats.peakDay);
  resetHighlight();
  renderNumericalComparisonChart(beta, gamma, N, S0, I0, R0_init, days, dt);

  updateDerivativesLive(S0, I0, R0_init, N, beta, gamma);

  document.getElementById('stat-r0').textContent = r0.toFixed(2);
  document.getElementById('stat-peak-day').textContent = Math.round(stats.peakDay);
  document.getElementById('stat-peak-val').textContent = Math.round(stats.peakVal);
  document.getElementById('stat-recovered').textContent = Math.round(stats.totalRecovered);
  document.getElementById('stat-affected').textContent = stats.pctAffected.toFixed(1) + '%';

  if (stepState.active) {
    stepState.currentIdx = 0;
    updateStepDisplay();
    renderStepCursor();
  }
}

// ==================== 12. FUNCIÓN loadScenario ====================

function loadScenario(scenarioKey) {
  const preset = scenarioPresets[scenarioKey];
  if (!preset) return;

  const mappings = {
    N: { slider: 'slider-N', val: 'val-N', fmt: (v) => String(v) },
    I0: { slider: 'slider-I0', val: 'val-I0', fmt: (v) => String(v) },
    beta: { slider: 'slider-beta', val: 'val-beta', fmt: (v) => v.toFixed(2) },
    gamma: { slider: 'slider-gamma', val: 'val-gamma', fmt: (v) => v.toFixed(2) },
    vac: { slider: 'slider-vac', val: 'val-vac', fmt: (v) => v + '%' },
    days: { slider: 'slider-days', val: 'val-days', fmt: (v) => String(v) },
    dt: { slider: 'slider-dt', val: 'val-dt', fmt: (v) => v.toFixed(1) }
  };

  Object.keys(mappings).forEach(key => {
    const m = mappings[key];
    const value = preset[key];
    document.getElementById(m.slider).value = value;
    document.getElementById(m.val).textContent = m.fmt(value);
    appState[key] = value;
  });

  // Mark active scenario card
  document.querySelectorAll('.scenario-card').forEach(card => {
    card.classList.remove('active');
  });
  const activeCard = document.querySelector('.scenario-card[data-scenario="' + scenarioKey + '"]');
  if (activeCard) activeCard.classList.add('active');

  // Turn off step mode if active
  if (stepState.active) {
    document.getElementById('btn-step-toggle').checked = false;
    toggleStepMode();
  }

  resetHighlight();

  runSimulation();
}

// ==================== 13. EVENT LISTENERS ====================

function setupEventListeners() {
  // Slider input events
  const sliderConfig = [
    { id: 'slider-N', key: 'N', format: (v) => String(v) },
    { id: 'slider-I0', key: 'I0', format: (v) => String(v) },
    { id: 'slider-beta', key: 'beta', format: (v) => parseFloat(v).toFixed(2) },
    { id: 'slider-gamma', key: 'gamma', format: (v) => parseFloat(v).toFixed(2) },
    { id: 'slider-vac', key: 'vac', format: (v) => v + '%' },
    { id: 'slider-days', key: 'days', format: (v) => String(v) },
    { id: 'slider-dt', key: 'dt', format: (v) => parseFloat(v).toFixed(1) }
  ];

  sliderConfig.forEach(cfg => {
    const slider = document.getElementById(cfg.id);
    slider.addEventListener('input', function () {
      const val = this.value;
      const numericVal = cfg.key === 'beta' || cfg.key === 'gamma' || cfg.key === 'dt' ? parseFloat(val) : Number(val);
      appState[cfg.key] = numericVal;

      // Map slider id to value display id
      const valMap = {
        'slider-N': 'val-N',
        'slider-I0': 'val-I0',
        'slider-beta': 'val-beta',
        'slider-gamma': 'val-gamma',
        'slider-vac': 'val-vac',
        'slider-days': 'val-days',
        'slider-dt': 'val-dt'
      };
      document.getElementById(valMap[cfg.id]).textContent = cfg.format(val);

      // Clear active scenario card
      document.querySelectorAll('.scenario-card').forEach(card => card.classList.remove('active'));

      runSimulation();
    });
  });

  // Run button
  document.getElementById('btn-run').addEventListener('click', function () {
    this.style.transform = 'scale(0.95)';
    setTimeout(() => { this.style.transform = ''; }, 150);
    runSimulation();
  });

  // Reset button
  document.getElementById('btn-reset').addEventListener('click', function () {
    const defaults = { N: 10000, I0: 5, beta: 0.30, gamma: 0.10, vac: 0, days: 160, dt: 1.0 };
    Object.keys(defaults).forEach(key => {
      appState[key] = defaults[key];
    });

    document.getElementById('slider-N').value = 10000;
    document.getElementById('val-N').textContent = '10000';
    document.getElementById('slider-I0').value = 5;
    document.getElementById('val-I0').textContent = '5';
    document.getElementById('slider-beta').value = '0.30';
    document.getElementById('val-beta').textContent = '0.30';
    document.getElementById('slider-gamma').value = '0.10';
    document.getElementById('val-gamma').textContent = '0.10';
    document.getElementById('slider-vac').value = 0;
    document.getElementById('val-vac').textContent = '0%';
    document.getElementById('slider-days').value = 160;
    document.getElementById('val-days').textContent = '160';
    document.getElementById('slider-dt').value = '1.0';
    document.getElementById('val-dt').textContent = '1.0';

    document.querySelectorAll('.scenario-card').forEach(card => card.classList.remove('active'));

    runSimulation();
  });

  // Scenario cards
  document.querySelectorAll('.scenario-card').forEach(card => {
    card.addEventListener('click', function () {
      const key = this.getAttribute('data-scenario');
      loadScenario(key);
    });
  });

  // Step mode toggle
  document.getElementById('btn-step-toggle').addEventListener('change', function () {
    toggleStepMode();
  });

  // Step navigation
  document.getElementById('btn-step-next').addEventListener('click', stepNext);
  document.getElementById('btn-step-prev').addEventListener('click', stepPrev);

  // Equation click handlers
  setupEquationClickHandlers();
}

// ==================== 14. ECUACIONES CLIQUEABLES (FEATURE 6) ====================

let _activeCurve = null;

function resetHighlight() {
  _activeCurve = null;
  document.querySelectorAll('.equation').forEach(function(e) { e.classList.remove('active'); });
  Plotly.restyle('mainChart', { opacity: 1, 'line.width': 2.5 });
}

function setupEquationClickHandlers() {
  document.querySelectorAll('.equation').forEach(function(eq) {
    eq.addEventListener('click', function() {
      const curve = this.getAttribute('data-curve');

      if (_activeCurve === curve) {
        resetHighlight();
        return;
      }

      _activeCurve = curve;
      document.querySelectorAll('.equation').forEach(function(e) { e.classList.remove('active'); });
      this.classList.add('active');

      const traceMap = { S: 0, I: 1, R: 2 };
      const activeIdx = traceMap[curve];

      const opacityVals = [0.15, 0.15, 0.15];
      const widthVals = [1.5, 1.5, 1.5];
      opacityVals[activeIdx] = 1;
      widthVals[activeIdx] = 3.5;

      Plotly.restyle('mainChart', 'opacity', opacityVals);
      Plotly.restyle('mainChart', 'line.width', widthVals);
    });
  });
}

// ==================== 15. INICIALIZACIÓN DE PESTAÑAS (NAVBAR) ====================

function initTabs() {
  document.querySelectorAll('.tab-btn').forEach(function(btn) {
    btn.addEventListener('click', function() {
      const tabId = this.getAttribute('data-tab');

      document.querySelectorAll('.tab-btn').forEach(function(b) { b.classList.remove('active'); });
      document.querySelectorAll('.tab-content').forEach(function(c) { c.classList.remove('active'); });

      this.classList.add('active');
      document.getElementById('tab-' + tabId).classList.add('active');
    });
  });
}

// ==================== 16. INICIALIZACIÓN ====================

document.addEventListener('DOMContentLoaded', function () {
  // Set initial slider values from appState
  document.getElementById('slider-N').value = appState.N;
  document.getElementById('val-N').textContent = appState.N;
  document.getElementById('slider-I0').value = appState.I0;
  document.getElementById('val-I0').textContent = appState.I0;
  document.getElementById('slider-beta').value = appState.beta;
  document.getElementById('val-beta').textContent = appState.beta.toFixed(2);
  document.getElementById('slider-gamma').value = appState.gamma;
  document.getElementById('val-gamma').textContent = appState.gamma.toFixed(2);
  document.getElementById('slider-vac').value = appState.vac;
  document.getElementById('val-vac').textContent = appState.vac + '%';
  document.getElementById('slider-days').value = appState.days;
  document.getElementById('val-days').textContent = appState.days;
  document.getElementById('slider-dt').value = appState.dt;
  document.getElementById('val-dt').textContent = appState.dt.toFixed(1);

  initTabs();
  setupEventListeners();
  runSimulation();
});
