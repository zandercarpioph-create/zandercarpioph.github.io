/* ===========================================================================
   Portfolio behaviour: theme, reveals, counters, and two live visualisations.
   No dependencies. The epidemic lab integrates the same systems shown in the
   R figures, so the page is a working demo rather than a picture of one.
   =========================================================================== */
(function () {
  'use strict';

  var $  = function (s, r) { return (r || document).querySelector(s); };
  var $$ = function (s, r) { return Array.prototype.slice.call((r || document).querySelectorAll(s)); };

  /* ------------------------------------------------------------- theme --- */
  var root = document.documentElement;
  try {
    var saved = localStorage.getItem('theme');
    if (saved === 'light' || saved === 'dark') root.setAttribute('data-theme', saved);
  } catch (e) { /* private mode, stay on the system preference */ }

  var toggle = $('#theme-toggle');
  if (toggle) {
    toggle.addEventListener('click', function () {
      var dark = root.getAttribute('data-theme') === 'dark' ||
        (!root.hasAttribute('data-theme') &&
          window.matchMedia('(prefers-color-scheme: dark)').matches);
      var next = dark ? 'light' : 'dark';
      root.setAttribute('data-theme', next);
      try { localStorage.setItem('theme', next); } catch (e) {}
      redrawAll();
    });
  }

  /* --------------------------------------------------------------- nav --- */
  var nav = $('.nav');
  var onScroll = function () {
    if (nav) nav.classList.toggle('stuck', window.scrollY > 8);
  };
  window.addEventListener('scroll', onScroll, { passive: true });
  onScroll();

  /* ----------------------------------------------------------- reveals --- */
  if ('IntersectionObserver' in window) {
    var io = new IntersectionObserver(function (entries) {
      entries.forEach(function (en) {
        if (en.isIntersecting) { en.target.classList.add('in'); io.unobserve(en.target); }
      });
    }, { rootMargin: '0px 0px -8% 0px', threshold: .06 });
    $$('.reveal').forEach(function (el) { io.observe(el); });
  } else {
    $$('.reveal').forEach(function (el) { el.classList.add('in'); });
  }

  /* ---------------------------------------------------------- counters --- */
  function runCounter(el) {
    var target = parseFloat(el.getAttribute('data-count'));
    var suffix = el.getAttribute('data-suffix') || '';
    var dur = 1100, t0 = null;
    function step(ts) {
      if (t0 === null) t0 = ts;
      var p = Math.min((ts - t0) / dur, 1);
      var eased = 1 - Math.pow(1 - p, 3);
      var v = target * eased;
      el.textContent = (target >= 100 ? Math.round(v) : v.toFixed(v < 10 && target % 1 ? 1 : 0)) + suffix;
      if (p < 1) requestAnimationFrame(step);
    }
    requestAnimationFrame(step);
  }
  if ('IntersectionObserver' in window) {
    var cio = new IntersectionObserver(function (entries) {
      entries.forEach(function (en) {
        if (en.isIntersecting) { runCounter(en.target); cio.unobserve(en.target); }
      });
    }, { threshold: .5 });
    $$('[data-count]').forEach(function (el) { cio.observe(el); });
  }

  /* ------------------------------------------------------- svg helpers --- */
  var NS = 'http://www.w3.org/2000/svg';
  function el(name, attrs) {
    var n = document.createElementNS(NS, name);
    for (var k in attrs) if (attrs[k] !== null && attrs[k] !== undefined) n.setAttribute(k, attrs[k]);
    return n;
  }
  function css(name) {
    return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
  }
  function pathFrom(pts) {
    return pts.map(function (p, i) { return (i ? 'L' : 'M') + p[0].toFixed(2) + ' ' + p[1].toFixed(2); }).join(' ');
  }

  /* ==========================================================================
     Epidemic lab: RK4 on SIR or SEIR, redrawn on every slider move.
     ========================================================================== */
  var LAB = {
    host: $('#lab'),
    model: 'sir',
    N: 1000,
    T: 160,
    dt: 0.5
  };

  function deriv(model, y, p) {
    var S = y[0];
    if (model === 'sir') {
      var I = y[1];
      return [-p.beta * S * I, p.beta * S * I - p.gamma * I, p.gamma * I];
    }
    var E = y[1], I2 = y[2];
    return [
      -p.beta * S * E,
      p.beta * S * E - p.sigma * E,
      p.sigma * E - p.gamma * I2,
      p.gamma * I2
    ];
  }

  function rk4(model, p) {
    var y = model === 'sir' ? [LAB.N - 1, 1, 0] : [LAB.N - 1, 1, 0, 0];
    var out = [y.slice()], h = LAB.dt, steps = Math.round(LAB.T / h);
    for (var s = 0; s < steps; s++) {
      var k1 = deriv(model, y, p);
      var y2 = y.map(function (v, i) { return v + h / 2 * k1[i]; });
      var k2 = deriv(model, y2, p);
      var y3 = y.map(function (v, i) { return v + h / 2 * k2[i]; });
      var k3 = deriv(model, y3, p);
      var y4 = y.map(function (v, i) { return v + h * k3[i]; });
      var k4 = deriv(model, y4, p);
      y = y.map(function (v, i) { return v + h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]); });
      out.push(y.slice());
    }
    return out;
  }

  function drawLab() {
    if (!LAB.host) return;
    var plot = $('#lab-plot');
    var p = {
      beta:  parseFloat($('#p-beta').value),
      gamma: parseFloat($('#p-gamma').value),
      sigma: parseFloat($('#p-sigma').value)
    };
    $('#v-beta').textContent  = p.beta.toFixed(5);
    $('#v-gamma').textContent = p.gamma.toFixed(3);
    $('#v-sigma').textContent = p.sigma.toFixed(2);

    var series = rk4(LAB.model, p);
    var W = 760, H = 300, m = { t: 16, r: 16, b: 30, l: 46 };
    var iw = W - m.l - m.r, ih = H - m.t - m.b;
    var x = function (d) { return m.l + (d / LAB.T) * iw; };
    var y = function (v) { return m.t + ih - (v / LAB.N) * ih; };

    var svg = el('svg', { viewBox: '0 0 ' + W + ' ' + H, role: 'img',
      'aria-label': 'Simulated epidemic curves for the selected model and parameters' });

    // grid
    for (var g = 0; g <= 4; g++) {
      var gy = m.t + (ih / 4) * g;
      svg.appendChild(el('line', { x1: m.l, x2: W - m.r, y1: gy, y2: gy,
        stroke: css('--line'), 'stroke-width': 1 }));
      var lab = el('text', { x: m.l - 9, y: gy + 4, 'text-anchor': 'end',
        fill: css('--faint'), 'font-size': 10, 'font-family': 'var(--mono)' });
      lab.textContent = Math.round(LAB.N - (LAB.N / 4) * g);
      svg.appendChild(lab);
    }
    for (var d = 0; d <= LAB.T; d += 40) {
      var tx = el('text', { x: x(d), y: H - 9, 'text-anchor': 'middle',
        fill: css('--faint'), 'font-size': 10, 'font-family': 'var(--mono)' });
      tx.textContent = d;
      svg.appendChild(tx);
    }

    var names = LAB.model === 'sir'
      ? [['Susceptible', '--accent-3', 0], ['Infected', '--accent', 1], ['Recovered', '--accent-2', 2]]
      : [['Susceptible', '--accent-3', 0], ['Exposed', '--faint', 1], ['Infected', '--accent', 2], ['Recovered', '--accent-2', 3]];

    var infIdx = LAB.model === 'sir' ? 1 : 2;

    // shaded area under the infected curve
    var areaPts = series.map(function (row, i) { return [x(i * LAB.dt), y(row[infIdx])]; });
    var area = 'M' + x(0) + ' ' + y(0) + ' ' + areaPts.map(function (q) {
      return 'L' + q[0].toFixed(2) + ' ' + q[1].toFixed(2); }).join(' ') +
      ' L' + x(LAB.T) + ' ' + y(0) + ' Z';
    svg.appendChild(el('path', { d: area, fill: css('--accent'), opacity: .1 }));

    names.forEach(function (n) {
      var pts = series.map(function (row, i) { return [x(i * LAB.dt), y(row[n[2]])]; });
      svg.appendChild(el('path', { d: pathFrom(pts), fill: 'none',
        stroke: css(n[1]), 'stroke-width': 2.1, 'stroke-linejoin': 'round' }));
    });

    plot.innerHTML = '';
    plot.appendChild(svg);

    // legend
    var lg = $('#lab-legend');
    lg.innerHTML = '';
    names.forEach(function (n) {
      var s = document.createElement('span');
      var i = document.createElement('i');
      i.style.background = css(n[1]);
      s.appendChild(i);
      s.appendChild(document.createTextNode(n[0]));
      lg.appendChild(s);
    });

    // readout
    var peak = 0, peakDay = 0;
    series.forEach(function (row, i) {
      if (row[infIdx] > peak) { peak = row[infIdx]; peakDay = i * LAB.dt; }
    });
    var r0 = p.beta * LAB.N / p.gamma;
    $('#r-r0').textContent = r0.toFixed(2);
    $('#r-peak').textContent = Math.round(peak);
    $('#r-day').textContent = Math.round(peakDay);
    var finalR = series[series.length - 1][LAB.model === 'sir' ? 2 : 3];
    $('#r-attack').textContent = Math.round(finalR / LAB.N * 100) + '%';
  }

  if (LAB.host) {
    $$('#lab .seg button').forEach(function (b) {
      b.addEventListener('click', function () {
        LAB.model = b.getAttribute('data-model');
        $$('#lab .seg button').forEach(function (o) {
          o.setAttribute('aria-pressed', String(o === b));
        });
        $('#ctl-sigma').style.visibility = LAB.model === 'seir' ? 'visible' : 'hidden';
        drawLab();
      });
    });
    ['#p-beta', '#p-gamma', '#p-sigma'].forEach(function (id) {
      $(id).addEventListener('input', drawLab);
    });
  }

  /* ==========================================================================
     Dengue chart: two views over the aggregates exported by make_figures.R
     ========================================================================== */
  var DEN = { data: null, view: 'weeks', host: $('#dengue-plot') };

  function drawDengue() {
    if (!DEN.data || !DEN.host) return;
    var d = DEN.data;
    var W = 760, H = 320, m = { t: 18, r: 18, b: 78, l: 52 };
    var iw = W - m.l - m.r, ih = H - m.t - m.b;
    var svg = el('svg', { viewBox: '0 0 ' + W + ' ' + H, role: 'img' });
    var tip = $('#dengue-tip');

    var labels, values, unit;
    if (DEN.view === 'weeks') {
      labels = d.weeks; values = d.cases; unit = 'cases';
      svg.setAttribute('aria-label', 'Reported dengue cases by surveillance week');
    } else {
      labels = d.barangays; values = d.means; unit = 'mean cases/week';
      svg.setAttribute('aria-label', 'Mean weekly dengue cases by barangay');
    }

    var max = Math.max.apply(null, values) * 1.08;
    // Space points by their real date where offsets are supplied, so a
    // four-month hole in the record reads as a hole and not as a straight line.
    var days = (DEN.view === 'weeks' && d.days && d.days.length === values.length) ? d.days : null;
    var span = days ? (days[days.length - 1] - days[0]) || 1 : 1;
    var x = days
      ? function (i) { return m.l + ((days[i] - days[0]) / span) * iw; }
      : function (i) { return m.l + (iw / (labels.length - 1 || 1)) * i; };
    var xb = function (i) { return m.l + (iw / labels.length) * (i + .5); };
    var y = function (v) { return m.t + ih - (v / max) * ih; };

    for (var g = 0; g <= 4; g++) {
      var gy = m.t + (ih / 4) * g;
      svg.appendChild(el('line', { x1: m.l, x2: W - m.r, y1: gy, y2: gy,
        stroke: css('--line'), 'stroke-width': 1 }));
      var t = el('text', { x: m.l - 9, y: gy + 4, 'text-anchor': 'end',
        fill: css('--faint'), 'font-size': 10, 'font-family': 'var(--mono)' });
      t.textContent = Math.round(max - (max / 4) * g);
      svg.appendChild(t);
    }

    if (DEN.view === 'weeks') {
      var pts = values.map(function (v, i) { return [x(i), y(v)]; });

      // break the series wherever more than one week separates two readings
      var runs = [[0]];
      for (var k = 1; k < values.length; k++) {
        if (days && days[k] - days[k - 1] > 7) runs.push([k]);
        else runs[runs.length - 1].push(k);
      }
      runs.forEach(function (run) {
        var rp = run.map(function (i) { return pts[i]; });
        var first = rp[0], last = rp[rp.length - 1];
        var areaD = 'M' + first[0].toFixed(2) + ' ' + y(0) + ' ' +
          rp.map(function (q) { return 'L' + q[0].toFixed(2) + ' ' + q[1].toFixed(2); }).join(' ') +
          ' L' + last[0].toFixed(2) + ' ' + y(0) + ' Z';
        svg.appendChild(el('path', { d: areaD, fill: css('--accent'), opacity: .12 }));
        if (rp.length > 1) {
          svg.appendChild(el('path', { d: pathFrom(rp), fill: 'none', stroke: css('--accent'),
            'stroke-width': 2.2, 'stroke-linejoin': 'round' }));
        }
      });

      pts.forEach(function (q, i) {
        svg.appendChild(el('circle', { cx: q[0], cy: q[1], r: 3, fill: css('--accent') }));
        var hit = el('rect', { x: q[0] - 9, y: m.t, width: 18, height: ih,
          fill: 'transparent', class: 'hitline' });
        hit.addEventListener('mouseenter', function () { showTip(q[0], q[1], labels[i], values[i], unit); });
        hit.addEventListener('mouseleave', hideTip);
        svg.appendChild(hit);
      });
    } else {
      values.forEach(function (v, i) {
        var bw = (iw / labels.length) * .68;
        var bx = xb(i) - bw / 2;
        var r = el('rect', { x: bx, y: y(v), width: bw, height: m.t + ih - y(v),
          fill: css('--accent'), opacity: .86, rx: 3 });
        r.addEventListener('mouseenter', function () { showTip(xb(i), y(v), labels[i], v, unit); });
        r.addEventListener('mouseleave', hideTip);
        svg.appendChild(r);
      });
    }

    // x labels, thinned so they never collide
    var every = DEN.view === 'weeks' ? 3 : 1;
    labels.forEach(function (l, i) {
      if (i % every) return;
      var cx = DEN.view === 'weeks' ? x(i) : xb(i);
      var short = l.length > 15 ? l.slice(0, 14) + '.' : l;
      var t = el('text', { x: cx, y: m.t + ih + 14, fill: css('--faint'),
        'font-size': 9.5, 'font-family': 'var(--mono)',
        transform: 'rotate(38 ' + cx + ' ' + (m.t + ih + 14) + ')' });
      t.textContent = short;
      svg.appendChild(t);
    });

    DEN.host.innerHTML = '';
    DEN.host.appendChild(svg);

    function showTip(px, py, label, val, u) {
      if (!tip) return;
      var box = DEN.host.getBoundingClientRect();
      tip.textContent = label + ': ' + val + ' ' + u;
      tip.style.left = (px / W * box.width) + 'px';
      tip.style.top = (py / H * box.height) + 'px';
      tip.classList.add('on');
    }
    function hideTip() { if (tip) tip.classList.remove('on'); }
  }

  $$('#dengue-seg button').forEach(function (b) {
    b.addEventListener('click', function () {
      DEN.view = b.getAttribute('data-view');
      $$('#dengue-seg button').forEach(function (o) {
        o.setAttribute('aria-pressed', String(o === b));
      });
      drawDengue();
    });
  });

  if (DEN.host) {
    fetch('assets/data/dengue.json')
      .then(function (r) { return r.json(); })
      .then(function (j) { DEN.data = j; drawDengue(); })
      .catch(function () {
        DEN.host.innerHTML = '<p style="color:var(--faint);font-size:.85rem;padding:20px">' +
          'Chart data unavailable. The static figures below carry the same information.</p>';
      });
  }

  /* ------------------------------------------------------- hero ribbon --- */
  function drawHero() {
    var host = $('#hero-curve');
    if (!host) return;
    var W = 1200, H = 260;
    var svg = el('svg', { viewBox: '0 0 ' + W + ' ' + H, preserveAspectRatio: 'none',
      'aria-hidden': 'true' });
    // an SIR infected curve, solved with the same integrator the lab uses
    var series = rk4('sir', { beta: 0.00032, gamma: 0.1, sigma: 0.2 });
    var peak = Math.max.apply(null, series.map(function (r) { return r[1]; }));
    var pts = series.map(function (r, i) {
      return [(i / (series.length - 1)) * W, H - (r[1] / peak) * (H * .82) - 10];
    });
    var d = pathFrom(pts);
    svg.appendChild(el('path', { class: 'fill', d: d + ' L' + W + ' ' + H + ' L0 ' + H + ' Z',
      fill: css('--accent'), stroke: 'none' }));
    svg.appendChild(el('path', { class: 'stroke', d: d }));
    host.innerHTML = '';
    host.appendChild(svg);
  }

  function redrawAll() { drawHero(); drawLab(); drawDengue(); }

  /* ---------------------------------------------------------- kickoff --- */
  drawHero();
  drawLab();
  window.addEventListener('resize', function () {
    clearTimeout(window.__rt);
    window.__rt = setTimeout(drawDengue, 160);
  });
})();
