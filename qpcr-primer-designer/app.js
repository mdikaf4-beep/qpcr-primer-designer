/* qPCR Primer Designer by DF — UI wiring (vanilla JS) */
(function () {
  'use strict';

  const E = window.PrimerEngine;
  const $ = sel => document.querySelector(sel);
  const $$ = sel => Array.from(document.querySelectorAll(sel));

  let currentResult = null;
  let currentPair = null;
  let sortKey = 'rank';
  let sortDir = 1;
  let evalPair = null;

  function esc(s) {
    return String(s).replace(/[<>&"']/g, ch => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;', '"': '&quot;', "'": '&#39;' }[ch]));
  }

  function num(id) {
    const v = parseFloat($(id).value);
    return isNaN(v) ? 0 : v;
  }

  function showError(msg) {
    const box = $('#errorBox');
    if (!msg) { box.classList.add('hidden'); box.textContent = ''; return; }
    box.textContent = msg;
    box.classList.remove('hidden');
  }

  function makeExampleSequence(n) {
    let seed = 42;
    const out = [];
    for (let i = 0; i < (n || 600); i++) {
      seed = (seed * 1664525 + 1013904223) >>> 0;
      const r = seed / 4294967296;
      out.push(r < 0.25 ? 'A' : r < 0.5 ? 'T' : r < 0.75 ? 'G' : 'C');
    }
    return out.join('');
  }

  function loadExample() {
    $('#sequenceInput').value = '>synthetic_target_600bp\n' + makeExampleSequence(600);
    $('#snpText').value = '';
    $('#maskText').value = '';
    $('#junctionText').value = '';
    $('#offtargetText').value = '';
    showError('');
  }

  function getConfig() {
    return {
      assayType: $('#assayType').value,
      minLen: num('#minLen'),
      maxLen: num('#maxLen'),
      hardTmMin: num('#tmMin'),
      hardTmMax: num('#tmMax'),
      maxDeltaTm: num('#deltaTm'),
      maxDeltaTmHard: Math.max(5, num('#deltaTm') + 3),
      gcMin: num('#gcMin'),
      gcMax: num('#gcMax'),
      ampliconMin: num('#ampMin'),
      ampliconMax: num('#ampMax'),
      topN: num('#topN') || 10,
      snps: E.parsePositions($('#snpText').value),
      mask: E.parseIntervals($('#maskText').value),
      junctions: E.parsePositions($('#junctionText').value),
      requireJunction: $('#requireJunction').checked,
      offtarget: E.parseFasta($('#offtargetText').value)
    };
  }

  function readSequence() {
    const recs = E.parseFasta($('#sequenceInput').value);
    if (!recs.length) throw new Error('Sekuens kosong atau tidak valid.');
    if (recs.length > 1) showError(`Ada ${recs.length} record FASTA; hanya record pertama yang dipakai.`);
    return recs[0];
  }

  function runDesign() {
    showError('');
    try {
      const rec = readSequence();
      const cfg = getConfig();
      const t0 = performance.now();
      currentResult = E.design(rec.seq, cfg);
      currentResult._elapsedUi = Math.round(performance.now() - t0);
      currentResult._targetName = rec.name;
      currentPair = currentResult.pairs[0] || null;
      renderSummary();
      renderTable();
      renderDetail();
      switchTab('design');
    } catch (e) {
      showError('Error: ' + e.message);
    }
  }

  function renderSummary() {
    const s = currentResult.summary;
    const m = currentResult.meta;
    const html = [
      stat(s.forwardCandidates, 'Kandidat FWD'),
      stat(s.reverseCandidates, 'Kandidat REV'),
      stat(s.pairsBeforeScoring, 'Pasangan lolos filter'),
      stat(s.pairsReturned, `Pasangan ditampilkan (top ${m.config.topN})`),
      stat((m.elapsedMs || 0) + ' ms', 'Waktu engine')
    ].join('');
    $('#summaryBar').innerHTML = html;
  }

  function stat(value, label) {
    return `<div class="stat"><b>${esc(value)}</b><span>${esc(label)}</span></div>`;
  }

  function sortedPairs() {
    const arr = currentResult.pairs.slice();
    arr.sort((a, b) => {
      let va, vb;
      switch (sortKey) {
        case 'score': va = a.score; vb = b.score; break;
        case 'fwdTm': va = a.fwd.tm; vb = b.fwd.tm; break;
        case 'gc': va = a.fwd.gc; vb = b.fwd.gc; break;
        case 'amp': va = a.ampliconLength; vb = b.ampliconLength; break;
        default: va = a.rank; vb = b.rank;
      }
      return (va - vb) * sortDir;
    });
    return arr;
  }

  function qcPills(qc) {
    return qc.map(q => {
      const cls = q.severity === 'FAIL' ? 'fail' : q.severity === 'WARN' ? 'warn' : 'ok';
      return `<span class="pill ${cls}" title="${esc(q.message || '')}">${esc(q.code)}</span>`;
    }).join(' ');
  }

  function renderTable() {
    const body = $('#resultsBody');
    if (!currentResult.pairs.length) {
      body.innerHTML = '<tr><td colspan="8">Tidak ada pasangan yang lolos filter. Longgarkan Tm/GC atau ukuran amplicon.</td></tr>';
      return;
    }
    body.innerHTML = sortedPairs().map(p => `
      <tr data-rank="${p.rank}" class="${currentPair && currentPair.rank === p.rank ? 'selected' : ''}">
        <td>${p.rank}</td>
        <td><b>${p.score}</b></td>
        <td class="mono">${esc(p.fwd.sequence)}</td>
        <td class="mono">${esc(p.rev.sequence)}</td>
        <td>${p.fwd.tm} / ${p.rev.tm}</td>
        <td>${p.fwd.gc} / ${p.rev.gc}</td>
        <td>${p.ampliconLength} bp</td>
        <td>${qcPills(p.qc)}</td>
      </tr>`).join('');

    $$('#resultsBody tr').forEach(tr => {
      tr.addEventListener('click', () => {
        const rank = parseInt(tr.dataset.rank, 10);
        currentPair = currentResult.pairs.find(p => p.rank === rank);
        $$('#resultsBody tr').forEach(x => x.classList.remove('selected'));
        tr.classList.add('selected');
        renderDetail();
      });
    });
  }

  function recommendProbe(pair, cfg) {
    if (cfg.assayType !== 'TaqMan') return null;
    const s = pair.amplicon;
    for (let len = 20; len <= 30; len++) {
      for (let i = 3; i + len <= s.length - 3; i++) {
        const sub = s.slice(i, i + len);
        if (/N/.test(sub)) continue;
        if (sub[0] === 'G') continue;
        const tm = E.nearestNeighborTm(sub);
        if (tm >= 65 && tm <= 70) {
          const gc = E.gcPercent(sub);
          if (gc >= 30 && gc <= 80) {
            return { sequence: sub, start: i, length: len, tm, gc };
          }
        }
      }
    }
    return null;
  }

  function detailHTML(pair, cfg) {
    const probe = recommendProbe(pair, cfg);
    const fStart = pair.fwd.start + 1;
    const fEnd = pair.fwd.end + 1;
    const rStart = pair.rev.start + 1;
    const rEnd = pair.rev.end + 1;
    const kv = `
      <div class="kv">
        <div><span>Peringkat</span><b>#${pair.rank || 1}</b></div>
        <div><span>Skor</span><b>${pair.score}</b></div>
        <div><span>Tm F / R</span><b>${pair.fwd.tm} / ${pair.rev.tm} °C</b></div>
        <div><span>ΔTm</span><b>${pair.deltaTm} °C</b></div>
        <div><span>GC F / R</span><b>${pair.fwd.gc} / ${pair.rev.gc}%</b></div>
        <div><span>Amplicon</span><b>${pair.ampliconLength} bp (${pair.ampliconGc}%)</b></div>
        <div><span>FWD posisi</span><b>${fStart}–${fEnd}</b></div>
        <div><span>REV posisi</span><b>${rStart}–${rEnd}</b></div>
        <div><span>Cross-dimer ΔG</span><b>${pair.crossDimer} kcal/mol</b></div>
        <div><span>3' dimer ΔG</span><b>${pair.threePrimeDimer} kcal/mol</b></div>
      </div>`;

    const map = `
<pre>5' ...[${esc(pair.fwd.sequence)}]..........[...].........[${esc(pair.rev.sequence)}] ... 3'
      FWD ${fStart}-${fEnd}      amplicon ${pair.ampliconLength} bp      REV ${rStart}-${rEnd}

Amplicon (5'→3'):
${esc(pair.amplicon)}
${probe ? `\nProbe (TaqMan): 5'-${esc(probe.sequence)}-3'  Tm ${probe.tm} °C  GC ${probe.gc}%  @${probe.start + 1}` : ''}</pre>`;

    const breakdown = pair.breakdown.length
      ? `<table><thead><tr><th>Rule</th><th>Penalty</th><th>Detail</th></tr></thead><tbody>` +
        pair.breakdown.map(b => `<tr><td class="mono">${esc(b.code)}</td><td>${b.penalty}</td><td>${esc(b.detail)}</td></tr>`).join('') +
        `</tbody></table>`
      : '<div class="hint">Tidak ada penalti — semua parameter berada di zona ideal.</div>';

    const qcList = pair.qc.map(q =>
      `<span class="pill ${q.severity === 'FAIL' ? 'fail' : q.severity === 'WARN' ? 'warn' : 'ok'}" title="${esc(q.message || '')}">${esc(q.code)}: ${esc(q.message || '')}</span>`
    ).join(' ');

    return kv + map +
      `<h3>QC</h3><div class="row">${qcList}</div>` +
      `<h3>Score breakdown</h3>${breakdown}`;
  }

  function renderDetail() {
    const card = $('#detailCard');
    if (!currentPair) { card.classList.add('hidden'); return; }
    const cfg = getConfig();
    card.innerHTML = `<h3>Detail pasangan #${currentPair.rank || 1}</h3>` + detailHTML(currentPair, cfg);
    card.classList.remove('hidden');
  }

  function switchTab(which) {
    const design = which === 'design';
    $('#tabDesign').classList.toggle('active', design);
    $('#tabEval').classList.toggle('active', !design);
    $('#designView').classList.toggle('hidden', !design);
    $('#evalView').classList.toggle('hidden', design);
  }

  function runEval() {
    showError('');
    try {
      const rec = readSequence();
      const pair = E.evaluatePair(rec.seq, $('#evalFwd').value, $('#evalRev').value, getConfig());
      evalPair = pair;
      $('#evalResult').innerHTML = `<h3>Evaluasi pasangan</h3>` + detailHTML(pair, getConfig());
      switchTab('eval');
    } catch (e) {
      showError('Error evaluasi: ' + e.message);
    }
  }

  function download(filename, text, type) {
    const blob = new Blob([text], { type: type || 'text/plain;charset=utf-8' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    a.remove();
    URL.revokeObjectURL(a.href);
  }

  function requireResult() {
    if (!currentResult || !currentResult.pairs.length) {
      showError('Jalankan desain terlebih dahulu.');
      return false;
    }
    return true;
  }

  function bind() {
    $('#versionBadge').textContent = 'engine v' + E.VERSION;
    loadExample();

    $('#exampleBtn').addEventListener('click', loadExample);
    $('#fileInput').addEventListener('change', e => {
      const file = e.target.files[0];
      if (!file) return;
      const reader = new FileReader();
      reader.onload = ev => { $('#sequenceInput').value = ev.target.result; showError(''); };
      reader.readAsText(file);
    });
    $('#assayType').addEventListener('change', () => {
      $('#probeHint').textContent = $('#assayType').value === 'TaqMan'
        ? 'Probe direkomendasikan dari dalam amplicon (Tm 65–70 °C, tidak overlap primer).'
        : 'Mode SYBR: hanya pasangan primer yang dihitung.';
    });

    $('#runBtn').addEventListener('click', runDesign);
    $('#clearBtn').addEventListener('click', () => {
      $('#sequenceInput').value = '';
      currentResult = null; currentPair = null; evalPair = null;
      $('#resultsBody').innerHTML = '<tr><td colspan="8">Belum ada hasil.</td></tr>';
      $('#summaryBar').innerHTML = '';
      $('#detailCard').classList.add('hidden');
      $('#evalResult').textContent = 'Belum ada evaluasi.';
      showError('');
    });
    $('#evalBtn').addEventListener('click', runEval);

    $$('th[data-sort]').forEach(th => {
      th.addEventListener('click', () => {
        const key = th.dataset.sort;
        if (sortKey === key) sortDir *= -1;
        else { sortKey = key; sortDir = 1; }
        renderTable();
      });
    });

    $('#tabDesign').addEventListener('click', () => switchTab('design'));
    $('#tabEval').addEventListener('click', () => switchTab('eval'));

    $('#exportJson').addEventListener('click', () => {
      if (!requireResult()) return;
      download('qpcr-design.json', JSON.stringify(currentResult, null, 2), 'application/json');
    });
    $('#exportCsv').addEventListener('click', () => {
      if (!requireResult()) return;
      download('qpcr-design.csv', E.exportCsv(currentResult), 'text/csv');
    });
    $('#exportBed').addEventListener('click', () => {
      if (!requireResult()) return;
      download('qpcr-design.bed', E.exportBed(currentResult), 'text/plain');
    });
    $('#exportRdml').addEventListener('click', () => {
      if (!requireResult()) return;
      download('qpcr-design.rdml.xml', E.exportRdml(currentResult, {
        description: 'Desain primer qPCR otomatis',
        targetName: currentResult._targetName || 'target'
      }), 'application/xml');
    });
    $('#copyTop').addEventListener('click', async () => {
      if (!requireResult()) return;
      const p = currentResult.pairs[0];
      const text = `FWD\t${p.fwd.sequence}\nREV\t${p.rev.sequence}\nAmplicon\t${p.ampliconLength} bp\nTm\t${p.fwd.tm}/${p.rev.tm} °C`;
      try {
        await navigator.clipboard.writeText(text);
        $('#copyTop').textContent = 'Tersalin!';
        setTimeout(() => { $('#copyTop').textContent = 'Salin primer #1'; }, 1200);
      } catch (_) {
        showError('Clipboard tidak tersedia; silakan salin manual dari tabel.');
      }
    });
  }

  document.addEventListener('DOMContentLoaded', bind);
})();
