/*!
 * qPCR Primer Designer by DF — core engine (vanilla JS, no dependencies)
 * Version 1.0.0
 *
 * Implementasi: nearest-neighbor Tm (SantaLucia 1998), screening kandidat,
 * assembly pasangan, scoring, QC, dan eksport RDML/CSV/BED.
 *
 * Catatan: hairpin/dimer memakai heuristik complementarity. Untuk produksi
 * diagnostik, ganti dengan ViennaRNA / primer3-py dan tambahkan BLAST.
 */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  if (root) root.PrimerEngine = api;
})(typeof window !== 'undefined' ? window : globalThis, function () {
  'use strict';

  const VERSION = '1.0.0';

  const NN = {
    AA: { h: -7.9, s: -22.2 },
    AT: { h: -7.2, s: -20.4 },
    TA: { h: -7.2, s: -21.3 },
    CA: { h: -8.5, s: -22.7 },
    GT: { h: -8.4, s: -22.4 },
    CT: { h: -7.8, s: -21.0 },
    GA: { h: -8.2, s: -22.2 },
    CG: { h: -10.6, s: -27.2 },
    GC: { h: -9.8, s: -24.4 },
    GG: { h: -8.0, s: -19.9 }
  };

  const NN_KEY = {
    AA: 'AA', TT: 'AA',
    AT: 'AT', TA: 'TA',
    CA: 'CA', TG: 'CA',
    GT: 'GT', AC: 'GT',
    CT: 'CT', AG: 'CT',
    GA: 'GA', TC: 'GA',
    CG: 'CG', GC: 'GC',
    GG: 'GG', CC: 'GG'
  };

  const DEFAULT_CONFIG = {
    minLen: 18,
    maxLen: 25,
    hardGcMin: 30,
    hardGcMax: 70,
    gcMin: 40,
    gcMax: 60,
    hardTmMin: 50,
    hardTmMax: 70,
    tmIdeal: 60,
    tmSoftMin: 58,
    tmSoftMax: 62,
    maxDeltaTm: 2,
    maxDeltaTmHard: 5,
    ampliconMin: 70,
    ampliconMax: 200,
    ampliconSoftMax: 150,
    maxHomopolymer: 4,
    maxRepeat: 3,
    clampMin: 1,
    clampMax: 2,
    na: 0.05,
    strandConc: 5e-7,
    minDimerDgHard: -12,
    minCrossDimerDgHard: -12,
    minAmpliconDgHard: -20,
    requireJunction: false,
    assayType: 'SYBR',
    topN: 10,
    maxPairsToScore: 8000,
    mask: [],
    snps: [],
    junctions: [],
    offtarget: []
  };

  function cleanSeq(input) {
    return String(input || '')
      .toUpperCase()
      .replace(/U/g, 'T')
      .replace(/[^ACGTN]/g, '');
  }

  function revcomp(seq) {
    const map = { A: 'T', T: 'A', G: 'C', C: 'G', N: 'N' };
    return cleanSeq(seq).split('').reverse().map(b => map[b] || 'N').join('');
  }

  function parseFasta(text) {
    const records = [];
    let name = null;
    let chunks = [];
    for (const raw of String(text || '').split(/\r?\n/)) {
      const line = raw.trim();
      if (!line) continue;
      if (line.startsWith('>')) {
        if (name !== null) records.push({ name, seq: cleanSeq(chunks.join('')) });
        name = line.slice(1).trim().split(/\s+/)[0] || `seq_${records.length + 1}`;
        chunks = [];
      } else if (name !== null) {
        chunks.push(line);
      }
    }
    if (name !== null) records.push({ name, seq: cleanSeq(chunks.join('')) });
    if (!records.length && String(text || '').trim()) {
      records.push({ name: 'input', seq: cleanSeq(text) });
    }
    return records.filter(r => r.seq.length > 0);
  }

  function parseIntervals(text) {
    if (!text || !String(text).trim()) return [];
    return String(text).split(/[,;\s]+/).map(part => {
      const m = part.match(/^(\d+)\s*[-:]\s*(\d+)$/);
      if (!m) return null;
      const a = parseInt(m[1], 10) - 1;
      const b = parseInt(m[2], 10) - 1;
      return [Math.min(a, b), Math.max(a, b)];
    }).filter(Boolean);
  }

  function parsePositions(text) {
    if (!text || !String(text).trim()) return [];
    return Array.from(new Set(String(text).split(/[,;\s]+/).map(x => parseInt(x, 10) - 1).filter(n => !isNaN(n) && n >= 0))).sort((a, b) => a - b);
  }

  function gcPercent(seq) {
    const s = cleanSeq(seq);
    if (!s.length) return 0;
    const gc = (s.match(/[GC]/g) || []).length;
    return +(100 * gc / s.length).toFixed(2);
  }

  function homopolymerMax(seq) {
    let best = 0, run = 0, prev = '';
    for (const b of cleanSeq(seq)) {
      run = (b === prev) ? run + 1 : 1;
      if (run > best) best = run;
      prev = b;
    }
    return best;
  }

  function maxRepeat(seq) {
    const s = cleanSeq(seq);
    let best = 1;
    // k=1 adalah homopolymer dan sudah ditangani homopolymerMax().
    for (let k = 2; k <= 3; k++) {
      for (let i = 0; i + 2 * k <= s.length; i++) {
        const a = s.slice(i, i + k), b = s.slice(i + k, i + 2 * k);
        if (a === b) {
          let n = 2, j = i + 2 * k;
          while (j + k <= s.length && s.slice(j, j + k) === a) { n++; j += k; }
          if (n > best) best = n;
        }
      }
    }
    return best;
  }

  function gcClampCount(seq, window) {
    const tail = cleanSeq(seq).slice(-(window || 5));
    return (tail.match(/[GC]/g) || []).length;
  }

  function nearestNeighborTm(seq, opts) {
    const o = Object.assign({}, DEFAULT_CONFIG, opts || {});
    const s = cleanSeq(seq);
    if (s.length < 2) return 0;
    if (/N/.test(s)) return NaN;

    let h = 0, sumS = 0;
    for (let i = 0; i < s.length - 1; i++) {
      const key = NN_KEY[s.slice(i, i + 2)];
      if (!key) return NaN;
      h += NN[key].h;
      sumS += NN[key].s;
    }

    const initH = b => (b === 'A' || b === 'T') ? 2.2 : 0.2;
    const initS = b => (b === 'A' || b === 'T') ? 6.9 : -2.8;
    h += initH(s[0]) + initH(s[s.length - 1]);
    sumS += initS(s[0]) + initS(s[s.length - 1]);

    sumS += 0.368 * (s.length - 1) * Math.log(Math.max(o.na, 1e-6));
    const R = 1.987;
    const ct = Math.max(o.strandConc, 1e-12);
    const tm = (h * 1000) / (sumS + R * Math.log(ct / 4)) - 273.15;
    return +tm.toFixed(2);
  }

  function pairScore(a, b) {
    if ((a === 'G' || a === 'C') && (b === 'G' || b === 'C')) return 1.5;
    if ((a === 'A' || a === 'T') && (b === 'A' || b === 'T')) return 1.0;
    return 0;
  }

  /** Skor complementarity lokal antara a dan anti-paralel b. */
  function maxComplementarity(a, b) {
    const s1 = cleanSeq(a);
    const s2 = revcomp(b);
    if (!s1.length || !s2.length) return 0;
    let prev = new Array(s2.length + 1).fill(0);
    let best = 0;
    for (let i = 1; i <= s1.length; i++) {
      const cur = new Array(s2.length + 1).fill(0);
      for (let j = 1; j <= s2.length; j++) {
        const sc = pairScore(s1[i - 1], s2[j - 1]);
        if (sc > 0) {
          cur[j] = prev[j - 1] + sc;
          if (cur[j] > best) best = cur[j];
        }
      }
      prev = cur;
    }
    return +best.toFixed(2);
  }

  function dgFromComplementarity(score) {
    return +(-1.4 * score).toFixed(2);
  }

  function hairpinDg(seq, minLoop) {
    const s = cleanSeq(seq);
    const loop = Math.max(3, minLoop || 3);
    let best = 0;
    for (let i = loop; i <= s.length - loop; i++) {
      const sc = maxComplementarity(s.slice(0, i), s.slice(i));
      const dg = dgFromComplementarity(sc);
      if (dg < best) best = dg;
    }
    return best;
  }

  /** Hairpin cepat untuk primer pendek: pasangkan ujung 5' dan 3'. */
  function primerHairpinDg(seq) {
    const s = cleanSeq(seq);
    const f = Math.max(4, Math.floor(s.length / 2) - 3);
    return dgFromComplementarity(maxComplementarity(s.slice(0, f), s.slice(-f)));
  }

  /** Approximasi struktur sekunder amplicon dengan memasangkan ujung 5' dan 3'. */
  function ampliconHairpinDg(amp, flank) {
    const s = cleanSeq(amp);
    const f = Math.max(15, Math.min(flank || 40, Math.floor(s.length / 2) - 5));
    if (f < 15 || s.length < 40) return 0;
    return dgFromComplementarity(maxComplementarity(s.slice(0, f), s.slice(-f)));
  }

  function dimerDg(a, b, threePrimeOnly) {
    const s1 = cleanSeq(a);
    const s2 = cleanSeq(b);
    if (threePrimeOnly) {
      return dgFromComplementarity(maxComplementarity(s1.slice(-8), s2.slice(-8)));
    }
    return dgFromComplementarity(maxComplementarity(s1, s2));
  }

  function overlaps(start, end, intervals) {
    return (intervals || []).some(([a, b]) => start <= b && end >= a);
  }

  function snpNearThreePrime(start, end, orientation, snps) {
    const win = orientation === 'FWD' ? [end - 4, end] : [start, start + 4];
    return (snps || []).some(p => p >= win[0] && p <= win[1]);
  }

  function makeCandidate(sequence, start, orientation, cfg) {
    return {
      orientation,
      sequence,
      start,
      end: start + sequence.length - 1,
      length: sequence.length,
      gc: gcPercent(sequence),
      tm: nearestNeighborTm(sequence, cfg),
      homopolymer: homopolymerMax(sequence),
      repeat: maxRepeat(sequence),
      clamp: gcClampCount(sequence, 5),
      terminalBase: sequence[sequence.length - 1],
      hairpinDg: primerHairpinDg(sequence),
      selfDimerDg: dimerDg(sequence, sequence)
    };
  }

  function generateCandidates(seq, cfg) {
    const c = Object.assign({}, DEFAULT_CONFIG, cfg || {});
    const s = cleanSeq(seq);
    const out = [];
    for (let len = c.minLen; len <= c.maxLen; len++) {
      for (let i = 0; i + len <= s.length; i++) {
        const raw = s.slice(i, i + len);
        if (/N/.test(raw)) continue;
        if (overlaps(i, i + len - 1, c.mask)) continue;
        if (homopolymerMax(raw) > c.maxHomopolymer) continue;
        if (maxRepeat(raw) > c.maxRepeat) continue;
        const gc = gcPercent(raw);
        if (gc < c.hardGcMin || gc > c.hardGcMax) continue;
        const tm = nearestNeighborTm(raw, c);
        if (!isFinite(tm) || tm < c.hardTmMin || tm > c.hardTmMax) continue;

        if (!snpNearThreePrime(i, i + len - 1, 'FWD', c.snps)) {
          out.push(makeCandidate(raw, i, 'FWD', c));
        }
        if (!snpNearThreePrime(i, i + len - 1, 'REV', c.snps)) {
          out.push(makeCandidate(revcomp(raw), i, 'REV', c));
        }
      }
    }
    return out;
  }

  function spansJunction(pair, junctions) {
    if (!junctions || !junctions.length) return true;
    return junctions.some(j =>
      (j >= pair.fwd.start && j <= pair.fwd.end) ||
      (j >= pair.rev.start && j <= pair.rev.end) ||
      (j > pair.fwd.end && j < pair.rev.start)
    );
  }

  function assemblePairs(seq, fwd, rev, cfg) {
    const c = Object.assign({}, DEFAULT_CONFIG, cfg || {});
    const s = cleanSeq(seq);
    const cheap = [];

    // Fase 1: indeks reverse berdasarkan panjang + posisi, lalu binary search
    // rentang amplicon. Menghindari O(FWD × REV) yang sangat mahal.
    const byLen = new Map();
    for (const r of rev) {
      if (!byLen.has(r.length)) byLen.set(r.length, []);
      byLen.get(r.length).push(r);
    }
    byLen.forEach(arr => arr.sort((a, b) => a.start - b.start));

    function lowerBound(arr, value) {
      let lo = 0, hi = arr.length;
      while (lo < hi) {
        const mid = (lo + hi) >> 1;
        if (arr[mid].start < value) lo = mid + 1; else hi = mid;
      }
      return lo;
    }

    for (const f of fwd) {
      for (const [revLen, arr] of byLen) {
        const lo = Math.max(f.end + 1, f.start + c.ampliconMin - revLen);
        const hi = f.start + c.ampliconMax - revLen;
        if (lo > hi) continue;

        let idx = lowerBound(arr, lo);
        while (idx < arr.length && arr[idx].start <= hi) {
          const r = arr[idx++];
          const ampLen = r.end - f.start + 1;
          if (ampLen < c.ampliconMin || ampLen > c.ampliconMax) continue;
          const deltaTm = +Math.abs(f.tm - r.tm).toFixed(2);
          if (deltaTm > c.maxDeltaTmHard) continue;

          const pair = {
            fwd: f,
            rev: r,
            amplicon: s.slice(f.start, r.end + 1),
            ampliconLength: ampLen,
            ampliconGc: 0,
            ampliconDg: 0,
            deltaTm,
            crossDimer: 0,
            threePrimeDimer: 0
          };
          if (c.requireJunction && !spansJunction(pair, c.junctions)) continue;
          pair.spansJunction = spansJunction(pair, c.junctions);
          pair._heuristic =
            Math.abs(f.tm - c.tmIdeal) + Math.abs(r.tm - c.tmIdeal) + deltaTm +
            Math.abs(f.gc - 50) / 10 + Math.abs(r.gc - 50) / 10;
          cheap.push(pair);
        }
      }
    }

    // Fase 2: batasi jumlah pasangan sebelum evaluasi struktur yang lebih mahal.
    cheap.sort((a, b) => a._heuristic - b._heuristic);
    const shortlist = cheap.slice(0, c.maxPairsToScore);

    // Fase 3: evaluasi mahal hanya untuk shortlist.
    const pairs = [];
    for (const p of shortlist) {
      p.crossDimer = dimerDg(p.fwd.sequence, p.rev.sequence);
      if (p.crossDimer < c.minCrossDimerDgHard) continue;
      p.threePrimeDimer = dimerDg(p.fwd.sequence, p.rev.sequence, true);
      p.ampliconGc = gcPercent(p.amplicon);
      p.ampliconDg = ampliconHairpinDg(p.amplicon, 40);
      if (p.ampliconDg < c.minAmpliconDgHard) continue;
      delete p._heuristic;
      pairs.push(p);
    }
    return pairs;
  }

  function penalty(list, code, value, detail) {
    if (value > 0) list.push({ code, penalty: +value.toFixed(2), detail: detail || '' });
  }

  function scorePair(pair, cfg) {
    const c = Object.assign({}, DEFAULT_CONFIG, cfg || {});
    const b = [];
    const qc = [];

    for (const [label, p] of [['FWD', pair.fwd], ['REV', pair.rev]]) {
      penalty(b, `TM_DEV_${label}`, Math.max(0, Math.abs(p.tm - c.tmIdeal) - (c.tmSoftMax - c.tmIdeal)) * 2, `${p.tm} °C`);
      penalty(b, `GC_${label}`, Math.max(0, c.gcMin - p.gc) + Math.max(0, p.gc - c.gcMax), `${p.gc}%`);
      penalty(b, `LEN_${label}`, Math.max(0, 20 - p.length) + Math.max(0, p.length - 23), `${p.length} nt`);
      if (p.clamp === 0) penalty(b, `CLAMP_${label}`, 8, '3\' clamp = 0');
      if (p.clamp > c.clampMax) penalty(b, `CLAMP_${label}`, 6, `3\' clamp = ${p.clamp}`);
      if (p.terminalBase === 'T') penalty(b, `TERMINAL_T_${label}`, 4, '3\' terminal T');
      if (p.homopolymer === 4) penalty(b, `HOMOPOLYMER_${label}`, 10, 'run = 4');
      if (p.hairpinDg < -9) penalty(b, `HAIRPIN_${label}`, 12, `${p.hairpinDg} kcal/mol`);
      if (p.selfDimerDg < -9) penalty(b, `SELF_DIMER_${label}`, 15, `${p.selfDimerDg} kcal/mol`);
    }

    penalty(b, 'DELTA_TM', Math.max(0, pair.deltaTm - c.maxDeltaTm) * 3, `${pair.deltaTm} °C`);
    if (pair.crossDimer < -9) penalty(b, 'CROSS_DIMER', 15, `${pair.crossDimer} kcal/mol`);
    if (pair.threePrimeDimer < -9) penalty(b, 'THREE_PRIME_DIMER', 10, `${pair.threePrimeDimer} kcal/mol`);
    if (pair.ampliconLength > c.ampliconSoftMax) penalty(b, 'AMPLICON_LONG', 5, `${pair.ampliconLength} bp`);
    if (pair.ampliconLength < c.ampliconMin) penalty(b, 'AMPLICON_SHORT', 8, `${pair.ampliconLength} bp`);

    if (pair.threePrimeDimer < -11) qc.push({ code: 'DIMER_RISK', severity: 'WARN', message: 'Dimer 3\' kuat terdeteksi' });
    if (!c.offtarget || !c.offtarget.length) {
      qc.push({ code: 'NO_OFFTARGET_CHECK', severity: 'WARN', message: 'Spesifisitas genome-wide belum dijalankan' });
    } else {
      const hits = checkOfftarget(c.offtarget, pair);
      hits.forEach(h => qc.push(h));
    }
    if (!qc.some(q => q.severity === 'FAIL')) qc.unshift({ code: 'PASS', severity: 'INFO', message: 'Lolos screening internal' });

    const total = +b.reduce((sum, x) => sum + x.penalty, 0).toFixed(2);
    return { total, breakdown: b, qc };
  }

  function hamming(a, b) {
    let d = 0;
    for (let i = 0; i < a.length; i++) if (a[i] !== b[i]) d++;
    return d;
  }

  function checkOfftarget(offtargetSeqs, pair) {
    const flags = [];
    const seeds = [pair.fwd.sequence.slice(-8), pair.rev.sequence.slice(-8)];
    for (const rec of offtargetSeqs) {
      const s = typeof rec === 'string' ? cleanSeq(rec) : cleanSeq(rec.seq || '');
      const name = typeof rec === 'string' ? 'offtarget' : (rec.name || 'offtarget');
      for (const seed of seeds) {
        for (let i = 0; i + seed.length <= s.length; i++) {
          if (hamming(seed, s.slice(i, i + seed.length)) <= 1) {
            flags.push({ code: 'OFFTARGET_3PRIME', severity: 'FAIL', message: `3'-seed ${seed} cocok di ${name}:${i + 1}` });
            return flags;
          }
        }
      }
    }
    flags.push({ code: 'OFFTARGET_CLEAR', severity: 'INFO', message: "Tidak ada 3'-seed yang cocok pada sekuens off-target" });
    return flags;
  }

  function design(seq, cfg) {
    const c = Object.assign({}, DEFAULT_CONFIG, cfg || {});
    const started = Date.now();
    const fwdAll = generateCandidates(seq, c);
    const revAll = fwdAll.filter(x => x.orientation === 'REV');
    const fwd = fwdAll.filter(x => x.orientation === 'FWD');
    const prePairs = assemblePairs(seq, fwd, revAll, c);

    const scored = prePairs.map(p => {
      const sc = scorePair(p, c);
      return Object.assign({}, p, { score: sc.total, breakdown: sc.breakdown, qc: sc.qc });
    }).sort((a, b) => a.score - b.score || a.ampliconLength - b.ampliconLength);

    scored.forEach((p, i) => { p.rank = i + 1; });

    return {
      meta: {
        engineVersion: VERSION,
        assayType: c.assayType,
        config: c,
        sequenceLength: cleanSeq(seq).length,
        elapsedMs: Date.now() - started
      },
      summary: {
        forwardCandidates: fwd.length,
        reverseCandidates: revAll.length,
        pairsBeforeScoring: prePairs.length,
        pairsReturned: scored.length
      },
      pairs: scored.slice(0, c.topN)
    };
  }

  function evaluatePair(seq, fwdSeq, revSeq, cfg) {
    const c = Object.assign({}, DEFAULT_CONFIG, cfg || {});
    const s = cleanSeq(seq);
    const fRaw = cleanSeq(fwdSeq);
    const rRaw = cleanSeq(revSeq);
    const i = s.indexOf(fRaw);
    const j = s.indexOf(revcomp(rRaw));
    if (i < 0 || j < 0) throw new Error('Primer tidak ditemukan pada sekuens (cek orientasi/reverse primer).');
    const f = makeCandidate(fRaw, i, 'FWD', c);
    const r = makeCandidate(revcomp(s.slice(j, j + rRaw.length)), j, 'REV', c);
    const ampLen = r.end - f.start + 1;
    const pair = {
      fwd: f, rev: r,
      amplicon: s.slice(f.start, r.end + 1),
      ampliconLength: ampLen,
      ampliconGc: gcPercent(s.slice(f.start, r.end + 1)),
      ampliconDg: ampliconHairpinDg(s.slice(f.start, r.end + 1), 40),
      deltaTm: +Math.abs(f.tm - r.tm).toFixed(2),
      crossDimer: dimerDg(f.sequence, r.sequence),
      threePrimeDimer: dimerDg(f.sequence, r.sequence, true),
      spansJunction: spansJunction({ fwd: f, rev: r }, c.junctions)
    };
    const sc = scorePair(pair, c);
    return Object.assign({}, pair, { rank: 1, score: sc.total, breakdown: sc.breakdown, qc: sc.qc });
  }

  function exportCsv(result) {
    const head = ['rank', 'score', 'fwd_seq', 'rev_seq', 'fwd_tm', 'rev_tm', 'delta_tm', 'fwd_gc', 'rev_gc', 'amplicon_bp', 'amplicon_gc', 'qc'];
    const rows = result.pairs.map(p => [
      p.rank, p.score, p.fwd.sequence, p.rev.sequence, p.fwd.tm, p.rev.tm,
      p.deltaTm, p.fwd.gc, p.rev.gc, p.ampliconLength, p.ampliconGc,
      p.qc.map(q => q.code).join('|')
    ]);
    return [head, ...rows].map(r => r.map(v => `"${String(v).replace(/"/g, '""')}"`).join(',')).join('\r\n');
  }

  function exportBed(result) {
    return result.pairs.map(p =>
      `target\t${p.fwd.start}\t${p.rev.end + 1}\t${p.fwd.sequence}/${p.rev.sequence}\t${p.score}\t+`
    ).join('\n');
  }

  function esc(s) {
    return String(s).replace(/[<>&"']/g, ch => ({ '<': '&lt;', '>': '&gt;', '&': '&amp;', '"': '&quot;', "'": '&apos;' }[ch]));
  }

  function exportRdml(result, meta) {
    const m = meta || {};
    const lines = [];
    lines.push('<?xml version="1.0" encoding="UTF-8"?>');
    lines.push('<rdml version="1.2" xmlns="http://www.rdml.org">');
    lines.push('  <date><time>' + new Date().toISOString() + '</time></date>');
    lines.push('  <experiment id="primer-design-' + Date.now() + '">');
    lines.push('    <description>' + esc(m.description || 'Desain primer qPCR otomatis') + '</description>');
    lines.push('    <documentation><miqeRequirements>' + esc(m.targetName || 'target') + '</miqeRequirements></documentation>');
    lines.push('    <primers>');
    result.pairs.forEach((p, i) => {
      lines.push('      <primer id="pair_' + (i + 1) + '">');
      lines.push('        <name>' + esc('pair_' + (i + 1)) + '</name>');
      lines.push('        <sequence>' + esc(p.fwd.sequence) + '</sequence>');
      lines.push('        <type>forward</type>');
      lines.push('        <tm>' + esc(p.fwd.tm) + '</tm>');
      lines.push('        <gc>' + esc(p.fwd.gc) + '</gc>');
      lines.push('      </primer>');
      lines.push('      <primer id="pair_' + (i + 1) + '_rev">');
      lines.push('        <name>' + esc('pair_' + (i + 1) + '_rev') + '</name>');
      lines.push('        <sequence>' + esc(p.rev.sequence) + '</sequence>');
      lines.push('        <type>reverse</type>');
      lines.push('        <tm>' + esc(p.rev.tm) + '</tm>');
      lines.push('        <gc>' + esc(p.rev.gc) + '</gc>');
      lines.push('      </primer>');
    });
    lines.push('    </primers>');
    lines.push('    <amplicons>');
    result.pairs.forEach((p, i) => {
      lines.push('      <amplicon id="amp_' + (i + 1) + '">');
      lines.push('        <sequence>' + esc(p.amplicon) + '</sequence>');
      lines.push('        <length>' + p.ampliconLength + '</length>');
      lines.push('      </amplicon>');
    });
    lines.push('    </amplicons>');
    lines.push('  </experiment>');
    lines.push('</rdml>');
    return lines.join('\n');
  }

  return {
    VERSION,
    DEFAULT_CONFIG,
    cleanSeq, revcomp, parseFasta, parseIntervals, parsePositions,
    gcPercent, homopolymerMax, maxRepeat, gcClampCount,
    nearestNeighborTm, hairpinDg, dimerDg, maxComplementarity,
    generateCandidates, assemblePairs, scorePair, checkOfftarget,
    design, evaluatePair,
    exportCsv, exportBed, exportRdml
  };
});
