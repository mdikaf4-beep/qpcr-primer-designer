# qPCR Primer Designer by DF

Perangkat lunak **nyata** untuk merancang dan mengevaluasi primer real-time PCR (qPCR). Berjalan langsung di browser dengan JavaScript murni — **tidak perlu Python, Node, npm, atau instalasi apa pun**.

## Cara menjalankan

1. Buka folder `qpcr-primer-designer`.
2. Klik dua kali **`jalankan.cmd`** (atau langsung **`index.html`** dengan Chrome/Edge/Firefox).
3. Klik **“Jalankan desain”** — contoh sudah termuat otomatis. Tombol **“Muat contoh”** tersedia jika Anda menghapusnya.

Untuk menjalankan pengujian, buka **`tests.html`** di browser yang sama.

Ingin aplikasi ini bisa dibuka semua orang lewat website? Lihat **[DEPLOYMENT.md](DEPLOYMENT.md)** — sudah disiapkan workflow GitHub Pages, Netlify, Cloudflare Pages, dan Vercel.

File yang bisa langsung dipakai:
- `index.html` — aplikasi utama (UI + results)
- `engine.js` — inti algoritma (bisa dipakai ulang dari JS/Node)
- `app.js` — logika UI
- `styles.css` — tampilan
- `tests.html` — 16 test otomatis
- `examples/synthetic_target_600bp.fa` — sekuens contoh

## Fitur

- Input FASTA / sekuens mentah, upload file .fa/.fasta/.txt
- Mode **SYBR** dan **TaqMan** (rekomendasi probe 20–30 nt, Tm 65–70 °C)
- Tm **nearest-neighbor** (SantaLucia 1998) + koreksi garam
- Screening: panjang, GC, Tm, ΔTm, homopolimer, repeat, 3' GC clamp, terminal T
- Hairpin & dimer (heuristik complementarity) + cross-dimer & 3'-dimer
- Masking SNP/repeat via posisi & interval; dukungan exon-junction
- Pemeriksaan off-target 3'-seed (masukkan FASTA pseudogen/off-target)
- Scoring berbobot + breakdown penalti + flag QC (PASS/WARN/FAIL)
- Evaluasi primer existing (mode batch bisa ditambahkan)
- Ekspor **JSON, CSV, BED, RDML (XML)**
- Deterministik: konfigurasi + sekuens sama → hasil sama

## Contoh penggunaan lewat konsol / Node

```js
const E = require('./engine.js');

const seq = 'TAGATATACCCTGCGAACCCGAGTCC...'; // atau E.parseFasta(fasta)[0].seq
const res = E.design(seq, {
  assayType: 'SYBR',
  ampliconMin: 70,
  ampliconMax: 200,
  topN: 5
});

console.log(res.pairs[0]);
console.log(E.exportCsv(res));

// Evaluasi primer yang sudah ada
const p = E.evaluatePair(seq, 'AAAA...', 'TTTT...');
```

## Struktur hasil

```json
{
  "meta": { "engineVersion": "1.0.0", "config": {}, "elapsedMs": 12 },
  "summary": { "forwardCandidates": 210, "reverseCandidates": 205, "pairsBeforeScoring": 320 },
  "pairs": [
    {
      "rank": 1,
      "score": 4.5,
      "fwd": { "sequence": "...", "tm": 60.1, "gc": 52, "start": 120, "end": 140 },
      "rev": { "sequence": "...", "tm": 60.4, "gc": 48, "start": 210, "end": 230 },
      "ampliconLength": 111,
      "deltaTm": 0.3,
      "breakdown": [],
      "qc": [{ "code": "PASS", "severity": "INFO" }]
    }
  ]
}
```

## Batasan penting

- **Spesifisitas genome-wide belum otomatis** jika Anda tidak mengisi kolom “Sekuens off-target”. Untuk diagnostik, pipeline produksi harus menghubungkan engine ini ke **BLAST/Primer-BLAST** terhadap genome referensi (lihat `spesifikasi-software-desain-primer-qpcr.md` bab Specificity Engine).
- Hairpin/dimer memakai **heuristik**, bukan nearest-neighbor termodinamik penuh. Untuk validasi regulated, ganti dengan **ViennaRNA** atau **primer3-py**.
- Efisiensi amplifikasi dan LOQ tetap harus diverifikasi secara eksperimen (MIQE).

## Status

Ini adalah **versi 1.0 fungsional** yang sesuai dengan spesifikasi yang saya buat sebelumnya. Jika Anda mau, langkah berikutnya adalah menambahkan:
1. server API + penyimpanan hasil,
2. integrasi BLAST nyata,
3. modul RDML lengkap dengan thermocycling profile,
4. mode batch 100 primer existing.

## Referensi

- SantaLucia J Jr. *PNAS* 1998 — nearest-neighbor parameters
- Untergasser A et al. *NAR* 2012 — Primer3
- Ye J et al. *NAR* 2012 — Primer-BLAST
- Bustin SA et al. *Clinical Chemistry* 2009 — MIQE
- RDML Consortium — RDML 1.2
