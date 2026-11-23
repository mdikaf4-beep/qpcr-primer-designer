# Deploy qPCR Primer Designer by DF ke Website Publik

Aplikasi ini **100% static** (HTML + CSS + JavaScript), jadi tidak perlu server backend, database, atau runtime Python/Node. Cukup di-host sebagai static site.

Alamat akhir nanti akan seperti:

```
https://<username>.github.io/<nama-repo>/
https://<nama-app>.netlify.app/
https://<nama-app>.pages.dev/
```

Karena semua perhitungan terjadi di browser, **data sekuens pengguna tidak dikirim ke server**.

---

## Opsi 1 — GitHub Pages (Paling Mudah & Gratis)

Saya sudah menyiapkan file deploy di:

```
.github/workflows/deploy-qpcr.yml
```

### Langkah

1. Buat repository di GitHub, misalnya `qpcr-primer-designer`.
2. Di folder workspace ini, jalankan:

   ```bash
   git init
   git add .
   git commit -m "Deploy qPCR Primer Designer by DF"
   git branch -M main
   git remote add origin https://github.com/<USERNAME>/<REPO>.git
   git push -u origin main
   ```

3. Di GitHub, buka repository → **Settings → Pages**.
4. Bagian **Build and deployment → Source**, pilih **GitHub Actions**.
5. Tunggu workflow `Deploy qPCR Primer Designer by DF` selesai (tab **Actions**).
6. Buka URL yang muncul, misalnya:

   ```
   https://<USERNAME>.github.io/<REPO>/
   ```

### Update aplikasi

Setiap kali Anda edit file lalu:

```bash
git add .
git commit -m "Update primer engine"
git push
```

GitHub otomatis deploy ulang.

> File `deploy-qpcr.yml` sudah mengarah ke `./qpcr-primer-designer`, jadi seluruh folder aplikasi akan jadi root website.

---

## Opsi 2 — Netlify (Paling Cepat, Tanpa Git)

Saya sudah menyiapkan `netlify.toml` di root workspace.

### Cara drag & drop

1. Buka https://app.netlify.com/drop
2. Tarik folder **`qpcr-primer-designer`** ke halaman tersebut.
3. Netlify langsung memberi URL seperti `https://random-name-123.netlify.app`.
4. Ganti nama di **Site settings → Change site name**.

### Cara lewat GitHub

1. Import repository yang sama seperti Opsi 1.
2. Netlify akan membaca `netlify.toml`.
3. Klik **Deploy site**.

---

## Opsi 3 — Cloudflare Pages

1. Buka https://dash.cloudflare.com → **Workers & Pages → Create → Pages → Connect to Git**.
2. Pilih repository GitHub Anda.
3. Set:
   - **Build command:** kosongkan
   - **Build output directory:** `qpcr-primer-designer`
4. Deploy.

---

## Opsi 4 — Vercel

1. Buka https://vercel.com/new → Import repository.
2. Set **Root Directory** = `qpcr-primer-designer`.
3. Build command: kosongkan.
4. Output directory: kosongkan atau `.`.
5. Deploy.

---

## Custom Domain + HTTPS

Semua opsi di atas menyediakan **HTTPS gratis** secara otomatis.

Jika punya domain sendiri, misalnya `qpcr.df.id`:

1. Di provider hosting, tambahkan **Custom domain**.
2. Di pengaturan DNS domain, buat record:

   | Type | Name | Value |
   |---|---|---|
   | CNAME | `qpcr` | `username.github.io` / `xxx.netlify.app` / `xxx.pages.dev` |

3. Centang **Force HTTPS**.

Untuk GitHub Pages dengan domain kustom, tambahkan file bernama **`CNAME`** di dalam folder `qpcr-primer-designer` berisi satu baris:

```
qpcr.df.id
```

---

## Hal yang sudah aman

- ✅ Semua path relatif (`engine.js`, `app.js`, `styles.css`) → aman di subdirectory GitHub Pages.
- ✅ Tidak ada API key / token di kode.
- ✅ Tidak ada backend → tidak ada biaya server.
- ✅ Data pengguna diproses lokal di browser.
- ✅ `tests.html` ikut terdeploy, jadi bisa dibuka siapa saja untuk verifikasi.

---

## Jika nanti ingin fitur server

Kalau Anda mau menambahkan:

- simpan hasil desain ke akun user,
- login,
- riwayat proyek,
- integrasi BLAST nyata di server,
- audit trail & e-signature,

maka arsitekturnya perlu ditambah backend. Spesifikasi untuk itu sudah ada di:

```
spesifikasi-software-desain-primer-qpcr.md
```

Tapi untuk “bisa dibuka semua orang dari website”, **static hosting di atas sudah cukup**.
