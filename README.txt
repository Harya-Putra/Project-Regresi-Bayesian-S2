PROJECT: Regresi Bayesian + Audit Reproducibility & Data Quality
(Studi Kasus: Yacht Hydrodynamics)

Penulis:
Yaafi Ferdian Syahputra (140720260018)
I Gusti Bagus Harya Putra (140720260016)


ISI FOLDER
────────────────────────────────────────────────────────────────

- Laporan_Project_Regresi_Bayesian_Yacht_final8_fixed_download.Rmd
  Laporan utama dalam format R Markdown yang memuat:
  cover, BAB 1–5, daftar pustaka, acknowledgement, dan lampiran
  (termasuk Lampiran C: Dashboard Analisis dan Lampiran F: Source
  Code Dashboard).

- Laporan_Project_Regresi_Bayesian_Yacht_final8_fixed_download.pdf
  Hasil proses knit laporan dalam format PDF.

- Laporan_Project_Regresi_Bayesian_Yacht_final8_fixed_download.tex
  File LaTeX perantara yang dihasilkan saat knit ke PDF.

- Laporan_Project_Regresi_Bayesian_Yacht.html
  Hasil proses knit dalam format HTML.
  Untuk menampilkan rumus MathJax dan logo dengan sempurna,
  buka file ini menggunakan koneksi internet.

- logo-unpad1.png
  Logo yang digunakan pada cover laporan.

- analisis_yacht_final.R
  Seluruh pipeline komputasi analisis, termasuk:
  audit data, preprocessing, pemodelan Bayesian, validasi,
  diagnostik, posterior predictive check, PSIS-LOO,
  analisis sensitivitas, dan audit reproducibility.
  Bagian 17 mengekspor data untuk dashboard
  (dashboard/dashboard_data.rds).
  File ini dipanggil secara otomatis oleh laporan R Markdown.

- yacht_hydrodynamics.data
  Data mentah Yacht Hydrodynamics dari UCI Machine Learning
  Repository (UCI ID 243).

- dashboard/app.R
  Dashboard R Shiny interaktif (satu file).

- dashboard/dashboard_data.rds
  Hasil analisis yang dibaca oleh dashboard. Angkanya sama dengan
  angka pada laporan.

- fits/
  Objek model brms yang sudah di-fit (M0–M6, sensitivitas prior,
  dan contoh prior salah skala).

- cache/
  Hasil komputasi berat yang disimpan (PSIS-LOO, cross-validation,
  fit multi-seed, prior predictive, Gibbs, dan HMC).

- cv_results.csv
  Ringkasan hasil cross-validation.

- model_comparison_loo.csv
  Ringkasan perbandingan model berdasarkan PSIS-LOO.

- reproducibility_manifest_R.json
  Manifest reproducibility yang mencatat informasi penting
  terkait lingkungan dan proses komputasi.


CARA MENJALANKAN PROJECT
────────────────────────────────────────────────────────────────

1. Buka folder project ini di RStudio dan jadikan folder tersebut
   sebagai working directory.

2. Pastikan seluruh file dalam folder project tetap berada pada
   struktur folder yang sama.

3. Pada Windows, pasang Rtools yang sesuai dengan versi R
   (R 4.5 -> Rtools45, R 4.4 -> Rtools44), lalu restart RStudio.

4. Install package R yang diperlukan dengan menjalankan:

   install.packages(c(
     "brms",
     "rstan",
     "loo",
     "posterior",
     "dplyr",
     "tidyr",
     "kableExtra",
     "patchwork",
     "bayesplot",
     "digest",
     "jsonlite",
     "rmarkdown"
   ))

5. Buka file:

   Laporan_Project_Regresi_Bayesian_Yacht_final8_fixed_download.Rmd

6. Klik tombol "Knit" pada RStudio untuk menghasilkan laporan.


DASHBOARD R SHINY
────────────────────────────────────────────────────────────────

Dashboard hanya membaca hasil yang sudah dihitung
(dashboard/dashboard_data.rds). Karena itu dashboard TIDAK
memerlukan brms, Stan, maupun Rtools, dan dapat langsung dibuka
dalam beberapa detik.

Cara menjalankan:

1. Install package dashboard (cukup sekali):

   install.packages(c("shiny", "bslib", "plotly"))

2. Buka file dashboard/app.R di RStudio, lalu klik "Run App".
   Alternatif: jalankan di Console (working directory = folder
   project ini):

   shiny::runApp("dashboard")

3. File app.R dan dashboard_data.rds harus tetap berada di folder
   yang sama (folder dashboard/).

4. Deploy online (opsional):
   install.packages("rsconnect"), buat akun di shinyapps.io,
   hubungkan token akun, lalu jalankan:

   rsconnect::deployApp("dashboard")

Fitur umum:

- Seluruh teks dashboard ditulis dalam bahasa Inggris.
- Semua grafik interaktif (plotly): arahkan kursor untuk melihat
  nilai, drag untuk zoom, klik legenda untuk menyembunyikan atau
  menampilkan seri, tombol kamera untuk menyimpan grafik (PNG),
  dan tombol di pojok kartu untuk tampilan layar penuh.
- Terdiri atas tiga tab utama. Setiap tab memiliki banner
  pertanyaan utama, empat kotak ringkasan (KPI), dan beberapa
  sub-tab.

Tab 1 – Data-Quality Summary
(Apakah data layak untuk regresi Bayesian?)

- KPI: 308 observasi (22 hull x 14 Froude number), missing value
  dan duplikat 0/0, 14 baris per hull, skewness target 1.75 -> 0.47
  setelah transformasi log(1 + y).
- Resistance curves  : kurva resistansi 22 hull; pilihan skala
                       linear/log, pewarnaan kurva menurut variabel
                       geometri, dan sorot satu hull.
- Hull design space  : sebaran 22 desain hull untuk dua variabel
                       geometri pilihan beserta tabel geometri.
- Target distribution: histogram resistansi mentah vs
                       log(1 + resistansi); jumlah bin dapat diatur.
- Audit tables       : ringkasan audit, jumlah outlier IQR per
                       variabel (beserta penjelasan artefak desain),
                       dan rentang setiap variabel.

Tab 2 – Preprocessing Impact
(Seberapa besar keputusan preprocessing mengubah hasil?)

- KPI: RMSE leave-hull-out M0 -> M3, ELPD-LOO, persentase replikasi
  negatif pada M0, dan model terbaik untuk hull baru.
- Model comparison            : metrik CV leave-hull-out (RMSE, MAE,
                                ELPD, coverage, lebar PI) untuk
                                model pilihan, perbandingan PSIS-LOO
                                terkoreksi Jacobian, tabel CV, dan
                                perbandingan berpasangan (SE per hull).
- Predictions vs actual       : prediksi out-of-sample vs nilai
                                aktual beserta PI 95% per model dan
                                skema validasi (grouped / row-random).
- Posterior predictive checks : densitas data vs replikasi,
                                histogram statistik diskrepansi
                                (mean, sd, max, skewness, proporsi
                                < 0), dan tabel M0 vs M6.
- Prior predictive check      : simulasi resistansi dari prior M0,
                                M6 (prior revisi), dan M6 (prior awal
                                yang terlalu lebar), serta kasus
                                koreksi prior salah skala.
- Posterior inference         : efek geometri M6 (median dan CrI
                                95%), kurva heteroskedastisitas, VIF,
                                korelasi posterior koefisien (M2),
                                dan tabel ringkasan efek.
- New-hull scenario           : SIMULATOR INTERAKTIF. Pilih model
                                (M0/M2/M3/M6), atur Froude number,
                                ambang t, dan lima rasio geometri
                                dengan slider. Keluaran: median
                                prediksi, PI 95%, P(R > t | data),
                                kurva prediktif dengan pita
                                ketidakpastian, serta peringatan
                                ekstrapolasi bila desain jauh dari
                                hull yang pernah diuji. Tombol reset
                                mengembalikan desain median.

Tab 3 – Reproducibility Checklist
(Apakah hasil stabil dan dapat direproduksi?)

- KPI: jumlah item checklist yang lulus, validasi komputasi
  (dalam satuan MCSE), rentang P(R > 20 | Fr = 0.40) antar lima
  seed, dan optimisme CV row-random.
- Checklist                    : status PASS/FAIL setiap item dan
                                 metadata run (seed, SHA-256 data,
                                 versi R/brms/rstan/loo, setelan
                                 sampler).
- Computational validation     : OLS vs posterior analitik vs Gibbs
                                 vs HMC, selisih dalam satuan MCSE,
                                 pengecekan log-likelihood manual vs
                                 Stan, dan diagnostik MCMC per model
                                 (R-hat, ESS, MCSE, waktu sampling).
- Seed and row-order stability : selisih posterior terhadap run
                                 acuan (seed berbeda / urutan data
                                 diacak), audit lima seed (SD
                                 antar-run / MCSE), dan kesimpulan
                                 antar seed.
- Validation scheme            : CV row-random vs leave-hull-out
                                 untuk M4 dan M6 beserta tabel
                                 optimisme.
- Sensitivity                  : sensitivitas prior M6 (x0.5, x1,
                                 x2.5), sensitivitas likelihood
                                 (Student-t), dan simulasi kompetisi
                                 prior vs data untuk 1, 3, dan 22 hull.

Prediksi pada simulator dihitung langsung dari 1.000 draw posterior
dengan rumus yang sama dengan model brms (sudah diverifikasi),
sehingga dashboard tetap ringan.


PROSES KNIT DAN CACHE
────────────────────────────────────────────────────────────────

Pada knit pertama, seluruh pipeline analisis dapat dijalankan dari
awal. Proses ini mencakup kompilasi Stan, fitting model Bayesian,
MCMC, cross-validation, dan berbagai analisis diagnostik.

Waktu komputasi knit pertama diperkirakan sekitar 30–60 menit,
tergantung pada spesifikasi CPU, RAM, versi R, dan versi package
yang digunakan.

Hasil komputasi yang membutuhkan waktu lama disimpan dalam folder:

- fits/
- cache/

Dengan demikian, proses knit berikutnya dapat berjalan jauh lebih
cepat karena hasil komputasi yang sudah tersedia akan digunakan
kembali.

Setiap knit juga memperbarui dashboard/dashboard_data.rds secara
otomatis, sehingga angka pada dashboard selalu sama dengan angka
pada laporan.

Untuk menjalankan seluruh analisis dari awal, hapus folder:

- fits/
- cache/

kemudian lakukan knit ulang.


REPRODUCIBILITY
────────────────────────────────────────────────────────────────

Pipeline dirancang agar analisis dapat direproduksi secara otomatis
tanpa langkah komputasi manual. Seed, preprocessing, model,
validasi, dan hasil analisis dikendalikan melalui kode dalam
analisis_yacht_final.R.

Hasil pada komputer lain dapat berbeda pada beberapa digit terakhir
karena Monte Carlo error, versi R, versi package, atau konfigurasi
komputasi yang berbeda. Perbedaan tersebut merupakan variasi
numerik yang wajar dalam komputasi Bayesian dan tidak diharapkan
mengubah kesimpulan substantif analisis.


CATATAN PENTING
────────────────────────────────────────────────────────────────

- Jangan memindahkan atau mengganti nama file utama karena
  analisis_yacht_final.R, file R Markdown, dan dashboard
  menggunakan struktur folder project ini.

- Untuk reproduksi hasil secara konsisten, gunakan seluruh isi folder
  project, termasuk data mentah, file CSV, manifest JSON, serta
  folder fits/, cache/, dan dashboard/.

- Seluruh angka yang ditampilkan di dalam laporan dihasilkan secara
  otomatis melalui kode R dan inline R expression, bukan diketik
  secara manual.

- Jika folder fits/ dan cache/ tersedia, knit ulang akan menggunakan
  hasil yang telah disimpan sehingga waktu komputasi menjadi lebih
  singkat.

- Jika dashboard menampilkan pesan error, app.R akan menyebutkan
  package mana yang belum terpasang atau gagal dimuat beserta cara
  memperbaikinya. Pesan "dashboard_data.rds was not found" berarti
  app.R dan dashboard_data.rds tidak berada di folder yang sama.
