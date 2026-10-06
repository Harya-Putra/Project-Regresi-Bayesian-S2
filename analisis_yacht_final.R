# =============================================================================
# analisis_yacht.R
# Regresi Bayesian + audit reproducibility & data quality - Yacht Hydrodynamics
# Seluruh pipeline dapat dijalankan ulang tanpa AI:  source("analisis_yacht.R")
# Hasil fit disimpan di folder fits/ dan cache/ agar knit ulang cepat & identik.
# =============================================================================

suppressPackageStartupMessages({
  library(brms)      # regresi Bayesian (Stan/HMC)
  library(rstan)     # model Stan konjugat untuk validasi komputasi
  library(loo)       # PSIS-LOO dan WAIC
  library(posterior) # ringkasan draws, R-hat, ESS, MCSE
  library(dplyr)
  library(tidyr)
  library(splines)
})

SEED <- 2026
set.seed(SEED)
N_CORES <- min(4, parallel::detectCores())
options(mc.cores = N_CORES)
rstan_options(auto_write = TRUE)
dir.create("fits", showWarnings = FALSE)
dir.create("cache", showWarnings = FALSE)

# Hasil mahal disimpan sebagai .rds; hapus folder cache/ dan fits/ untuk menghitung ulang dari nol
cache_rds <- function(name, expr) {
  f <- file.path("cache", paste0(name, ".rds"))
  if (file.exists(f)) return(readRDS(f))
  val <- expr
  saveRDS(val, f)
  val
}

SAMPLER    <- list(chains = 4, iter = 2000, warmup = 1000)   # 4 x 1000 draws
CV_SAMPLER <- list(chains = 2, iter = 1000, warmup = 500)    # konfigurasi ringan khusus CV
CONTROL    <- list(adapt_delta = 0.95)

# -----------------------------------------------------------------------------
# 1. DATA: load, nama kolom sesuai dokumentasi UCI 243, checksum
# -----------------------------------------------------------------------------
COLS <- c("LCB", "Prismatic_Coef", "LD_Ratio", "BD_Ratio", "LB_Ratio", "Froude_No", "Resistance")
PRED <- COLS[1:6]
GEOM <- COLS[1:5]
RAW_PATH <- "yacht_hydrodynamics.data"
df <- read.table(RAW_PATH, col.names = COLS)
RAW_SHA256 <- digest::digest(file = RAW_PATH, algo = "sha256")

geom_key <- do.call(paste, df[GEOM])
df$hull_id <- as.integer(factor(geom_key, levels = unique(geom_key)))

sdp  <- function(x) sqrt(mean((x - mean(x))^2))                   # SD populasi (ddof = 0)
skew <- function(x) mean((x - mean(x))^3) / sdp(x)^3

# -----------------------------------------------------------------------------
# 2. AUDIT KUALITAS DATA
# -----------------------------------------------------------------------------
rows_per_hull <- table(df$hull_id)
iqr_out <- sapply(COLS, function(v) {
  q <- quantile(df[[v]], c(.25, .75)); w <- diff(q)
  sum(df[[v]] < q[1] - 1.5 * w | df[[v]] > q[2] + 1.5 * w)
})
AUDIT <- list(
  n_rows = nrow(df), n_cols = length(COLS), n_missing = sum(is.na(df[COLS])),
  n_dup = sum(duplicated(df[COLS])), n_hull = length(unique(df$hull_id)),
  rows_per_hull = unique(as.integer(rows_per_hull)), n_froude = length(unique(df$Froude_No)),
  skew_raw = skew(df$Resistance), skew_log = skew(log1p(df$Resistance)),
  shapiro_p = shapiro.test(df$Resistance)$p.value, iqr_out = iqr_out,
  range = t(sapply(COLS, function(v) c(min = min(df[[v]]), max = max(df[[v]]),
                                       mean = mean(df[[v]]), sd = sd(df[[v]]))))
)

# -----------------------------------------------------------------------------
# 3. PREPROCESSING EKSPLISIT (statistik dihitung dari data TRAIN saja)
# -----------------------------------------------------------------------------
make_frames <- function(train, test = NULL) {
  sc <- lapply(PRED, function(v) c(mean(train[[v]]), sdp(train[[v]])))
  names(sc) <- PRED
  sc$y  <- c(mean(train$Resistance), sdp(train$Resistance))
  ly    <- log1p(train$Resistance)
  sc$ly <- c(mean(ly), sdp(ly))
  basis <- bs((train$Froude_No - sc$Froude_No[1]) / sc$Froude_No[2], df = 6)  # B-spline kubik, knot dari train
  tf <- function(d) {
    o <- as.data.frame(lapply(setNames(PRED, paste0("z_", PRED)),
                              function(v) (d[[v]] - sc[[v]][1]) / sc[[v]][2]))
    o$zF2 <- o$z_Froude_No^2
    o$zF3 <- o$z_Froude_No^3
    B <- predict(basis, o$z_Froude_No)
    colnames(B) <- paste0("F", 1:6)
    o <- cbind(o, as.data.frame(B))
    o$y_z  <- (d$Resistance - sc$y[1]) / sc$y[2]
    o$ly_z <- (log1p(d$Resistance) - sc$ly[1]) / sc$ly[2]
    o$hull <- sprintf("H%02d", d$hull_id)
    o$hull_id <- d$hull_id
    o$Froude_No <- d$Froude_No
    o$Resistance <- d$Resistance
    o
  }
  list(train = tf(train), test = if (!is.null(test)) tf(test), sc = sc, tf = tf)
}

full <- df[c(PRED, "Resistance", "hull_id")]
FR <- make_frames(full)
mdl <- FR$train
SCALE <- FR$sc

to_raw <- function(z, resp, sc) if (resp == "raw") z * sc$y[2] + sc$y[1] else expm1(z * sc$ly[2] + sc$ly[1])
jacobian <- function(frame, resp, sc) {
  if (resp == "raw") rep(-log(sc$y[2]), nrow(frame)) else -log(sc$ly[2]) - log1p(frame$Resistance)
}

# -----------------------------------------------------------------------------
# 4. SPESIFIKASI MODEL & PRIOR
# -----------------------------------------------------------------------------
GEO <- paste0("z_", GEOM, collapse = " + ")
LIN <- paste(GEO, "+ z_Froude_No")
SPL <- paste0("F", 1:6, collapse = " + ")

SPECS <- list(
  M0 = list(label = "M0 Baseline: Gaussian, raw target (Variant 1)", f = paste("y_z ~", LIN),
            resp = "raw", fam = "gaussian", froude = "z_Froude_No", spline = FALSE, hetero = FALSE, hier = FALSE),
  M1 = list(label = "M1 Log target, linear (Variant 2)", f = paste("ly_z ~", LIN),
            resp = "log", fam = "gaussian", froude = "z_Froude_No", spline = FALSE, hetero = FALSE, hier = FALSE),
  M2 = list(label = "M2 Log target + Froude^2, ^3 (Variant 3)", f = paste("ly_z ~", LIN, "+ zF2 + zF3"),
            resp = "log", fam = "gaussian", froude = c("z_Froude_No", "zF2", "zF3"), spline = FALSE, hetero = FALSE, hier = FALSE),
  M3 = list(label = "M3 Log target + Froude spline", f = paste("ly_z ~", GEO, "+", SPL),
            resp = "log", fam = "gaussian", froude = character(0), spline = TRUE, hetero = FALSE, hier = FALSE),
  M4 = list(label = "M4 M3 + heteroskedastic sigma(Froude)", f = paste("ly_z ~", GEO, "+", SPL),
            resp = "log", fam = "gaussian", froude = character(0), spline = TRUE, hetero = TRUE, hier = FALSE),
  M5 = list(label = "M5 M3 with Student-t likelihood", f = paste("ly_z ~", GEO, "+", SPL),
            resp = "log", fam = "student", froude = character(0), spline = TRUE, hetero = FALSE, hier = FALSE),
  M6 = list(label = "M6 M4 + hull random intercept (final)", f = paste("ly_z ~", GEO, "+", SPL, "+ (1 | hull)"),
            resp = "log", fam = "gaussian", froude = character(0), spline = TRUE, hetero = TRUE, hier = TRUE)
)

make_formula <- function(sp) {
  if (sp$hetero) bf(as.formula(sp$f), sigma ~ z_Froude_No) else bf(as.formula(sp$f))
}

# s = pengali skala SD prior (1 = dasar, 0.5 lebih ketat, 2.5 lebih lebar) untuk analisis sensitivitas
make_prior <- function(sp, s = 1) {
  p <- prior_string(sprintf("normal(0, %g)", 1 * s), class = "Intercept") +
       prior_string(sprintf("normal(0, %g)", 0.5 * s), class = "b")                  # geometri
  for (v in sp$froude) p <- p + prior_string(sprintf("normal(0, %g)", 1 * s), class = "b", coef = v)
  if (sp$spline) for (k in 1:6) p <- p + prior_string(sprintf("normal(0, %g)", 1.5 * s), class = "b", coef = paste0("F", k))
  if (sp$hetero) {
    p <- p + prior_string(sprintf("normal(-1.5, %g)", 0.75 * s), class = "Intercept", dpar = "sigma") +
             prior_string(sprintf("normal(0, %g)", 0.5 * s), class = "b", dpar = "sigma")
  } else {
    p <- p + prior_string(sprintf("exponential(%g)", 1 / s), class = "sigma")          # rate = 1/s
  }
  if (sp$fam == "student") p <- p + prior(gamma(2, 0.1), class = "nu")               # shape-rate
  if (sp$hier) p <- p + prior_string(sprintf("normal(0, %g)", 0.25 * s), class = "sd")  # half-normal
  p
}

fit_model <- function(name, data = mdl, s = 1, seed = SEED, tag = name, sampler = SAMPLER) {
  sp <- SPECS[[name]]
  invisible(gc())                                   # bebaskan memori sebelum kompilasi Stan
  brm(make_formula(sp), data = data,
      family = if (sp$fam == "student") student() else gaussian(),
      prior = make_prior(sp, s), chains = sampler$chains, iter = sampler$iter,
      warmup = sampler$warmup, cores = N_CORES, seed = seed, control = CONTROL,
      refresh = 0, silent = 2, file = file.path("fits", tag), file_refit = "on_change")
}

MAIN_PAR <- function(fit) grep("^b_|^sigma$|^sd_|^nu$", variables(fit), value = TRUE)

diag_row <- function(fit, name) {
  t0 <- sum(rstan::get_elapsed_time(fit$fit))
  s <- summarise_draws(as_draws_array(fit, variable = MAIN_PAR(fit)),
                       "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")
  np <- nuts_params(fit)
  data.frame(model = name, sampling_s = round(t0, 1), max_rhat = max(s$rhat),
             min_ess_bulk = min(s$ess_bulk), min_ess_tail = min(s$ess_tail),
             max_mcse_over_sd = max(s$mcse_mean / s$sd),
             divergences = sum(np$Value[np$Parameter == "divergent__"]))
}

loo_raw <- function(fit, name, data = mdl) {
  ll <- log_lik(fit)
  ll <- sweep(ll, 2, jacobian(data, SPECS[[name]]$resp, SCALE), "+")   # pindah ke satuan asli Resistance
  chain_id <- rep(seq_len(nchains(fit)), each = ndraws(fit) / nchains(fit))
  list(loo = loo(ll, r_eff = relative_eff(exp(ll), chain_id = chain_id)), waic = waic(ll))
}

# -----------------------------------------------------------------------------
# 5. PRIOR PREDICTIVE CHECK (simulasi eksplisit dari prior, materi Bag. 5.3)
# -----------------------------------------------------------------------------
prior_predictive <- function(name, S = 1000, wide = FALSE) {
  sp <- SPECS[[name]]
  set.seed(SEED)
  Xg <- scale(as.matrix(mdl[paste0("z_", GEOM)]), scale = FALSE)          # brms: intercept pada prediktor terpusat
  zF <- mdl$z_Froude_No
  out <- matrix(NA_real_, S, nrow(mdl))
  for (s in 1:S) {
    ic  <- rnorm(1, 0, if (wide) 2.5 else 1)
    mu  <- ic + Xg %*% rnorm(length(GEOM), 0, if (wide) 1 else 0.5)
    if (sp$spline) {
      Fm <- scale(as.matrix(mdl[paste0("F", 1:6)]), scale = FALSE)
      mu <- mu + Fm %*% rnorm(6, 0, if (wide) 2.5 else 1.5)
    } else {
      mu <- mu + (zF - mean(zF)) * rnorm(1, 0, 1)
    }
    if (sp$hier) {
      u <- rnorm(AUDIT$n_hull, 0, abs(rnorm(1, 0, if (wide) 1 else 0.25)))
      mu <- mu + u[mdl$hull_id]
    }
    sig <- if (sp$hetero) exp(rnorm(1, if (wide) 0 else -1.5, if (wide) 2 else 0.75) +
                              rnorm(1, 0, if (wide) 1 else 0.5) * zF) else rexp(1, 1)
    out[s, ] <- rnorm(nrow(mdl), mu, sig)
  }
  to_raw(out, sp$resp, SCALE)
}
PRIOR_PRED <- cache_rds("prior_predictive", list(
  M0 = prior_predictive("M0"), M6 = prior_predictive("M6"), M6_initial = prior_predictive("M6", wide = TRUE)))

pp_summary <- function(v) {
  v <- v[is.finite(v)]
  c(share_negative = mean(v < 0), share_below_minus_0.5 = mean(v < -0.5), share_gt_1000 = mean(v > 1000),
    q01 = unname(quantile(v, .01)), q50 = median(v), q95 = unname(quantile(v, .95)), q99 = unname(quantile(v, .99)))
}

# -----------------------------------------------------------------------------
# 6. PRIOR SALAH SKALA (reproduksi kesalahan versi awal, materi Bag. 7/10A)
# -----------------------------------------------------------------------------
fit_old <- brm(as.formula(paste("Resistance ~", paste(PRED, collapse = " + "))), data = df,
               prior = prior_string(sprintf("normal(%g, %g)", mean(df$Resistance), 2.5 * sd(df$Resistance)), class = "Intercept") +
                       prior(normal(0, 1), class = "b") +
                       prior_string(sprintf("normal(0, %g)", sd(df$Resistance)), class = "sigma"),
               chains = CV_SAMPLER$chains, iter = CV_SAMPLER$iter, warmup = CV_SAMPLER$warmup, cores = N_CORES,
               seed = SEED, refresh = 0, silent = 2, file = "fits/old_misscaled_prior", file_refit = "on_change")
Xraw <- cbind(1, as.matrix(df[PRED]))
b_ols_raw <- drop(solve(crossprod(Xraw), crossprod(Xraw, df$Resistance)))
OLD <- list(ols = b_ols_raw[["Froude_No"]], old = fixef(fit_old)["Froude_No", ],
            prior_span95 = 1.96 * diff(range(df$Froude_No)))

# -----------------------------------------------------------------------------
# 7. VALIDASI KOMPUTASI: analitik konjugat vs Gibbs manual vs HMC (Stan)
# -----------------------------------------------------------------------------
X <- cbind(Intercept = 1, as.matrix(mdl[paste0("z_", PRED)]))
y <- mdl$y_z
n <- nrow(X); p <- ncol(X)
m0 <- rep(0, p); V0 <- diag(p); a0 <- 2; b0 <- 1
V0i <- solve(V0)
Vn <- solve(V0i + crossprod(X))
mn <- drop(Vn %*% (V0i %*% m0 + crossprod(X, y)))
an <- a0 + n / 2
bn <- b0 + 0.5 * drop(crossprod(y) + t(m0) %*% V0i %*% m0 - t(mn) %*% solve(Vn) %*% mn)
sd_analytic <- sqrt(diag(Vn) * bn / (an - 1))
b_ols <- drop(solve(crossprod(X), crossprod(X, y)))

GIBBS <- cache_rds("gibbs", {
  set.seed(SEED)
  B <- 12000; burn <- 2000
  beta <- matrix(0, B, p); s2 <- rep(1, B)
  L <- t(chol(Vn))
  for (it in 2:B) {
    beta[it, ] <- mn + sqrt(s2[it - 1]) * drop(L %*% rnorm(p))
    r <- y - X %*% beta[it, ]
    dp <- beta[it, ] - m0
    s2[it] <- 1 / rgamma(1, shape = a0 + (n + p) / 2,
                         rate = b0 + 0.5 * drop(crossprod(r) + t(dp) %*% V0i %*% dp))  # 1/sigma2 ~ Gamma(shape, rate)
  }
  list(beta = beta[-(1:burn), ], s2 = s2[-(1:burn)])
})

NIG_STAN <- "
data { int<lower=1> N; int<lower=1> P; matrix[N, P] X; vector[N] y; real a0; real b0; }
parameters { real<lower=0> sigma2; vector[P] beta_raw; }
transformed parameters { vector[P] beta = sqrt(sigma2) * beta_raw; }
model { sigma2 ~ inv_gamma(a0, b0); beta_raw ~ std_normal(); y ~ normal(X * beta, sqrt(sigma2)); }
"
HMC_NIG <- cache_rds("hmc_nig", {
  sm <- stan_model(model_code = NIG_STAN)
  fit <- sampling(sm, data = list(N = n, P = p, X = X, y = y, a0 = a0, b0 = b0),
                  chains = 4, iter = 3000, warmup = 1000, seed = SEED, refresh = 0)
  summarise_draws(as_draws_array(fit), "mean", "sd", "mcse_mean", "rhat", "ess_bulk")
})
hmc_beta <- HMC_NIG[grepl("^beta\\[", HMC_NIG$variable), ]
COMP <- data.frame(
  OLS = b_ols, Analytic_mean = mn, Analytic_SD = sd_analytic, Gibbs_mean = colMeans(GIBBS$beta),
  HMC_mean = hmc_beta$mean, HMC_MCSE = hmc_beta$mcse_mean, row.names = colnames(X))
COMP$HMC_minus_Analytic_in_MCSE <- (COMP$HMC_mean - COMP$Analytic_mean) / COMP$HMC_MCSE
SIGMA2 <- c(Analytic = bn / (an - 1), Gibbs = mean(GIBBS$s2),
            HMC = HMC_NIG$mean[HMC_NIG$variable == "sigma2"])

# log-likelihood manual (materi Bag. 4) untuk diverifikasi terhadap Stan
loglik_gaussian <- function(beta, sigma, X, y) sum(dnorm(y, drop(X %*% beta), sigma, log = TRUE))

# -----------------------------------------------------------------------------
# 8. FIT MODEL M0 - M6
# -----------------------------------------------------------------------------
FITS <- lapply(setNames(names(SPECS), names(SPECS)), fit_model)
DIAG <- bind_rows(lapply(names(FITS), function(m) diag_row(FITS[[m]], m)))
DIAG$verdict <- ifelse(DIAG$max_rhat <= 1.01 & DIAG$min_ess_bulk >= 400 & DIAG$min_ess_tail >= 400 &
                         DIAG$divergences == 0, "OK", "CHECK")

# verifikasi likelihood manual pada satu posterior draw M0
d0 <- as_draws_df(FITS$M0)
LL_CHECK <- c(
  manual = loglik_gaussian(unlist(d0[1, paste0("b_", c("Intercept", paste0("z_", PRED)))]), d0$sigma[1], X, y) +
    n * (-log(SCALE$y[2])),
  stan = sum(log_lik(FITS$M0, draw_ids = 1)) + n * (-log(SCALE$y[2])))

# -----------------------------------------------------------------------------
# 9. LOO / WAIC (koreksi Jacobian ke satuan asli Resistance)
# -----------------------------------------------------------------------------
# Kompatibilitas cache loo:
# Object loo yang tersimpan dari versi package lama dapat memiliki struktur
# diagnostik Pareto-k yang berbeda. Ambil Pareto-k dari struktur yang tersedia,
# dan bila cache lama tidak menyimpannya, hitung ulang LOOS dari fit saat ini.
extract_pareto_k <- function(x) {
  if (is.null(x)) return(NULL)

  # Struktur utama pada loo modern
  if (!is.null(x$diagnostics) && !is.null(x$diagnostics$pareto_k)) {
    return(as.numeric(x$diagnostics$pareto_k))
  }

  # Fallback yang tetap tersedia pada objek loo dengan pointwise diagnostics
  if (!is.null(x$pointwise) &&
      "influence_pareto_k" %in% colnames(x$pointwise)) {
    return(as.numeric(x$pointwise[, "influence_pareto_k"]))
  }

  # Fallback terakhir bila psis object tersimpan
  if (!is.null(x$psis_object)) {
    out <- tryCatch(loo::pareto_k_values(x$psis_object),
                    error = function(e) NULL)
    if (!is.null(out)) return(as.numeric(out))
  }

  NULL
}

build_loos <- function() {
  lapply(setNames(names(FITS), names(FITS)),
         function(m) loo_raw(FITS[[m]], m))
}

LOOS <- cache_rds("loos", build_loos())

# Cache loo lama mungkin sah untuk ELPD/WAIC tetapi tidak menyimpan
# diagnostik Pareto-k yang dikenali oleh versi loo yang sedang dipakai.
# Dalam kondisi ini, cache khusus LOO dibangun ulang dari posterior yang ada.
has_pareto_k <- vapply(LOOS, function(z) {
  pk <- extract_pareto_k(z$loo)
  length(pk) > 0 && any(is.finite(pk))
}, logical(1))

if (!all(has_pareto_k)) {
  unlink(file.path("cache", "loos.rds"))
  LOOS <- cache_rds("loos", build_loos())
}

loo_cmp <- loo_compare(lapply(LOOS, `[[`, "loo"))
LOO_TAB <- data.frame(model = rownames(loo_cmp), loo_cmp[, c("elpd_loo", "se_elpd_loo", "p_loo", "elpd_diff", "se_diff")],
                      row.names = NULL)
LOO_TAB$elpd_waic <- sapply(LOO_TAB$model, function(m) LOOS[[m]]$waic$estimates["elpd_waic", "Estimate"])

LOO_TAB$n_pareto_k_gt_0.7 <- sapply(LOO_TAB$model, function(m) {
  pk <- extract_pareto_k(LOOS[[m]]$loo)
  if (is.null(pk)) return(NA_integer_)
  sum(pk > 0.7, na.rm = TRUE)
})

LOO_TAB$diff_over_se <- ifelse(LOO_TAB$se_diff > 0, LOO_TAB$elpd_diff / LOO_TAB$se_diff, NA)

# -----------------------------------------------------------------------------
# 10. PPC & RESIDUAL (satuan asli)
# -----------------------------------------------------------------------------
ppc_raw <- function(name, ndraws = 1000) {
  set.seed(SEED)
  to_raw(posterior_predict(FITS[[name]], ndraws = ndraws), SPECS[[name]]$resp, SCALE)
}
STATS <- list(mean = mean, sd = sd, max = max, skewness = skew, `share < 0` = function(v) mean(v < 0))
ppc_table <- function(yrep, yobs) {
  do.call(rbind, lapply(names(STATS), function(k) {
    rep <- apply(yrep, 1, STATS[[k]]); obs <- STATS[[k]](yobs)
    data.frame(statistic = k, observed = obs, rep_median = median(rep), rep_2.5 = unname(quantile(rep, .025)),
               rep_97.5 = unname(quantile(rep, .975)), P_rep_ge_obs = mean(rep >= obs))
  }))
}
YOBS <- mdl$Resistance
YREP <- cache_rds("yrep", list(M0 = ppc_raw("M0"), M6 = ppc_raw("M6")))
PPC <- list(M0 = ppc_table(YREP$M0, YOBS), M6 = ppc_table(YREP$M6, YOBS))

MU_M0 <- colMeans(to_raw(posterior_linpred(FITS$M0), "raw", SCALE))
MU6 <- colMeans(posterior_linpred(FITS$M6))
SIG6 <- colMeans(posterior_linpred(FITS$M6, dpar = "sigma", transform = TRUE))
RES6 <- (mdl$ly_z - MU6) / SIG6

# VIF dan korelasi posterior (materi Bag. 16, 24.4)
vif_table <- function(M) sapply(seq_len(ncol(M)), function(i) 1 / (1 - summary(lm(M[, i] ~ M[, -i]))$r.squared))
VIF <- list(geometry = setNames(vif_table(as.matrix(mdl[paste0("z_", PRED)])), PRED),
            poly = setNames(vif_table(as.matrix(mdl[c("z_Froude_No", "zF2", "zF3")])), c("Froude", "Froude^2", "Froude^3")))
POST_COR <- cor(as.data.frame(as_draws_df(FITS$M2))[, paste0("b_", c("z_LD_Ratio", "z_BD_Ratio", "z_LB_Ratio", "z_Froude_No", "zF2", "zF3"))])

# -----------------------------------------------------------------------------
# 11. POSTERIOR MANAGEMENT
# -----------------------------------------------------------------------------
hdi <- function(v, prob = .95) {
  v <- sort(v); k <- floor(prob * length(v)); w <- v[(k + 1):length(v)] - v[1:(length(v) - k)]
  i <- which.min(w); c(v[i], v[i + k])
}
kde_mode <- function(v) { d <- density(v); d$x[which.max(d$y)] }
manage <- function(v) {
  h <- hdi(v)
  c(mean = mean(v), median = median(v), mode_kde = kde_mode(v), sd = sd(v),
    eti_2.5 = unname(quantile(v, .025)), eti_97.5 = unname(quantile(v, .975)),
    hdi_low = h[1], hdi_high = h[2], P_gt_0 = mean(v > 0))
}
M0_RAW <- t(sapply(PRED, function(v) manage(d0[[paste0("b_z_", v)]] * SCALE$y[2] / SCALE[[v]][2])))
M0_RAW <- cbind(M0_RAW, OLS_raw_units = b_ols_raw[PRED])
BF <- d0$b_z_Froude_No * SCALE$y[2] / SCALE$Froude_No[2]

d6 <- as_draws_df(FITS$M6)
sl <- SCALE$ly[2]
EFF <- rbind(
  t(sapply(GEOM, function(v) { e <- 100 * (exp(d6[[paste0("b_z_", v)]] * sl) - 1); c(manage(e), P_abs_gt_1pct = mean(abs(e) > 1)) })),
  `Residual SD multiplier per +1 SD Froude` = { e <- exp(d6$b_sigma_z_Froude_No); c(manage(e), P_abs_gt_1pct = mean(e > 1.01)) },
  `Between-hull variation (%)` = { e <- 100 * (exp(d6$sd_hull__Intercept * sl) - 1); c(manage(e), P_abs_gt_1pct = mean(e > 1)) })

# kurva sigma(Froude) untuk hull tipikal
grid_df <- data.frame(Froude_No = seq(min(df$Froude_No), max(df$Froude_No), length.out = 60))
for (v in GEOM) grid_df[[v]] <- mean(df[[v]])
grid_df$Resistance <- 0; grid_df$hull_id <- 999
GRID <- FR$tf(grid_df)
SIG_GRID <- posterior_linpred(FITS$M6, newdata = GRID, dpar = "sigma", transform = TRUE, re_formula = NA)

# -----------------------------------------------------------------------------
# 12. PREDIKSI HULL BARU (posterior predictive, materi Bag. 24.6)
# -----------------------------------------------------------------------------
hull_geom <- df[!duplicated(df$hull_id), GEOM]
FROUDE_GRID <- sort(unique(df$Froude_No))
scen <- data.frame(Froude_No = FROUDE_GRID)
for (v in GEOM) scen[[v]] <- median(hull_geom[[v]])
scen$Resistance <- 0; scen$hull_id <- 999
SCEN <- FR$tf(scen)
predict_new_hull <- function(fit) {
  set.seed(SEED)
  to_raw(posterior_predict(fit, newdata = SCEN, allow_new_levels = TRUE, sample_new_levels = "gaussian"), "log", SCALE)
}
PRED_NEW <- predict_new_hull(FITS$M6)
MEAN_TYP <- to_raw(posterior_linpred(FITS$M6, newdata = SCEN, re_formula = NA), "log", SCALE)
PRED_TAB <- data.frame(Froude = FROUDE_GRID, pred_median = apply(PRED_NEW, 2, median),
                       PI_2.5 = apply(PRED_NEW, 2, quantile, .025), PI_97.5 = apply(PRED_NEW, 2, quantile, .975),
                       P_R_gt_5 = colMeans(PRED_NEW > 5), P_R_gt_20 = colMeans(PRED_NEW > 20))
I40 <- which(abs(FROUDE_GRID - 0.40) < 1e-9)

# -----------------------------------------------------------------------------
# 13. SENSITIVITAS: lab prior-data (analitik) & prior pada M6 (MCMC)
# -----------------------------------------------------------------------------
sigma2_known <- bn / (an - 1)
hull_order <- unique(mdl$hull_id)
SUBSETS <- list(`1 hull (n=14)` = hull_order[1], `3 hulls (n=42)` = hull_order[1:3], `22 hulls (n=308)` = hull_order)
PRIORS_LAB <- list(`Weak N(0, 1)` = c(0, 1), `Strong-correct N(0.81, 0.02)` = c(0.81, 0.02),
                   `Strong-wrong N(0.40, 0.02)` = c(0.40, 0.02))
LAB <- do.call(rbind, lapply(names(SUBSETS), function(sn) {
  idx <- mdl$hull_id %in% SUBSETS[[sn]]
  Xs <- X[idx, ]; ys <- y[idx]
  do.call(rbind, lapply(names(PRIORS_LAB), function(pn) {
    pr <- PRIORS_LAB[[pn]]
    m0v <- rep(0, p); m0v[p] <- pr[1]
    V0v <- diag(p); V0v[p, p] <- pr[2]^2
    Vni <- solve(V0v) + crossprod(Xs) / sigma2_known
    Vn_ <- solve(Vni)
    mn_ <- drop(Vn_ %*% (solve(V0v) %*% m0v + crossprod(Xs, ys) / sigma2_known))
    b_hat <- coef(lm(ys ~ Xs[, p]))[2]
    data.frame(data = sn, prior = pn, post_mean = mn_[p], post_sd = sqrt(Vn_[p, p]),
               data_only_slope = unname(b_hat), prior_pull = (1 / pr[2]^2) / Vni[p, p])
  }))
}))

SENS_FITS <- list(`Tighter (x0.5)` = fit_model("M6", s = 0.5, tag = "M6_prior_x0.5"),
                  `Base (x1)` = FITS$M6,
                  `Wider (x2.5)` = fit_model("M6", s = 2.5, tag = "M6_prior_x2.5"))
KEY <- c("b_z_LD_Ratio", "b_z_BD_Ratio", "b_z_LB_Ratio", "b_z_LCB", "b_z_Prismatic_Coef",
         "b_sigma_z_Froude_No", "sd_hull__Intercept")
SENS <- cache_rds("sensitivity", do.call(rbind, lapply(names(SENS_FITS), function(lab) {
  f <- SENS_FITS[[lab]]; dd <- as_draws_df(f)
  data.frame(prior = lab, param = KEY, median = sapply(KEY, function(k) median(dd[[k]])),
             lo = sapply(KEY, function(k) quantile(dd[[k]], .025)), hi = sapply(KEY, function(k) quantile(dd[[k]], .975)),
             P_R20_Fr040 = mean(predict_new_hull(f)[, I40] > 20),
             ELPD_LOO = loo_raw(f, "M6")$loo$estimates["elpd_loo", "Estimate"], row.names = NULL)
})))
NU <- as_draws_df(FITS$M5)$nu

# -----------------------------------------------------------------------------
# 14. AUDIT REPRODUCIBILITY: seed & urutan data
# -----------------------------------------------------------------------------
compare_runs <- function(fa, fb) {
  v <- MAIN_PAR(fa)
  sa <- summarise_draws(as_draws_array(fa, variable = v), "mean", "sd", "mcse_mean")
  sb <- summarise_draws(as_draws_array(fb, variable = v), "mean", "sd", "mcse_mean")
  data.frame(param = v, mean_a = sa$mean, mean_b = sb$mean, sd_a = sa$sd,
             z_mcse = (sb$mean - sa$mean) / sqrt(sa$mcse_mean^2 + sb$mcse_mean^2),
             pct_of_sd = 100 * abs(sb$mean - sa$mean) / sa$sd,
             expected_noise_pct = 100 * sqrt(2) * median(sa$mcse_mean / sa$sd))
}
SEEDS_EXTRA <- c(777, 999, 31, 123)
SEED_FITS <- c(list(`2026` = FITS$M6),
               lapply(setNames(SEEDS_EXTRA, SEEDS_EXTRA), function(s)
                 cache_rds(paste0("fit_M6_seed", s), update(FITS$M6, seed = s, refresh = 0, silent = 2, recompile = FALSE))))
FIT_SEED <- SEED_FITS[["777"]]
set.seed(SEED)
SHUFFLED <- mdl[sample(nrow(mdl)), ]
FIT_SHUF <- cache_rds("fit_M6_shuffled", update(FITS$M6, newdata = SHUFFLED, seed = SEED, refresh = 0,
                                                silent = 2, recompile = FALSE))
CMP_SEED <- compare_runs(FITS$M6, FIT_SEED)
CMP_ORDER <- compare_runs(FITS$M6, FIT_SHUF)
FLAGGED <- unique(c(CMP_SEED$param[abs(CMP_SEED$z_mcse) >= 3], CMP_ORDER$param[abs(CMP_ORDER$z_mcse) >= 3]))

# Audit lima seed: apakah MCSE yang dilaporkan sesuai dengan variasi nyata antar-run?
v_main <- MAIN_PAR(FITS$M6)
seed_sum <- lapply(SEED_FITS, function(f) summarise_draws(as_draws_array(f, variable = v_main), "mean", "mcse_mean"))
seed_means <- sapply(seed_sum, `[[`, "mean"); rownames(seed_means) <- v_main
MULTI_SEED <- data.frame(param = v_main, t(apply(seed_means, 1, range)),
                         between_run_sd = apply(seed_means, 1, sd),
                         mean_reported_mcse = rowMeans(sapply(seed_sum, `[[`, "mcse_mean")))
names(MULTI_SEED)[2:3] <- c("min_mean_5_seeds", "max_mean_5_seeds")
MULTI_SEED$ratio_sd_to_mcse <- MULTI_SEED$between_run_sd / MULTI_SEED$mean_reported_mcse
# Kesimpulan substantif antar seed: variasi antar-hull (%) dan probabilitas eksedans hull baru
SEED_CONCL <- do.call(rbind, lapply(names(SEED_FITS), function(s) {
  dd <- as_draws_df(SEED_FITS[[s]])
  data.frame(seed = s, between_hull_pct = median(100 * (exp(dd$sd_hull__Intercept * SCALE$ly[2]) - 1)),
             P_R20_Fr040 = mean(predict_new_hull(SEED_FITS[[s]])[, I40] > 20))
}))

# -----------------------------------------------------------------------------
# 15. CROSS-VALIDATION: per hull (grouped) vs acak per baris
# -----------------------------------------------------------------------------
logmeanexp <- function(x) { m <- max(x); m + log(mean(exp(x - m))) }
make_folds <- function(scheme) {
  if (scheme == "grouped") {
    hulls <- sort(unique(full$hull_id))
    fold_of_hull <- rep(1:5, length.out = length(hulls))            # 22 hull -> fold berisi 4-5 hull utuh
    fold_of_hull[match(full$hull_id, hulls)]
  } else {
    set.seed(SEED); sample(rep(1:5, length.out = nrow(full)))
  }
}
cv_eval <- function(name, scheme) cache_rds(paste0("cv_", name, "_", scheme), {
  sp <- SPECS[[name]]
  resp <- if (sp$resp == "raw") "y_z" else "ly_z"
  folds <- make_folds(scheme)
  recs <- list(); rhats <- c()
  for (k in 1:5) {
    fr <- make_frames(full[folds != k, ], full[folds == k, ])
    f <- update(FITS[[name]], newdata = fr$train, chains = CV_SAMPLER$chains, iter = CV_SAMPLER$iter,
                warmup = CV_SAMPLER$warmup, seed = SEED + k, refresh = 0, silent = 2, recompile = FALSE)
    rhats <- c(rhats, max(summarise_draws(as_draws_array(f, variable = MAIN_PAR(f)), "rhat")$rhat))
    ll <- log_lik(f, newdata = fr$test, allow_new_levels = TRUE, sample_new_levels = "gaussian")
    lpd <- apply(ll, 2, logmeanexp) + jacobian(fr$test, sp$resp, fr$sc)
    set.seed(SEED + k)
    yr <- to_raw(posterior_predict(f, newdata = fr$test, allow_new_levels = TRUE, sample_new_levels = "gaussian"),
                 sp$resp, fr$sc)
    recs[[k]] <- data.frame(fold = k, hull_id = fr$test$hull_id, Froude_No = fr$test$Froude_No,
                            y = fr$test$Resistance, pred_mean = colMeans(yr),
                            lo = apply(yr, 2, quantile, .025), hi = apply(yr, 2, quantile, .975), lpd = lpd)
  }
  r <- bind_rows(recs)
  list(pred = r, summary = data.frame(
    model = name, scheme = scheme, RMSE = sqrt(mean((r$y - r$pred_mean)^2)), MAE = mean(abs(r$y - r$pred_mean)),
    coverage_95 = mean(r$y >= r$lo & r$y <= r$hi), PI_width = mean(r$hi - r$lo),
    ELPD_per_obs = mean(r$lpd), ELPD_total = sum(r$lpd), max_rhat_folds = max(rhats)))
})
CV_PLAN <- list(c("M0", "grouped"), c("M1", "grouped"), c("M2", "grouped"), c("M3", "grouped"), c("M4", "grouped"),
                c("M6", "grouped"), c("M4", "row-random"), c("M6", "row-random"))
CV <- lapply(CV_PLAN, function(x) cv_eval(x[1], x[2]))
names(CV) <- sapply(CV_PLAN, paste, collapse = "_")
CV_TAB <- bind_rows(lapply(CV, `[[`, "summary"))

cv_key <- function(r) paste(r$hull_id, round(r$Froude_No, 3))
grouped_models <- CV_TAB$model[CV_TAB$scheme == "grouped"]
BEST_CV <- grouped_models[which.max(CV_TAB$ELPD_total[CV_TAB$scheme == "grouped"])]
ref <- CV[[paste0(BEST_CV, "_grouped")]]$pred
PAIRED <- bind_rows(lapply(grouped_models, function(m) {
  o <- CV[[paste0(m, "_grouped")]]$pred
  o <- o[match(cv_key(ref), cv_key(o)), ]
  d <- o$lpd - ref$lpd
  hs <- tapply(d, ref$hull_id, sum)
  data.frame(model = m, ELPD_total = sum(o$lpd), elpd_diff = sum(d),
             se_diff_obs = sd(d) * sqrt(length(d)), se_diff_hull = sd(hs) * sqrt(length(hs)))
})) %>% arrange(desc(ELPD_total)) %>% mutate(diff_over_se_hull = ifelse(se_diff_hull > 0, elpd_diff / se_diff_hull, NA))
OPTIMISM <- bind_rows(lapply(c("M4", "M6"), function(m) {
  g <- CV[[paste0(m, "_grouped")]]$pred; r <- CV[[paste0(m, "_row-random")]]$pred
  r <- r[match(cv_key(g), cv_key(r)), ]
  d <- r$lpd - g$lpd
  data.frame(model = m, ELPD_optimism = sum(d), se = sd(d) * sqrt(length(d)),
             RMSE_grouped = sqrt(mean((g$y - g$pred_mean)^2)), RMSE_row = sqrt(mean((r$y - r$pred_mean)^2)))
})) %>% mutate(optimism_over_se = ELPD_optimism / se)

# -----------------------------------------------------------------------------
# 16. CHECKLIST & MANIFEST REPRODUCIBILITY
# -----------------------------------------------------------------------------
CHECKS <- data.frame(Item = c(
  "Seed eksplisit tercatat", "Checksum SHA-256 data mentah tercatat", "Preprocessing scripted (make_frames)",
  "Standardisasi dihitung di dalam fold (tanpa kebocoran)", "Prior eksplisit + prior predictive check",
  "Komputasi divalidasi (analitik vs Gibbs vs HMC, < 3 MCSE)", "Likelihood manual = likelihood Stan",
  "Semua model: R-hat <= 1.01, ESS >= 400, 0 divergence", "Stabil terhadap seed (semua |z| < 3 MCSE)",
  "Stabil terhadap urutan data (semua |z| < 3 MCSE)",
  "Kesimpulan prediktif stabil antar 5 seed (rentang P(R>20 | Fr=0.40) < 0.05)",
  "Validasi per-hull dipakai untuk klaim prediktif", "Objek fit tersimpan (fits/*.rds)"),
  Pass = c(TRUE, nchar(RAW_SHA256) == 64, TRUE, TRUE, TRUE, max(abs(COMP$HMC_minus_Analytic_in_MCSE)) < 3,
           isTRUE(all.equal(LL_CHECK[["manual"]], LL_CHECK[["stan"]])), all(DIAG$verdict == "OK"),
           all(abs(CMP_SEED$z_mcse) < 3), all(abs(CMP_ORDER$z_mcse) < 3),
           diff(range(SEED_CONCL$P_R20_Fr040)) < 0.05, TRUE, file.exists("fits/M6.rds")))

# --------------------------------------------------------------------------
# Export CSV yang kompatibel dengan semua versi R:
# write.table()/write.csv() tidak dapat menulis kolom bertipe list.
# Hasil analisis tetap sama; bagian ini hanya mengubah cara file CSV disimpan.
# --------------------------------------------------------------------------
safe_csv <- function(x, file) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)

  for (j in seq_along(x)) {
    if (is.factor(x[[j]])) {
      x[[j]] <- as.character(x[[j]])
    } else if (is.list(x[[j]])) {
      x[[j]] <- vapply(
        x[[j]],
        function(v) {
          if (length(v) == 0L || all(is.na(v))) return("")
          paste(as.character(unlist(v, use.names = FALSE)), collapse = "; ")
        },
        character(1)
      )
    }
  }

  utils::write.csv(x, file, row.names = FALSE, na = "")
}

safe_csv(CV_TAB, "cv_results.csv")
safe_csv(LOO_TAB, "model_comparison_loo.csv")
MANIFEST <- list(
  seed = SEED, raw_data_sha256 = RAW_SHA256, R = R.version.string,
  packages = sapply(c("brms", "rstan", "loo", "posterior", "dplyr", "splines"), function(x) as.character(packageVersion(x))),
  sampler = SAMPLER, cv_sampler = CV_SAMPLER, control = CONTROL,
  scaling_full_data = SCALE, diagnostics = DIAG, loo = LOO_TAB[, c("model", "elpd_loo", "elpd_diff", "se_diff")],
  cv = CV_TAB, multi_seed = MULTI_SEED, seed_conclusions = SEED_CONCL, checks = CHECKS)
jsonlite::write_json(MANIFEST, "reproducibility_manifest_R.json", pretty = TRUE, auto_unbox = TRUE, digits = NA)
