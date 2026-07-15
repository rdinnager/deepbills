###############################################################################
## cv_score.R — held-out-tip cross-validation scorer (tip-cv-design.md §2, §6).
##
## For a MASKED fit (train_v3_param.R with MASK_FRAC>0), simulate each held-out
## tip forward from its parent ancestral state under
##   (a) the Riemannian manifold BM (Euler-Maruyama; diagonal metric G), and
##   (b) plain Euclidean BM  (analytic N(z_p, diag(sigma^2) * t)),
## and score the observed tip against the simulated predictive cloud with the
## ENERGY SCORE (strictly proper, needs only samples -> the intractable rho/lambda
## normaliser never appears). Score in latent space AND decoded beak-shape space.
##
## Lower energy score = better. The manifold-vs-Euclidean gap (esp. on long
## terminal branches) is the geometry test (R5). rho, lambda selected by picking
## the (rho,lambda) grid point with the best held-out score OUTSIDE this script.
##
## Run:  SCORE_TAG=<masked_out_tag> [N_DRAWS=500 N_STEPS=200 DRIFT=on] Rscript R/cv_score.R
###############################################################################
setwd(Sys.getenv("DEEPBILLS_DIR", unset = "~/scratch/deepbills"))
suppressMessages({library(tidyverse); library(ape); library(phyf); library(torch)})
select <- dplyr::select

TAG     <- Sys.getenv("SCORE_TAG");            stopifnot(nzchar(TAG))
N_DRAWS <- as.integer(Sys.getenv("N_DRAWS", "500"))
N_STEPS <- as.integer(Sys.getenv("N_STEPS", "200"))
## DRIFT: include the Onsager-Machlup drift term in the manifold Euler-Maruyama step.
## Default OFF: correct-with-drift is IMPRACTICAL in R-torch — even after collapsing the
## Jacobian to one d(logdetG) backward + a vectorised finite-diff diagonal, each step
## issues many small torch ops and the run is per-op-CALL-overhead bound (GPU sits at ~2%
## util; a 200-draw/25-step smoke case did not finish in 5+ min). Exact drift needs the
## simulator ported to Python/libtorch (no R per-op overhead), or use the analytic
## heat-kernel score (tip-cv-design.md §4), which carries the same 1/2 logdetG volume term
## WITHOUT stepwise simulation. DRIFT=off = driftless metric-scaled diffusion: fast,
## validated, and defensible on the short terminal branches this CV uses.
DRIFT   <- tolower(Sys.getenv("DRIFT", "off")) %in% c("on","1","true","yes")
dev     <- if (cuda_is_available()) "cuda" else "cpu"
rundir  <- file.path("data", TAG)
cat(sprintf("=== CV scoring %s  (draws=%d steps=%d drift=%s dev=%s) ===\n",
            TAG, N_DRAWS, N_STEPS, DRIFT, dev))

## ---- run-specific metric hyperparameters (same source as score_run.R) ----
conv   <- read.csv(file.path(rundir, "convergence.csv"))
lambda <- as.numeric(conv$lambda[1])
rho    <- as.numeric(readRDS(file.path(rundir, "estimated_rho_16dim_rho_schedule_v3.rds")))
cat(sprintf("rho=%.5f  lambda=%g\n", rho, lambda))

## ---- centroids / vars / latvars (replicate training pivots; from score_run.R) ----
bill_z  <- read_csv("data/bills_vae_latent_codes_16dim.csv", show_col_types = FALSE)
cent_df <- bill_z %>% select(Species, latent_dim, latent_mean) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_mean)
var_df  <- bill_z %>% select(Species, latent_dim, latent_var) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_var)
cent_order <- setdiff(names(cent_df), "Species")
C  <- as.matrix(cent_df[, cent_order]); V <- as.matrix(var_df[, cent_order])
bill_pf <- readRDS("data/bill_pf_16dim.rds")
bill_tips_mat <- bill_pf %>% filter(is_tip) %>% select(starts_with("latent_")) %>% as.matrix()
latvars <- apply(bill_tips_mat, 2, sd)[cent_order]

## torch constants on device
Ct  <- torch_tensor(C,  dtype = torch_float(), device = dev)   # 2021 x 16
Vt  <- torch_tensor(V,  dtype = torch_float(), device = dev)
invVt <- torch_tensor(1/V, dtype = torch_float(), device = dev)
lvt <- torch_tensor(latvars, dtype = torch_float(), device = dev)  # 16

## diagonal Riemannian metric G(Z): Z (m x 16) -> G (m x 16), matches score_run.R metric_G.
## Accumulate the Mahalanobis distance by looping over the 16 dims so the peak intermediate
## is m x 2021, NOT m x 2021 x 16 — the 3-D form OOM/thrashes at m = n_held * N_DRAWS (~10 GB).
metric_G <- function(Z) {
  m <- Z$size(1)
  dmah2 <- torch_zeros(c(m, nrow(C)), device = dev)     # m x 2021
  for (j in seq_len(16)) {
    Zj <- Z[, j]$unsqueeze(2)                            # m x 1
    diff2 <- (Zj - Ct[, j]$unsqueeze(1))$pow(2)          # m x 2021
    den   <- Zj$pow(2) + Vt[, j]$unsqueeze(1) + latvars[[j]]
    dmah2 <- dmah2 + diff2 / den
  }
  mh <- torch_exp(-dmah2 / (rho^2))                      # m x 2021
  1 / (lambda + torch_matmul(mh, invVt))                 # m x 16
}

## ---- held-out tips + parent ancestral states + sigma^2 ----
cv  <- readRDS(file.path(rundir, "cv_holdout.rds"))
r3  <- readRDS(file.path(rundir, "bill_vae_aces_16dim_v3_bayesian.rds"))
lat_cols <- paste0("latent_", 1:16)
edge_of  <- setNames(seq_len(nrow(r3)), r3$edge)

zp_of <- function(lab) {                       # parent ancestral state = path start
  as.numeric(as.matrix(r3$z_seqs[[edge_of[[lab]]]][1, lat_cols]))
}
Zp   <- t(vapply(cv$held_labels, zp_of, numeric(16)))   # n_held x 16 (cent_order == latent_1..16)
tvec <- cv$term_blen                                    # n_held
Yobs <- as.matrix(cv$obs_latent)[, lat_cols, drop = FALSE]  # n_held x 16 observed
n_held <- nrow(Zp)
stopifnot(ncol(Zp) == 16, length(tvec) == n_held)

## sigma^2 per dim from TRAINING (non-held) edge increments: incr = z_end - z_start, Var=sigma^2*blen
blens <- pf_mean_edge_features(bill_pf$phlo); names(blens) <- pf_edge_names(bill_pf$phlo)
blens <- blens[bill_pf$label]
held_set <- cv$held_labels
incr <- matrix(NA_real_, nrow(r3), 16); bl <- numeric(nrow(r3))
for (i in seq_len(nrow(r3))) {
  lab <- r3$edge[i]; if (lab %in% held_set) next
  zs <- as.matrix(r3$z_seqs[[i]][, lat_cols]); incr[i, ] <- zs[nrow(zs), ] - zs[1, ]
  bl[i] <- if (!is.na(blens[lab])) blens[lab] else NA_real_
}
ok <- is.finite(bl) & bl > 1e-8 & is.finite(rowSums(incr))
## PER-MODEL, metric-consistent sigma^2 (scalar; matches cv_heatkernel.R). Sharing a
## flat sigma^2 mis-calibrates the manifold model's metric-scaled predictive spread and
## unfairly penalises it. Euclidean: flat increments. Manifold: training-edge geodesic
## energies (fitted z_seqs as near-geodesic proxy). d=16, blen normalisation.
s2_euc <- mean(rowSums(incr[ok, ]^2) / (16 * bl[ok]))
train_idx <- which(!(r3$edge %in% held_set) & is.finite(bl) & bl > 1e-8)
metric_path_energy <- function(idx_chunk) {
  paths <- lapply(idx_chunk, function(i) as.matrix(r3$z_seqs[[i]][, lat_cols]))
  npt <- nrow(paths[[1]]); k <- length(paths)
  P <- torch_tensor(array(unlist(lapply(paths, t)), dim = c(16, npt, k)),
                    dtype = torch_float(), device = dev)$permute(c(3,1,2))
  vec  <- P[ , , 2:npt] - P[ , , 1:(npt-1)]
  ymid <- (P[ , , 2:npt] + P[ , , 1:(npt-1)]) / 2
  G <- metric_G(ymid$permute(c(1,3,2))$reshape(c(k*(npt-1), 16)))$reshape(c(k, npt-1, 16))$permute(c(1,3,2))
  as.numeric(((vec * G * vec)$sum(dim = 2)$sum(dim = 2) * (npt-1))$cpu())
}
e_tr <- numeric(0)
for (st in seq(1, length(train_idx), by = 500)) {
  en <- min(st + 499, length(train_idx)); e_tr <- c(e_tr, metric_path_energy(train_idx[st:en]))
}
s2_G <- mean(e_tr / (16 * bl[train_idx]))
cat(sprintf("sigma^2: euc=%.4g  manifold(G)=%.4g; n_edges=%d\n", s2_euc, s2_G, sum(ok)))
sig_euc <- sqrt(s2_euc); sig_man <- sqrt(s2_G)          # scalars

###############################################################################
## Simulate N draws for all held tips, in ROW CHUNKS (bounds GPU memory: metric_G
## and the drift autograd both materialise chunk x 2021 intermediates, which OOM at
## the full M = n_held*N_DRAWS). CHUNK rows per pass.
###############################################################################
CHUNK <- as.integer(Sys.getenv("CHUNK", if (DRIFT) "4000" else "8000"))
FD_EPS <- as.numeric(Sys.getenv("FD_EPS", "1e-3"))
eye16 <- torch_eye(16, device = dev)                # for batched finite-diff

## Riemannian-BM drift b^i = (1/G_ii)*(0.5*sum_k dlogG_k/dz_i - dlogG_i/dz_i).
## sum_k dlogG_k/dz_i == d(sum logG)/dz_i (rows independent) -> ONE backward.
## The diagonal dlogG_i/dz_i via batched central finite-difference (2 metric evals,
## no autograd graph). Replaces the old 16-backward-passes-per-step Jacobian.
drift_b <- function(Zc) {
  b <- Zc$size(1)
  Zg <- Zc$requires_grad_(TRUE)
  sumk <- autograd_grad(torch_log(metric_G(Zg))$sum(), Zg)[[1]]    # b x 16
  Zd <- Zc$detach()
  ## perturb dim i by +/-eps in block i: (16 x b x 16) -> (16b x 16)
  Zt <- Zd$unsqueeze(1)$expand(c(16, b, 16))
  E  <- (eye16 * FD_EPS)$view(c(16, 1, 16))
  lp <- torch_log(metric_G((Zt + E)$reshape(c(16*b, 16))))$reshape(c(16, b, 16))
  lm <- torch_log(metric_G((Zt - E)$reshape(c(16*b, 16))))$reshape(c(16, b, 16))
  ## diagonal: dp[j,i] = lp[i,j,i]  == torch_diagonal over the two size-16 axes (one op)
  dp <- torch_diagonal(lp, dim1 = 1, dim2 = 3)   # b x 16
  dm <- torch_diagonal(lm, dim1 = 1, dim2 = 3)
  diag_ik <- (dp - dm) / (2 * FD_EPS)
  Gd <- metric_G(Zd)
  list(b = (0.5 * sumk - diag_ik) / Gd, G = Gd)
}

## sg = scalar per-model BM rate (sig_man for manifold, sig_euc for Euclidean).
sim_chunk <- function(Zc, dtc, manifold, sg) {      # Zc: b x 16, dtc: b x 1
  sqrt_dt <- dtc$sqrt()
  for (s in seq_len(N_STEPS)) {
    xi <- torch_randn(Zc$size(), device = dev)
    if (!manifold) { Zc <- Zc + sg * sqrt_dt * xi; next }
    if (DRIFT) {
      db <- drift_b(Zc$detach())
      Zc <- Zc + 0.5 * (sg^2) * db$b * dtc + sg * torch_rsqrt(db$G) * sqrt_dt * xi
    } else {
      Zc <- Zc + sg * torch_rsqrt(metric_G(Zc)) * sqrt_dt * xi
    }
  }
  Zc$detach()
}

simulate <- function(manifold = TRUE) {
  sg <- if (manifold) sig_man else sig_euc
  Z0 <- Zp[rep(seq_len(n_held), each = N_DRAWS), , drop = FALSE]     # M x 16
  tt <- rep(tvec, each = N_DRAWS)                                     # M
  M  <- nrow(Z0); out <- matrix(0, M, 16)
  starts <- seq(1, M, by = CHUNK)
  for (st in starts) {
    en <- min(st + CHUNK - 1, M)
    Zc  <- torch_tensor(Z0[st:en, , drop = FALSE], dtype = torch_float(), device = dev)
    dtc <- torch_tensor(tt[st:en] / N_STEPS, dtype = torch_float(), device = dev)$unsqueeze(2)
    out[st:en, ] <- as.matrix(sim_chunk(Zc, dtc, manifold, sg)$cpu())
  }
  out                                               # M x 16
}

## energy score per tip: ES = mean_i||X_i - y|| - 0.5 * mean_pairs||X_a - X_b||
## (second term via N random pairings — unbiased, cheap)
energy_score_block <- function(X, Y) {   # X: (n_held*N) x d draws;  Y: n_held x d obs
  d <- ncol(Y); es <- numeric(n_held)
  for (h in seq_len(n_held)) {
    idx <- ((h-1)*N_DRAWS + 1):(h*N_DRAWS)
    Xi  <- X[idx, , drop = FALSE]; y <- Y[h, ]
    t1  <- mean(sqrt(rowSums((Xi - matrix(y, N_DRAWS, d, byrow = TRUE))^2)))
    p   <- sample.int(N_DRAWS)                       # random pairing
    t2  <- mean(sqrt(rowSums((Xi - Xi[p, , drop = FALSE])^2)))
    es[h] <- t1 - 0.5 * t2
  }
  es
}

set.seed(1)
cat("simulating manifold...\n"); Xm <- simulate(manifold = TRUE)
cat("simulating euclidean...\n"); Xe <- simulate(manifold = FALSE)

## ---- latent-space energy score ----
es_lat_m <- energy_score_block(Xm, Yobs)
es_lat_e <- energy_score_block(Xe, Yobs)

## ---- decoded beak-shape energy score (frozen decoder; primary space) ----
decode64 <- function(Z16) {                        # Z16: n x 16 -> n x 64 standardized shape codes
  active_dims <- readRDS("data/active_dims_16dim.rds")
  vae <- torch_load("data/bill_vae_w_trophic_v1.to"); dc <- vae$decoder; dc$eval()
  z64 <- matrix(0, nrow(Z16), 64); z64[, active_dims] <- Z16
  out <- with_no_grad(as.matrix(dc(torch_tensor(z64, dtype = torch_float()))$out_codes$cpu()))
  out
}
cat("decoding draws + observed to shape space...\n")
Sm  <- decode64(Xm); Se <- decode64(Xe); Sy <- decode64(Yobs)
es_shape_m <- energy_score_block(Sm, Sy)
es_shape_e <- energy_score_block(Se, Sy)

## ---- summary + per-tip (for branch-length stratification) ----
per_tip <- tibble(tag = TAG, label = cv$held_labels, term_blen = tvec,
                  es_lat_manifold = es_lat_m, es_lat_euclid = es_lat_e,
                  es_shape_manifold = es_shape_m, es_shape_euclid = es_shape_e)
write_csv(per_tip, file.path(rundir, "cv_per_tip.csv"))

summ <- tibble(
  tag = TAG, rho = rho, lambda = lambda, n_held = n_held,
  es_lat_manifold = mean(es_lat_m), es_lat_euclid = mean(es_lat_e),
  es_shape_manifold = mean(es_shape_m), es_shape_euclid = mean(es_shape_e),
  # win = fraction of held tips where manifold beats Euclidean (lower ES)
  win_lat = mean(es_lat_m < es_lat_e), win_shape = mean(es_shape_m < es_shape_e)
)
write_csv(summ, file.path(rundir, "cv_score.csv"))
cat("\n=== CV RESULT ===\n")
cat(sprintf("latent : manifold %.4f  vs euclid %.4f   (manifold wins %.0f%% of tips)\n",
            summ$es_lat_manifold, summ$es_lat_euclid, 100*summ$win_lat))
cat(sprintf("shape  : manifold %.4f  vs euclid %.4f   (manifold wins %.0f%% of tips)  <- primary\n",
            summ$es_shape_manifold, summ$es_shape_euclid, 100*summ$win_shape))
cat(sprintf("lower is better; manifold better if < euclid. wrote %s/cv_score.csv + cv_per_tip.csv\n", rundir))
