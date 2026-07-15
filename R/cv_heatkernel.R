###############################################################################
## cv_heatkernel.R — ANALYTIC held-out-tip CV score (tip-cv-design.md §4, §6).
##
## Short-time heat-kernel density (Minakshisundaram-Pleijel / Varadhan), the
## drift-aware alternative to the Monte-Carlo energy score that needs NO stepwise
## simulation (so no R-torch per-op-overhead wall):
##
##   log p_t(z_p, z_i) ~= -d2_G/(2 s2 t) - (d/2) log(2 pi s2 t) + 1/2 log det G(z_i)
##
## The 1/2 log det G VOLUME term is exactly the normaliser the naive energy
## comparison omits (why that comparison is circular). Van Vleck-Morette ~= 1 for
## short terminal branches. d = 16, s2 = scalar BM rate (mean of per-dim), t = blen.
##
## d2_G(z_p, z_i) = PINNED geodesic energy: optimise the cubic-path curvature (a,b)
## with BOTH endpoints fixed (parent, observed tip). Euclidean baseline uses the
## straight-line ||z_i - z_p||^2 and log det I = 0.
##
## Higher log p = better. Manifold beats Euclidean if logp_manifold > logp_euclid.
## Run:  SCORE_TAG=<masked_tag> [GEO_ITERS=400 N_SEGS=50] Rscript R/cv_heatkernel.R
###############################################################################
setwd(Sys.getenv("DEEPBILLS_DIR", unset = "~/scratch/deepbills"))
suppressMessages({library(tidyverse); library(phyf); library(torch)})
select <- dplyr::select

TAG       <- Sys.getenv("SCORE_TAG"); stopifnot(nzchar(TAG))
GEO_ITERS <- as.integer(Sys.getenv("GEO_ITERS", "400"))
NSEG      <- as.integer(Sys.getenv("N_SEGS", "50"))
dev       <- if (cuda_is_available()) "cuda" else "cpu"
rundir    <- file.path("data", TAG)
cat(sprintf("=== heat-kernel CV %s  (geo_iters=%d nseg=%d dev=%s) ===\n", TAG, GEO_ITERS, NSEG, dev))

conv   <- read.csv(file.path(rundir, "convergence.csv"))
lambda <- as.numeric(conv$lambda[1])
rho    <- as.numeric(readRDS(file.path(rundir, "estimated_rho_16dim_rho_schedule_v3.rds")))

## centroids / vars / latvars (same pivots as score_run.R / cv_score.R)
bill_z  <- read_csv("data/bills_vae_latent_codes_16dim.csv", show_col_types = FALSE)
cent_df <- bill_z %>% select(Species, latent_dim, latent_mean) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_mean)
var_df  <- bill_z %>% select(Species, latent_dim, latent_var) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_var)
cent_order <- setdiff(names(cent_df), "Species")
C <- as.matrix(cent_df[, cent_order]); V <- as.matrix(var_df[, cent_order])
bill_pf <- readRDS("data/bill_pf_16dim.rds")
bill_tips_mat <- bill_pf %>% filter(is_tip) %>% select(starts_with("latent_")) %>% as.matrix()
latvars <- apply(bill_tips_mat, 2, sd)[cent_order]

Ct <- torch_tensor(C, dtype = torch_float(), device = dev)
Vt <- torch_tensor(V, dtype = torch_float(), device = dev)
invVt <- torch_tensor(1/V, dtype = torch_float(), device = dev)

## diagonal metric G(Z): Z (m x 16) -> G (m x 16). Dim-looped (bounded memory).
metric_G <- function(Z) {
  m <- Z$size(1); dmah2 <- torch_zeros(c(m, nrow(C)), device = dev)
  for (j in seq_len(16)) {
    Zj <- Z[, j]$unsqueeze(2)
    diff2 <- (Zj - Ct[, j]$unsqueeze(1))$pow(2)
    den   <- Zj$pow(2) + Vt[, j]$unsqueeze(1) + latvars[[j]]
    dmah2 <- dmah2 + diff2 / den
  }
  1 / (lambda + torch_matmul(torch_exp(-dmah2 / (rho^2)), invVt))
}

## ---- held-out tips: parent state z_p, observed z_i, branch length t ----
cv <- readRDS(file.path(rundir, "cv_holdout.rds"))
r3 <- readRDS(file.path(rundir, "bill_vae_aces_16dim_v3_bayesian.rds"))
lat_cols <- paste0("latent_", 1:16); edge_of <- setNames(seq_len(nrow(r3)), r3$edge)
Zp <- t(vapply(cv$held_labels, function(l) as.numeric(as.matrix(r3$z_seqs[[edge_of[[l]]]][1, lat_cols])), numeric(16)))
Yi <- as.matrix(cv$obs_latent)[, lat_cols, drop = FALSE]
tvec <- cv$term_blen; n_held <- nrow(Zp)

## sigma^2 (scalar) from non-held training-edge increments: Var(incr)=s2*blen
blens <- pf_mean_edge_features(bill_pf$phlo); names(blens) <- pf_edge_names(bill_pf$phlo)
blens <- blens[bill_pf$label]; held_set <- cv$held_labels
incr <- matrix(NA_real_, nrow(r3), 16); bl <- numeric(nrow(r3))
for (i in seq_len(nrow(r3))) {
  lab <- r3$edge[i]; if (lab %in% held_set) next
  zs <- as.matrix(r3$z_seqs[[i]][, lat_cols]); incr[i, ] <- zs[nrow(zs), ] - zs[1, ]
  bl[i] <- if (!is.na(blens[lab])) blens[lab] else NA_real_
}
ok <- is.finite(bl) & bl > 1e-8 & is.finite(rowSums(incr))
## Euclidean rate: flat increments. s2_euc = E[||incr||^2 / (d * blen)]  (d=16)
s2_euc <- mean(rowSums(incr[ok, ]^2) / (16 * bl[ok]))

## Manifold rate s2_G: PER-MODEL, metric-consistent. Under Riemannian BM,
## E[d2_G] ~= s2_G * d * t, so s2_G = E[ energy(fitted path) / (d * blen) ] over the
## NON-held training edges. The model already fitted near-geodesic paths (it minimised
## manifold energy), so their z_seqs give d2_G cheaply — no per-edge re-solve. Without
## this the manifold model is scored with a flat rate that mis-calibrates its own
## metric-scaled predictive spread and is unfairly penalised (smoke test: it "lost").
train_idx <- which(!(r3$edge %in% held_set) & is.finite(bl) & bl > 1e-8)
metric_path_energy <- function(idx_chunk) {        # returns energy per edge (length k)
  paths <- lapply(idx_chunk, function(i) as.matrix(r3$z_seqs[[i]][, lat_cols]))
  npt <- nrow(paths[[1]]); k <- length(paths)
  P <- torch_tensor(array(unlist(lapply(paths, t)), dim = c(16, npt, k)),
                    dtype = torch_float(), device = dev)$permute(c(3,1,2))   # k x 16 x npt
  vec  <- P[ , , 2:npt] - P[ , , 1:(npt-1)]
  ymid <- (P[ , , 2:npt] + P[ , , 1:(npt-1)]) / 2
  G <- metric_G(ymid$permute(c(1,3,2))$reshape(c(k*(npt-1), 16)))$reshape(c(k, npt-1, 16))$permute(c(1,3,2))
  as.numeric(((vec * G * vec)$sum(dim = 2)$sum(dim = 2) * (npt-1))$cpu())     # E = Nseg * sum vec'Gvec
}
e_tr <- numeric(0)
for (st in seq(1, length(train_idx), by = 500)) {
  en <- min(st + 499, length(train_idx))
  e_tr <- c(e_tr, metric_path_energy(train_idx[st:en]))
}
s2_G <- mean(e_tr / (16 * bl[train_idx]))
cat(sprintf("rho=%.4f lambda=%g  sigma^2: euc=%.4g  manifold(G)=%.4g  n_held=%d\n",
            rho, lambda, s2_euc, s2_G, n_held))

###############################################################################
## Pinned geodesic: min over curvature (a,b) of path energy from z_p to z_i.
## y(s) = a s^3 + b s^2 + (len - a - b) s + z_p ;  y(0)=z_p, y(1)=z_i for any a,b.
###############################################################################
zp_t  <- torch_tensor(Zp, dtype = torch_float(), device = dev)      # n x 16
len_t <- torch_tensor(Yi - Zp, dtype = torch_float(), device = dev) # n x 16
s_grid <- torch_linspace(0, 1, NSEG + 1, device = dev)$view(c(1,1,-1))  # 1x1x(NSEG+1)
a <- torch_zeros(c(n_held, 16), device = dev, requires_grad = TRUE)
b <- torch_zeros(c(n_held, 16), device = dev, requires_grad = TRUE)

path_energy <- function() {                        # returns per-tip energy (n,)
  A <- a$unsqueeze(3); B <- b$unsqueeze(3); L <- len_t$unsqueeze(3); Z0 <- zp_t$unsqueeze(3)
  y <- A * s_grid^3 + B * s_grid^2 + (L - A - B) * s_grid + Z0      # n x 16 x (NSEG+1)
  vec <- y[ , , 2:(NSEG+1)] - y[ , , 1:NSEG]                        # n x 16 x NSEG
  ymid <- (y[ , , 2:(NSEG+1)] + y[ , , 1:NSEG]) / 2
  yflat <- ymid$permute(c(1,3,2))$reshape(c(n_held * NSEG, 16))
  G <- metric_G(yflat)$reshape(c(n_held, NSEG, 16))$permute(c(1,3,2))  # n x 16 x NSEG
  e_seg <- (vec * G * vec)$sum(dim = 2)            # n x NSEG  (sum over dims)
  e_seg$sum(dim = 2) * NSEG                         # n : d2_G = N_seg * sum_seg vec'Gvec
}

opt <- optim_adam(list(a, b), lr = 0.05)
for (it in seq_len(GEO_ITERS)) {
  opt$zero_grad(); E <- path_energy()$sum(); E$backward(); opt$step()
  if (it %% 100 == 0) cat(sprintf("  geo iter %d: mean d2_G=%.4f\n", it, as.numeric((path_energy()$mean())$cpu())))
}
d2_geo <- as.numeric(path_energy()$detach()$cpu())            # n : manifold geodesic^2
d2_euc <- rowSums((Yi - Zp)^2)                                # n : straight-line^2

## logdet G at the observed tip (volume term)
logdetG <- as.numeric((torch_log(metric_G(torch_tensor(Yi, dtype = torch_float(), device = dev)))$sum(dim = 2))$cpu())

## ---- heat-kernel log-densities (higher = better), d=16, PER-MODEL sigma^2 ----
logp_m <- -d2_geo / (2 * s2_G   * tvec) - 8 * log(2 * pi * s2_G   * tvec) + 0.5 * logdetG
logp_e <- -d2_euc / (2 * s2_euc * tvec) - 8 * log(2 * pi * s2_euc * tvec)   # log det I = 0
dlog   <- logp_m - logp_e                                     # >0 => manifold better

per_tip <- tibble(tag = TAG, label = cv$held_labels, term_blen = tvec,
                  d2_geo = d2_geo, d2_euc = d2_euc, logdetG = logdetG,
                  logp_manifold = logp_m, logp_euclid = logp_e, dlog = dlog)
write_csv(per_tip, file.path(rundir, "cv_heatkernel_per_tip.csv"))

summ <- tibble(tag = TAG, rho = rho, lambda = lambda, sigma2_euc = s2_euc, sigma2_G = s2_G, n_held = n_held,
               mean_logp_manifold = mean(logp_m), mean_logp_euclid = mean(logp_e),
               mean_dlog = mean(dlog), win_manifold = mean(dlog > 0),
               mean_d2_geo = mean(d2_geo), mean_d2_euc = mean(d2_euc))
write_csv(summ, file.path(rundir, "cv_heatkernel.csv"))
cat("\n=== HEAT-KERNEL CV RESULT ===\n")
cat(sprintf("mean log p : manifold %.3f  vs euclid %.3f   (Δ=%.3f; manifold wins %.0f%% of tips)\n",
            summ$mean_logp_manifold, summ$mean_logp_euclid, summ$mean_dlog, 100*summ$win_manifold))
cat(sprintf("mean d^2   : geodesic %.4f  vs straight %.4f\n", summ$mean_d2_geo, summ$mean_d2_euc))
cat(sprintf("wrote %s/cv_heatkernel.csv + cv_heatkernel_per_tip.csv\n", rundir))
