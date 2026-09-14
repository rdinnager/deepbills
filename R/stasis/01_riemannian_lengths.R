###############################################################################
## 01_riemannian_lengths.R  —  keystone data build for the stasis test (CPU)
##
## Reimplements the v3 diagonal metric G(z) in pure R (no torch/CUDA) and
## computes, per edge, three "rulers" of evolutionary change for R1/R2/R3:
##   L_str = ||z_end - z_start||           (naive straight-line displacement)
##   L_euc = sum ||dz_k||                  (Euclidean arc length of fitted path)
##   L_G   = sum sqrt(dz_k^T G(z_mid) dz_k) (Riemannian geodesic length)
## plus branch length T_e and node age (height above root).
##
## Metric (faithful to get_metric_tensor / mahalanobis_squared in
## R/.VAE_evo_model_v3_bayesian.R):
##   dMah^2(z, mu_i) = sum_j (z_j - mu_ij)^2 / (z_j^2 + var_ij + latvar_j)
##   mh_i(z)         = exp(-dMah^2_i / rho^2)
##   G_j(z)          = 1 / (lambda + sum_i mh_i(z) * (1/var_ij))
## Output: output/stasis/edge_lengths.csv
###############################################################################
## Two-repo split 2026-09-14: run from the inner repo root (deepbills/deepbills/, which holds _targets.R).
stopifnot(file.exists("_targets.R"))
suppressMessages({library(tidyverse)})
dir.create("output/stasis", showWarnings = FALSE, recursive = TRUE)

lambda <- 1e-2
rho    <- as.numeric(readRDS("data/v3_bayesian/estimated_rho_16dim_rho_schedule_v3.rds"))
cat("rho =", rho, "  lambda =", lambda, "\n")

## --- centroids / vars : replicate training pivot exactly, capture dim order ---
bill_z <- read_csv("data/bills_vae_latent_codes_16dim.csv", show_col_types = FALSE)
cent_df <- bill_z %>% select(Species, latent_dim, latent_mean) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_mean)
var_df  <- bill_z %>% select(Species, latent_dim, latent_var) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_var)
cent_order <- setdiff(names(cent_df), "Species")          # model's dim order
cat("centroid dim order:", paste(cent_order, collapse=", "), "\n")
C  <- as.matrix(cent_df[, cent_order])                    # 2021 x 16
V  <- as.matrix(var_df[,  cent_order])                    # 2021 x 16
invV <- 1 / V

## --- latvars : replicate lat_vars <- apply(bill_tips, 2, sd) ---
bill_pf <- readRDS("data/bill_pf_16dim.rds")
bill_tips <- bill_pf %>% filter(is_tip) %>% select(starts_with("latent_")) %>% as.matrix()
latvars <- apply(bill_tips, 2, sd)                        # named latent_1..16
latvars <- latvars[cent_order]                            # align to model dim order
cat("latvars:", paste(round(latvars,3), collapse=", "), "\n")

## --- metric: G_j at a matrix of points Z (m x 16), columns in cent_order ---
metric_G <- function(Z) {
  m <- nrow(Z)
  dmah2 <- matrix(0, m, nrow(C))                          # m x n_cent
  for (j in seq_len(ncol(Z))) {
    diff2 <- outer(Z[, j], C[, j], "-")^2                 # m x n_cent
    den   <- outer(Z[, j]^2, V[, j], "+") + latvars[j]    # m x n_cent
    dmah2 <- dmah2 + diff2 / den
  }
  mh <- exp(-dmah2 / rho^2)                               # m x n_cent
  denomG <- lambda + mh %*% invV                          # m x 16
  1 / denomG                                              # G_j(z), m x 16
}

## --- per-path length triple from a 50 x 16 matrix of latent points ---
path_lengths <- function(Zpath) {          # Zpath columns already in cent_order
  n <- nrow(Zpath)
  dz  <- Zpath[2:n, , drop=FALSE] - Zpath[1:(n-1), , drop=FALSE]   # (n-1) x 16
  zmid <- (Zpath[2:n, , drop=FALSE] + Zpath[1:(n-1), , drop=FALSE]) / 2
  G   <- metric_G(zmid)                                            # (n-1) x 16
  L_euc <- sum(sqrt(rowSums(dz^2)))
  L_G   <- sum(sqrt(rowSums(G * dz^2)))
  L_str <- sqrt(sum((Zpath[n, ] - Zpath[1, ])^2))
  c(L_str = L_str, L_euc = L_euc, L_G = L_G)
}

## --- times / tree structure (edge_trajs is row-aligned with r3) ---
et <- readRDS("data/bill_edge_trajs_16dim.rds")
times <- tibble(edge = et$end,
                start_time = et$start_time, end_time = et$end_time,
                T_e = et$end_time - et$start_time,
                node_age = et$end_time)          # height above root of child node

## R1 straight paths from edge_trajs BM predictions (pred_start -> pred_end)
r1_start <- as.matrix(et[, paste0(cent_order, "_pred_start")])
r1_end   <- as.matrix(et[, paste0(cent_order, "_pred_end")])
colnames(r1_start) <- cent_order; colnames(r1_end) <- cent_order

###############################################################################
## Compute lengths for R3 (curved), R2 (straight between R3 endpoints), R1
###############################################################################
r3 <- readRDS("data/v3_bayesian/bill_vae_aces_16dim_v3_bayesian.rds")
n_edge <- nrow(r3)
cat("computing lengths for", n_edge, "edges x 3 models ...\n")

R3 <- matrix(NA, n_edge, 3); R2 <- matrix(NA, n_edge, 3); R1 <- matrix(NA, n_edge, 3)
t0 <- Sys.time()
for (i in seq_len(n_edge)) {
  Zc <- as.matrix(r3$z_seqs[[i]][, cent_order])       # R3 curved path (50 x 16)
  R3[i, ] <- path_lengths(Zc)
  ## R2: straight path between the SAME endpoints (shared nodes), 50 pts
  z0 <- Zc[1, ]; z1 <- Zc[nrow(Zc), ]
  s  <- seq(0, 1, length.out = nrow(Zc))
  Zs <- outer(1 - s, z0) + outer(s, z1); colnames(Zs) <- cent_order
  R2[i, ] <- path_lengths(Zs)
  ## R1: straight path between R1's OWN BM endpoints, 50 pts
  z0r <- r1_start[i, ]; z1r <- r1_end[i, ]
  Zr  <- outer(1 - s, z0r) + outer(s, z1r); colnames(Zr) <- cent_order
  R1[i, ] <- path_lengths(Zr)
  if (i %% 500 == 0) cat("  ", i, "/", n_edge, " (", round(difftime(Sys.time(), t0, units="secs"),1), "s)\n")
}
cat("done in", round(difftime(Sys.time(), t0, units="mins"),2), "min\n")

mk <- function(M, model) tibble(edge = r3$edge, is_tip = r3$is_tip, model = model,
                                L_str = M[,1], L_euc = M[,2], L_G = M[,3])
edge_lengths <- bind_rows(mk(R1,"R1"), mk(R2,"R2"), mk(R3,"R3")) %>%
  left_join(times, by = "edge") %>%
  mutate(rate_str = L_str / T_e, rate_euc = L_euc / T_e, rate_G = L_G / T_e,
         Phi = L_G / L_str, phi_bend = L_euc / L_str, phi_metric = L_G / L_euc)

write_csv(edge_lengths, "output/stasis/edge_lengths.csv")
cat("wrote output/stasis/edge_lengths.csv  (", nrow(edge_lengths), "rows )\n")

## quick sanity print
edge_lengths %>% group_by(model) %>%
  summarise(across(c(L_str, L_euc, L_G, Phi, phi_bend, phi_metric),
                   ~median(.x, na.rm=TRUE))) %>% print()
