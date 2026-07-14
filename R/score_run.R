###############################################################################
## score_run.R — per-sweep-run scoring for the deepbills v3 sensitivity sweep.
##
## Computes the two quantities that ACTUALLY depend on the retrained evo model:
##   (#3) geodesic stasis change-vs-time slope  (port of stasis 01 + 04)
##   (#4) decoded-tip fidelity  (decode predicted tip latent -> 64-d DeepSDF code,
##        per-dim R2/RMSE vs the true code)
## Reads the run's OWN rho (estimated_rho...rds) and lambda (convergence.csv) so
## the metric is faithful to that run. Writes data/<TAG>/score.csv (one row).
##
## (#1 0.219 diet tie and #2 R5 0.164 do NOT read the aces — invariant across the
##  sweep by construction — so they are anchors, not per-run scores.)
##
## Run per tag:  SCORE_TAG=<out_tag> Rscript R/score_run.R
###############################################################################
setwd(Sys.getenv("DEEPBILLS_DIR", unset = "~/scratch/deepbills"))
suppressMessages({library(tidyverse); library(ape); library(phyf); library(torch)})
select <- dplyr::select

TAG <- Sys.getenv("SCORE_TAG")
stopifnot(nzchar(TAG))
rundir <- file.path("data", TAG)
cat("=== scoring run:", TAG, "===\n")

## ---- run-specific hyperparameters ----
conv <- read.csv(file.path(rundir, "convergence.csv"))
lambda <- as.numeric(conv$lambda[1])
rho    <- as.numeric(readRDS(file.path(rundir, "estimated_rho_16dim_rho_schedule_v3.rds")))
cat(sprintf("rho=%.5f  lambda=%g  (epochs_run=%s stopped_early=%s final_loss=%.4f)\n",
            rho, lambda, conv$epochs_run[1], conv$stopped_early[1], conv$final_loss[1]))

## ---- shared: centroids / vars / latvars (replicate training pivots) ----
bill_z <- read_csv("data/bills_vae_latent_codes_16dim.csv", show_col_types = FALSE)
cent_df <- bill_z %>% select(Species, latent_dim, latent_mean) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_mean)
var_df  <- bill_z %>% select(Species, latent_dim, latent_var) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_var)
cent_order <- setdiff(names(cent_df), "Species")
C  <- as.matrix(cent_df[, cent_order]); V <- as.matrix(var_df[, cent_order]); invV <- 1/V

bill_pf <- readRDS("data/bill_pf_16dim.rds")
bill_tips_mat <- bill_pf %>% filter(is_tip) %>% select(starts_with("latent_")) %>% as.matrix()
latvars <- apply(bill_tips_mat, 2, sd)[cent_order]

metric_G <- function(Z) {          # Z: m x 16 in cent_order
  m <- nrow(Z); dmah2 <- matrix(0, m, nrow(C))
  for (j in seq_len(ncol(Z))) {
    diff2 <- outer(Z[, j], C[, j], "-")^2
    den   <- outer(Z[, j]^2, V[, j], "+") + latvars[j]
    dmah2 <- dmah2 + diff2 / den
  }
  mh <- exp(-dmah2 / rho^2)
  1 / (lambda + mh %*% invV)
}
path_LG_Lstr <- function(Zpath) {
  n <- nrow(Zpath)
  dz   <- Zpath[2:n, , drop=FALSE] - Zpath[1:(n-1), , drop=FALSE]
  zmid <- (Zpath[2:n, , drop=FALSE] + Zpath[1:(n-1), , drop=FALSE]) / 2
  G <- metric_G(zmid)
  c(L_str = sqrt(sum((Zpath[n,] - Zpath[1,])^2)), L_G = sum(sqrt(rowSums(G * dz^2))))
}

###############################################################################
## (#3) geodesic stasis change-vs-time slope
###############################################################################
geo_slope <- NA_real_; eucpath_slope <- NA_real_; eucobs_slope <- NA_real_
n_pairs <- NA_integer_
res3 <- tryCatch({
  r3 <- readRDS(file.path(rundir, "bill_vae_aces_16dim_v3_bayesian.rds"))
  n_edge <- nrow(r3)
  LG <- numeric(n_edge); LST <- numeric(n_edge)
  for (i in seq_len(n_edge)) {
    Zc <- as.matrix(r3$z_seqs[[i]][, cent_order])
    pl <- path_LG_Lstr(Zc); LST[i] <- pl["L_str"]; LG[i] <- pl["L_G"]
  }
  el <- tibble(edge = r3$edge, L_G = LG, L_str = LST)
  Lg   <- setNames(el$L_G,   el$edge)
  Lstr <- setNames(el$L_str, el$edge)

  tree <- pf_as_phylo(bill_pf); Ntip <- length(tree$tip.label)
  tips <- bill_pf %>% filter(is_tip)
  Y <- as.matrix(tips %>% select(starts_with("latent_"))); rownames(Y) <- tips$label
  Y <- Y[tree$tip.label, ]
  child <- tree$edge[, 2]; child_lab <- character(length(child))
  tip_edge <- child <= Ntip
  child_lab[tip_edge]  <- tree$tip.label[child[tip_edge]]
  child_lab[!tip_edge] <- tree$node.label[child[!tip_edge] - Ntip]
  te_Lg <- Lg[child_lab]; te_Lstr <- Lstr[child_lab]

  tree_Lg <- tree; tree_Lg$edge.length   <- ifelse(is.na(te_Lg),  0, te_Lg)
  tree_Lstr <- tree; tree_Lstr$edge.length <- ifelse(is.na(te_Lstr),0, te_Lstr)
  D_time <- cophenetic(tree); D_geo <- cophenetic(tree_Lg); D_pstr <- cophenetic(tree_Lstr)

  set.seed(1)
  ut <- which(upper.tri(D_time), arr.ind = TRUE); tvals <- D_time[ut]
  dec <- cut(tvals, quantile(tvals, seq(0,1,.1), na.rm=TRUE), include.lowest=TRUE)
  idx <- unlist(lapply(split(seq_len(nrow(ut)), dec), function(g) sample(g, min(3000, length(g)))))
  ut <- ut[idx, ]; i <- ut[,1]; j <- ut[,2]
  lab_i <- rownames(D_time)[i]; lab_j <- colnames(D_time)[j]
  pw <- tibble(time = D_time[cbind(i,j)],
               euc_obs = sqrt(rowSums((Y[lab_i,,drop=FALSE] - Y[lab_j,,drop=FALSE])^2)),
               euc_path = D_pstr[cbind(i,j)], geo = D_geo[cbind(i,j)]) %>%
    filter(time > 0, euc_obs > 0, geo > 0)
  sl <- function(y,x) coef(lm(log10(y) ~ log10(x)))[2]
  geo_slope     <<- as.numeric(sl(pw$geo,      pw$time))
  eucpath_slope <<- as.numeric(sl(pw$euc_path, pw$time))
  eucobs_slope  <<- as.numeric(sl(pw$euc_obs,  pw$time))
  n_pairs       <<- nrow(pw)
  cat(sprintf("#3 geodesic change-slope = %.3f  (euc_path %.3f, euc_obs %.3f; n=%d)\n",
              geo_slope, eucpath_slope, eucobs_slope, n_pairs))
  TRUE
}, error = function(e) { cat("#3 FAILED:", conditionMessage(e), "\n"); FALSE })

###############################################################################
## (#4) decoded-tip fidelity: decode predicted tip latent -> 64-d code vs true
###############################################################################
tipR2_corr_mean <- NA_real_; tipR2_corr_med <- NA_real_
tipR2_proper_mean <- NA_real_; tipRMSE_std_mean <- NA_real_; n_tip_matched <- NA_integer_
res4 <- tryCatch({
  active_dims <- readRDS("data/active_dims_16dim.rds")
  mns <- readRDS("data/code_dat_means_sds.rds")
  code_means <- mns$means; code_sds <- mns$sds
  truec <- read_csv("data/true_codes_64.csv", show_col_types = FALSE)
  code_cols <- paste0("latent_code_", 1:64)

  r3 <- readRDS(file.path(rundir, "bill_vae_aces_16dim_v3_bayesian.rds"))
  tip_i <- which(r3$is_tip); tip_lab <- r3$edge[tip_i]
  lat_cols <- paste0("latent_", 1:16)
  tip_last16 <- t(vapply(tip_i, function(i){
    zs <- as.matrix(r3$z_seqs[[i]][, lat_cols]); zs[nrow(zs), ]
  }, numeric(16)))                                   # n_tip x 16

  ## decode (CPU): scatter 16 active dims into 64, run frozen decoder, take out_codes
  options(torch.serialization_version = 2)
  vae <- torch_load("data/bill_vae_w_trophic_v1.to")
  dec <- vae$decoder; dec$eval()
  z64 <- matrix(0, nrow(tip_last16), 64); z64[, active_dims] <- tip_last16
  pred_std <- with_no_grad({ as.matrix(dec(torch_tensor(z64, dtype = torch_float()))$out_codes$cpu()) }) # n_tip x 64 (standardized)

  ## true codes -> standardized space (pred is already standardized)
  tmatch <- truec[match(tip_lab, truec$label), ]
  keep <- !is.na(tmatch$label)
  pred_std <- pred_std[keep, , drop=FALSE]
  true_raw <- as.matrix(tmatch[keep, code_cols])
  true_std <- sweep(sweep(true_raw, 2, code_means, "-"), 2, code_sds, "/")
  n_tip_matched <<- nrow(true_std)

  per_dim_corr2 <- sapply(1:64, function(d) suppressWarnings(cor(pred_std[,d], true_std[,d]))^2)
  ss_res <- colSums((pred_std - true_std)^2)
  ss_tot <- colSums(sweep(true_std, 2, colMeans(true_std))^2)
  per_dim_R2 <- 1 - ss_res/ss_tot
  per_dim_rmse <- sqrt(colMeans((pred_std - true_std)^2))     # standardized units

  tipR2_corr_mean  <<- mean(per_dim_corr2, na.rm=TRUE)
  tipR2_corr_med   <<- median(per_dim_corr2, na.rm=TRUE)
  tipR2_proper_mean<<- mean(per_dim_R2, na.rm=TRUE)
  tipRMSE_std_mean <<- mean(per_dim_rmse, na.rm=TRUE)
  cat(sprintf("#4 tip fidelity: R2(corr) mean=%.3f med=%.3f | R2(proper) mean=%.3f | RMSE(std) mean=%.3f | n_tip=%d\n",
              tipR2_corr_mean, tipR2_corr_med, tipR2_proper_mean, tipRMSE_std_mean, n_tip_matched))
  TRUE
}, error = function(e) { cat("#4 FAILED:", conditionMessage(e), "\n"); FALSE })

###############################################################################
## write one-row score record
###############################################################################
g <- function(nm) if (nm %in% names(conv)) conv[[nm]][1] else NA   # tolerate older convergence.csv
row <- data.frame(
  out_tag = TAG,
  seed = g("seed"), epochs_run = g("epochs_run"), stopped_early = g("stopped_early"),
  final_loss = g("final_loss"), best_loss = g("best_loss"),
  trophic_w = g("trophic_w"), tip_w = g("tip_w"),
  rho_target_div = g("rho_target_div"), lambda = g("lambda"), rho = rho,
  geo_slope = geo_slope, eucpath_slope = eucpath_slope, eucobs_slope = eucobs_slope, n_pairs = n_pairs,
  tipR2_corr_mean = tipR2_corr_mean, tipR2_corr_med = tipR2_corr_med,
  tipR2_proper_mean = tipR2_proper_mean, tipRMSE_std_mean = tipRMSE_std_mean,
  n_tip_matched = n_tip_matched
)
write.csv(row, file.path(rundir, "score.csv"), row.names = FALSE)
cat("wrote", file.path(rundir, "score.csv"), "\n")
