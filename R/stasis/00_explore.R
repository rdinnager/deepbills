## Explore data structures for the stasis test (CPU, no torch)
## Two-repo split 2026-09-14: run from the inner repo root (deepbills/deepbills/, which holds _targets.R).
stopifnot(file.exists("_targets.R"))
suppressMessages({library(tidyverse)})

cat("=== R3 extract (curved) ===\n")
r3 <- readRDS("data/v3_bayesian/bill_vae_aces_16dim_v3_bayesian.rds")
cat("class:", class(r3), " dim:", paste(dim(r3), collapse="x"), "\n")
cat("names:", paste(names(r3), collapse=", "), "\n")
cat("nrow:", nrow(r3), "\n")
cat("first z_seqs class:", class(r3$z_seqs[[1]]), " dim:", paste(dim(r3$z_seqs[[1]]), collapse="x"), "\n")
print(head(r3$z_seqs[[1]], 3))
cat("is_tip table:\n"); print(table(r3$is_tip))
cat("edge head:\n"); print(head(r3$edge))

cat("\n=== R1 preds (bill_w_preds_16dim) ===\n")
r1 <- readRDS("data/bill_w_preds_16dim.rds")
cat("class:", class(r1), " dim:", paste(dim(r1), collapse="x"), "\n")
cat("names:", paste(head(names(r1),40), collapse=", "), "\n")

cat("\n=== edge trajs (times) ===\n")
et <- readRDS("data/bill_edge_trajs_16dim.rds")
cat("class:", class(et), " dim:", paste(dim(et), collapse="x"), "\n")
cat("names:", paste(names(et), collapse=", "), "\n")
if ("time_seqs" %in% names(et)) {
  cat("first time_seqs:", paste(round(head(et$time_seqs[[1]],5),4), collapse=", "), " ... len:", length(et$time_seqs[[1]]), "\n")
}

cat("\n=== latent codes csv (centroids/vars) ===\n")
bz <- read_csv("data/bills_vae_latent_codes_16dim.csv", show_col_types = FALSE)
cat("names:", paste(names(bz), collapse=", "), "\n")
cat("n distinct latent_dim:", length(unique(bz$latent_dim)), " levels:", paste(sort(unique(bz$latent_dim)), collapse=","), "\n")
cat("n species:", length(unique(bz$Species)), "\n")

cat("\n=== active dims ===\n")
ad <- readRDS("data/active_dims_16dim.rds"); print(ad)

cat("\n=== rho ===\n")
cat("target_rho:", readRDS("data/v3_bayesian/estimated_rho_16dim_rho_schedule_v3.rds"), "\n")

cat("\n=== bill_pf ===\n")
bpf <- readRDS("data/bill_pf_16dim.rds")
cat("names:", paste(head(names(bpf),40), collapse=", "), "\n")
cat("nrow:", nrow(bpf), " is_tip sum:", sum(bpf$is_tip), "\n")
