###############################################################################
## Riemannian Manifold Evolutionary Model — v3 (Bayesian / Energy formulation)
## PARAMETERIZED for hyperparameter sweeps.
##
## All tunable hyperparameters are read from environment variables, defaulting
## to the current v3 values. Fractions like "1/64" are accepted. See the block
## right after library() calls. Everything else is identical to
## .VAE_evo_model_v3_bayesian.R.
###############################################################################

library(torch)
library(tidyverse)
library(FNN)
library(ape)
library(phyf)
library(Matrix)
library(GPUmatrix)
library(conflicted)

conflicts_prefer(dplyr::filter)

###############################################################################
## Hyperparameters from environment variables (sweep knobs)
###############################################################################
env_num <- function(name, default) {
  v <- Sys.getenv(name, unset = NA)
  if (is.na(v) || v == "") return(default)
  as.numeric(eval(parse(text = v)))   # allows "1/64" etc.
}
tip_weight      <- env_num("TIP_W",           10)
manifold_weight <- env_num("MANIFOLD_W",       1)
code_weight     <- env_num("CODE_W",        1/64)
trophic_weight  <- env_num("TROPHIC_W",     1/10)
root_weight     <- env_num("ROOT_W",       1/100)
rho_start_mult  <- env_num("RHO_START_MULT",   3)
rho_target_div  <- env_num("RHO_TARGET_DIV",   3)
lambda_val      <- env_num("LAMBDA",        1e-2)
n_epoch         <- as.integer(env_num("N_EPOCHS", 2500))   # hard cap
out_tag         <- Sys.getenv("OUT_TAG", unset = "v3_bayesian")

## Early stopping (sweep efficiency). Stop once total loss has plateaued:
## no improvement in the best loss by more than TOL (relative) over the last
## PATIENCE epochs, with a MIN_EPOCHS floor (so rho annealing completes before
## stopping is allowed) and N_EPOCHS as the hard cap.
env_flag <- function(name, default_on = TRUE) {
  v <- tolower(Sys.getenv(name, unset = if (default_on) "on" else "off"))
  v %in% c("on", "true", "1", "yes", "t")
}
early_stop <- env_flag("EARLY_STOP", TRUE)
patience   <- as.integer(env_num("PATIENCE",   200))
tol        <- env_num("TOL", 0.005)   # 0.5% relative improvement threshold

## --- Schedules are decoupled from the epoch CAP (fixed 2026-07-14) -----------
## Previously the rho anneal used prop_epoch = epoch / N_EPOCHS(cap), so a run
## that early-stopped at ~half the cap ended MID-ANNEAL, at a rho it was never
## meant to be evaluated at. Worse, total_loss is not comparable across epochs
## while rho moves: as rho falls the metric sharpens and the SAME path costs more
## energy, so the loss rises even as the fit improves (measured: tip_loss -14.6%
## while manifold_energy +179% over the "plateau"). Early stopping on that is
## stopping on a non-stationary objective.
##
## Now: rho anneals over ANNEAL_EPOCHS and is then HELD at target. Early stopping
## is only permitted once the objective is stationary (anneal done + a grace
## period), so "converged" and "best" finally mean what they say. This is the
## standard continuation/graduated-non-convexity discipline: the schedule is a
## device for reaching a good basin of the FINAL objective; you anneal to target,
## then train at fixed target to convergence.
anneal_epochs <- as.integer(env_num("ANNEAL_EPOCHS", 800))
stop_grace    <- as.integer(env_num("STOP_GRACE",    200))  # epochs at fixed rho before stopping may fire
min_epochs    <- as.integer(env_num("MIN_EPOCHS", anneal_epochs + stop_grace))
## LR schedule horizon, also independent of the cap. The one-cycle is stepped over
## LR_EPOCHS and then held at its final LR. (Previously `scheduler` was created but
## `scheduler$step()` was NEVER called, so the LR sat constant at max_lr/25 = 8e-4
## and the one-cycle was dead code.)
lr_epochs <- as.integer(env_num("LR_EPOCHS", 1500))

## --- Speed knobs (2026-07-14) ------------------------------------------------
## Measured on Vulcan: GPU utilisation was 41%/0%/0% with only 12 GB of 48 GB used
## — the run was HOST-bound, not GPU-bound. Two culprits, both fixed here:
##  FAST_BATCH: keep the tree matrices resident on the GPU instead of slicing them
##    on the host and copying them across on every batch of every epoch.
##  SAVE_LATEST_EVERY: `checkpoint_latest.to` was written to Lustre EVERY epoch.
##    The best checkpoint is saved on improvement regardless, so this is only a
##    crash-resume convenience — it does not need to run every epoch.
## FAST_BATCH is OFF by default: the idea is sound (the host-side per-batch slicing
## is real, and the GPU measurably idles) but keeping the full matrices resident makes
## the process hold ~25 GB before the first forward pass and then OOM on a 24.3 GiB
## allocation. The `gpu.matrix()` conversion of the FULL sparse matrix is doing
## something far more expensive than the per-batch slices did — not yet diagnosed.
## Do NOT turn this on without re-testing memory. See compute-optimization-plan.md.
fast_batch        <- env_flag("FAST_BATCH", FALSE)
save_latest_every <- as.integer(env_num("SAVE_LATEST_EVERY", 100))

## LR range test (Smith 2017). Ramps LR geometrically LR_MIN -> LR_MAX_TEST after the
## anneal completes; read the resulting loss-vs-lr curve to choose max_lr.
lr_range_test <- env_flag("LR_RANGE_TEST", FALSE)
lr_min        <- env_num("LR_MIN", 1e-5)
lr_max_test   <- env_num("LR_MAX_TEST", 1)

## TOPK_K: 0 = exact metric (all 2021 centroids). >0 = top-K nearest centroids only.
## An APPROXIMATION — see get_metric_tensor_topk(). Default 0 (exact) until validated
## end-to-end against the geodesic slope.
topk_k <- as.integer(env_num("TOPK_K", 0))

## Optional RNG seed (for replicate runs / reproducibility). If unset, the run
## is non-deterministic (original behaviour). Seeds both base R and torch, which
## together drive the a/b/root init and the dataloader shuffle.
seed_env <- Sys.getenv("SEED", unset = NA)
if (!is.na(seed_env) && seed_env != "") {
  seed_int <- as.integer(seed_env)
  set.seed(seed_int)
  torch::torch_manual_seed(seed_int)
} else {
  seed_int <- NA_integer_
}

cat("=== HYPERPARAMETERS ===\n")
cat(sprintf("  TIP_W=%g  MANIFOLD_W=%g  CODE_W=%g  TROPHIC_W=%g  ROOT_W=%g\n",
            tip_weight, manifold_weight, code_weight, trophic_weight, root_weight))
cat(sprintf("  RHO_START_MULT=%g  RHO_TARGET_DIV=%g  LAMBDA=%g\n",
            rho_start_mult, rho_target_div, lambda_val))
cat(sprintf("  N_EPOCHS=%d (cap)  OUT_TAG=%s\n", n_epoch, out_tag))
cat(sprintf("  EARLY_STOP=%s  PATIENCE=%d  MIN_EPOCHS=%d  TOL=%g  SEED=%s\n",
            early_stop, patience, min_epochs, tol,
            ifelse(is.na(seed_int), "none", seed_int)))

## Create output directory (avoids overwriting existing data/)
out_dir <- file.path("data", out_tag)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## get decoder only from bill vae
options(torch.serialization_version = 2)
bill_vae <- torch_load("data/bill_vae_w_trophic_v1.to")
bill_decoder <- bill_vae$decoder
bill_decoder$parameters %>% purrr::walk(function(param) param$requires_grad_(FALSE))
bill_decoder <- bill_decoder$cuda()

active_dims <- read_rds("data/active_dims_16dim.rds")

bill_pf <- read_rds("data/bill_pf_16dim.rds")
init_rates <- read_csv("data/init_rates_16dim.csv")
bill_pf <- bill_pf %>%
  left_join(init_rates %>%
              rename_with(~ paste0("rate_", .x),
                          starts_with("latent_")),
            by = c("label" = "label"))
blens <- pf_mean_edge_features(bill_pf$phlo)
names(blens) <- pf_edge_names(bill_pf$phlo)
blens <- blens[bill_pf$label]

bill_tips <- bill_pf %>%
  filter(is_tip) %>%
  select(starts_with("latent_")) %>%
  as.matrix()

lat_vars <- apply(bill_tips, 2, sd)

bill_tree_mat <- pf_as_sparse(bill_pf$phlo)
bill_tree_mat <- bill_tree_mat[ , bill_pf$label]

start_mat <- pf_as_sparse(bill_pf$phlo)
start_mat <- start_mat[ , bill_pf$label]
for(i in 1:nrow(start_mat)) {
  start_mat[i, colnames(start_mat) == rownames(start_mat)[i]] <- 0
}

bill_tree_mat_tips <- bill_pf %>%
  filter(is_tip) %>%
  pf_as_sparse()
bill_tree_mat_tips <- bill_tree_mat_tips[ , bill_pf$label]

bill_dataset <- dataset(name = "bill_ds",
                           initialize = function(end_tree_mat,
                                                 start_tree_mat,
                                                 tip_tree_mat,
                                                 tip_dat,
                                                 blens) {
                             self$end_tree_mat <- end_tree_mat
                             self$start_tree_mat <- start_tree_mat
                             self$tip_tree_mat <- gpu.matrix(tip_tree_mat, dtype = "float32", device = "cuda")@gm
                             self$tip_dat <- torch_tensor(tip_dat, device = "cuda")
                             self$blens <- torch_tensor(blens, device = "cuda")
                             ## FAST_BATCH (2026-07-14): the tree matrices used to be sliced on the
                             ## HOST and pushed to the GPU on EVERY batch of EVERY epoch
                             ## (gpu.matrix(self$end_tree_mat[i, ], ...)). That host work serialised
                             ## against GPU compute and left the card idle ~60% of the time
                             ## (measured: 41%/0%/0% utilisation, 12 GB of 48 GB used).
                             ## Now the FULL matrices live on the GPU permanently; the forward pass
                             ## multiplies once and selects the batch's rows afterwards — identical
                             ## maths (row-slice then multiply == multiply then row-slice), no host
                             ## work, no transfers. The matrices are sparse and tiny, so the full
                             ## product is free.
                             self$fast_batch <- fast_batch
                             if (self$fast_batch) {
                               if (Sys.getenv("DBG_SIZES", "off") == "on") {
                                 message(sprintf("[dbg] end_tree_mat dim = %d x %d ; start = %d x %d ; tip = %d x %d",
                                                 nrow(end_tree_mat), ncol(end_tree_mat),
                                                 nrow(start_tree_mat), ncol(start_tree_mat),
                                                 nrow(tip_tree_mat), ncol(tip_tree_mat)))
                               }
                               self$end_full   <- gpu.matrix(end_tree_mat,   dtype = "float32", device = "cuda")@gm
                               self$start_full <- gpu.matrix(start_tree_mat, dtype = "float32", device = "cuda")@gm
                               if (Sys.getenv("DBG_SIZES", "off") == "on") {
                                 message(sprintf("[dbg] resident GPU matrices built OK (end rows=%d)",
                                                 self$end_full$size(1)))
                               }
                             }
                           },
                           .getbatch = function(i) {
                              if (self$fast_batch) {
                                ## no host slicing, no H2D copy: hand back the resident matrices and
                                ## let the forward pass do the row selection on-device.
                                return(list(self$end_full, self$start_full, self$tip_tree_mat,
                                            self$tip_dat, self$blens[i],
                                            torch_tensor(i, device = "cuda")))
                              }
                              end_tree_mat <- gpu.matrix(self$end_tree_mat[i, ], dtype = "float32", device = "cuda")@gm
                              start_tree_mat <- gpu.matrix(self$start_tree_mat[i, ], dtype = "float32", device = "cuda")@gm
                              list(end_tree_mat, start_tree_mat, self$tip_tree_mat, self$tip_dat, self$blens[i], torch_tensor(i, device = "cuda"))                           },
                           .length = function() {
                             nrow(self$end_tree_mat)
                           })

write_rds(list(bill_tree_mat = bill_tree_mat,
               start_mat = start_mat,
               bill_tree_mat_tips = bill_tree_mat_tips,
               bill_tips = bill_tips,
               blens = blens), file.path(out_dir, "dataloader_dat_v3.rds"))

batch_size <- 200

bill_ds <- bill_dataset(bill_tree_mat, start_mat, bill_tree_mat_tips, bill_tips, blens)
bill_dl <- dataloader(bill_ds, batch_size, shuffle = TRUE, drop_last = FALSE)

bill_init_rates <- bill_pf %>%
  select(starts_with("rate_")) %>%
  as.matrix()

bill_z <- read_csv("data/bills_vae_latent_codes_16dim.csv")
centroids_df <- bill_z %>%
  select(Species, latent_dim, latent_mean) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_mean)

centroids <- centroids_df %>%
  select(-Species) %>%
  as.matrix() %>%
  torch_tensor(device = "cuda")

vars <- bill_z %>%
  select(Species, latent_dim, latent_var) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_var) %>%
  select(-Species) %>%
  as.matrix() %>%
  torch_tensor(device = "cuda")

## N_SEGS: path points per edge for the energy integral. 50 by default; 25 halves the
## work in BOTH the metric kernel and the decoder. Changing it changes the
## discretisation of the path integral, so it must be validated (see
## compute-optimization-plan.md), not just assumed.
n_segs <- as.integer(env_num("N_SEGS", 50))

###############################################################################
## Decoder pass: decode latent path points to phenotype and trophic space
###############################################################################

decode_segments <- function(y, bill_decoder, active_dims, device = "cuda") {
  y_size_1 <- y$size()[1]
  y_size_3 <- y$size()[3]
  y_array <- torch_zeros(y_size_1, 64, y_size_3, device = device)
  indexer <- torch_tensor(active_dims, device = device)$unsqueeze(1)$unsqueeze(-1)$expand(c(y_size_1, 16, y_size_3))
  y_array <- y_array$scatter(2, indexer, y)
  y_mat <- y_array$movedim(3, 2)$flatten(end_dim = 2)
  z_mat <- bill_decoder(y_mat)
  code_array <- z_mat$out_codes$unflatten(1, c(y_size_1, y_size_3))$movedim(3, 2)
  code_diff <- torch_diff(code_array) + 1e-6
  trophic_array <- z_mat$out_trophic$unflatten(1, c(y_size_1, y_size_3))$movedim(3, 2)
  trophic_diff <- torch_diff(nnf_softmax(trophic_array, dim = 2)) + 1e-6
  list(code_diff, trophic_diff)
}

###############################################################################
## Cubic path parameterization: z(s) for s in [0, 1]
## z(s) = a*s^3 + b*s^2 + (lens - a - b)*s + z_start
###############################################################################

get_segments <- function(z_starts, z_ends, segs, a, b, bill_decoder, active_dims, device = "cuda") {
  lens <- z_ends - z_starts
  d <- lens$size()[2]
  n <- segs$size()
  segs <- segs$unsqueeze(1L)$`repeat`(c(lens$size()[1], 1))
  y <- (a$unsqueeze(-1)$`repeat`(c(1, 1, n)) * (segs^3)$unsqueeze(2) + b$unsqueeze(-1)$`repeat`(c(1, 1, n)) * (segs^2)$unsqueeze(2)) +
    (lens - a - b)$unsqueeze(-1) * segs$unsqueeze(2)$`repeat`(c(1, d, 1)) +
    z_starts$unsqueeze(-1)
  z <- decode_segments(y, bill_decoder, active_dims = active_dims, device = device)
  vec <- torch_diff(y) + 1e-6
  y_mid <- y[ , , 1:(y$size()[3] - 1)] + vec / 2
  list(z[[1]], z[[2]], vec, y_mid)
}

latvars_tens <- torch_tensor(lat_vars, device = "cuda")$unsqueeze(1L)

###############################################################################
## Mahalanobis distance (unchanged from v2)
###############################################################################

mahalanobis_squared_fun <- function(z1, z2, v2, latvars_tens) {

  n_1 <- z1$size(1)
  segs <- z1$size(3)
  n_2 <- z2$size(1)
  dim <- z1$size(2)

  mahalanobis_squared <- function(z1, z2, v2, latvars_tens) {

    expanded_1 = z1$unsqueeze(2)$expand(c(n_1, n_2, dim, segs))
    expanded_2 = z2$unsqueeze(1)$unsqueeze(-1)$expand(c(n_1, n_2, dim, segs))
    expanded_3 = (v2 + latvars_tens)$unsqueeze(1)$unsqueeze(-1)$expand(c(n_1, n_2, dim, segs))

    diff <- expanded_1 - expanded_2
    torch_sum((diff^2) / ((expanded_1^2) + expanded_3), 3)
  }

  mahalanobis_squared_tr <- jit_trace(mahalanobis_squared,
                                      z1, z2, v2, latvars_tens)

  mahalanobis_squared <- function(z1, z2, v2, latvars_tens) {

    n_1 <- z1$size(1)
    segs <- z1$size(3)
    n_2 <- z2$size(1)
    dim <- z1$size(2)

    expanded_1 = z1$unsqueeze(2)$expand(c(n_1, n_2, dim, segs))
    expanded_2 = z2$unsqueeze(1)$unsqueeze(-1)$expand(c(n_1, n_2, dim, segs))
    expanded_3 = (v2 + latvars_tens)$unsqueeze(1)$unsqueeze(-1)$expand(c(n_1, n_2, dim, segs))

    diff <- expanded_1 - expanded_2
    torch_sum((diff^2) / ((expanded_1^2) + expanded_3), 3)
  }

  list(mahalanobis_squared, mahalanobis_squared_tr)

}

mahalanobis_simple <- function(z1, z2, v2, latvars_tens) {

  n_1 <- z1$size(1)
  n_2 <- z2$size(1)
  dim <- z1$size(2)

  expanded_1 = z1$unsqueeze(2)$expand(c(n_1, n_2, dim))
  expanded_2 = z2$unsqueeze(1)$expand(c(n_1, n_2, dim))

  diff <- expanded_1 - expanded_2
  torch_sum((diff^2) / ((expanded_1^2) + v2 + latvars_tens), 3)

}

###############################################################################
## Riemannian metric tensor G(z) — data-dependent, diagonal
###############################################################################

get_metric_tensor <- function(z, centroids, vars, lambda = lambda_val, rho, latvars_tens) {
  if (topk_k > 0) return(get_metric_tensor_topk(z, centroids, vars, lambda, rho, latvars_tens, topk_k))
  mh <- torch_exp(-(mahalanobis_squared[[1]](z, centroids, vars, latvars_tens) / (rho^2)))
  1 / (((1 / vars)$unsqueeze(1)$unsqueeze(-1) * mh$unsqueeze(3))$sum(dim = 2) + lambda)
}

## TOP-K metric (TOPK_K > 0). The kernel exp(-d2_M/rho^2) decays, so all but the nearest
## few hundred centroids contribute ~nothing to the sum. Restricting to the K nearest
## shrinks the saved activations from (points x 2021 x dim) to (points x K x dim), which
## is what the 79%-of-epoch BACKWARD pass has to traverse.
##
## HONEST CAVEAT: this is an approximation, and a BIASED one — dropping positive kernel
## mass under-estimates the sum, hence OVER-estimates G and the manifold energy.
## Measured on real centroids at rho=0.841 (validate_topk_metric.R):
##   K=256 -> mean rel. error in energy 1.5e-3, max 1.8e-2 (kernel mass kept: 99.83%)
## Neighbour selection is done under no_grad using the EXACT Mahalanobis distance;
## gradients flow only through the K retained terms.
get_metric_tensor_topk <- function(z, centroids, vars, lambda, rho, latvars_tens, K) {
  n1 <- z$size(1); dm <- z$size(2); S <- z$size(3)
  pts <- z$permute(c(1, 3, 2))$reshape(c(n1 * S, dm))          # (P, dim)

  idx <- with_no_grad({
    diff <- pts$unsqueeze(2) - centroids$unsqueeze(1)          # (P, C, dim)
    den  <- pts$unsqueeze(2)^2 + (vars + latvars_tens)$unsqueeze(1)
    m2   <- torch_sum(diff^2 / den, 3)                         # (P, C)
    m2$topk(K, dim = 2, largest = FALSE)[[2]]                  # (P, K)
  })

  mu_k <- centroids[idx, ]                                     # (P, K, dim)
  v_k  <- vars[idx, ]                                          # (P, K, dim)
  diff <- pts$unsqueeze(2) - mu_k
  den  <- pts$unsqueeze(2)^2 + (v_k + latvars_tens)
  mh   <- torch_exp(-(torch_sum(diff^2 / den, 3)) / (rho^2))   # (P, K)
  Ginv <- ((1 / v_k) * mh$unsqueeze(3))$sum(dim = 2) + lambda  # (P, dim)
  (1 / Ginv)$reshape(c(n1, S, dm))$permute(c(1, 3, 2))         # (n1, dim, S)
}

###############################################################################
## Energy formulation (squared Riemannian norm, no sqrt)
###############################################################################

get_manifold_energy <- function(vel, metric) {
  (vel * metric * vel)$sum(dim = 2)
}

get_manifold_dist <- function(vel, metric) {
  (vel * metric * vel)$sum(dim = 2)$sqrt()
}

###############################################################################
## Loss functions
###############################################################################

energy_loss <- function(energies) {
  torch_mean(energies)
}

tip_loss <- function(tip_data, tip_recon) {
  tip_dists <- torch_square(tip_data - tip_recon)
  torch_mean(tip_dists)
}

root_loss <- function(root_values) {
  torch_sum(root_values^2)
}

###############################################################################
## Model module
###############################################################################

mani_evo_mod <- nn_module("mani_evo_bayesian",
                          initialize = function(n_rates, n_dim, centroids, vars,
                                                n_segs = 100, init_rates = NULL,
                                                bill_decoder, active_dims, latvars,
                                                device = "cuda") {

                            self$device <- device
                            self$n_rates <- n_rates
                            self$n_dim <- n_dim
                            self$n_segs <- n_segs
                            self$segs <- torch_arange(0, 1, 1 / n_segs, device = device)
                            self$centroids <- centroids$to(device = device)
                            self$vars <- vars$to(device = device)

                            self$bill_decoder <- bill_decoder
                            self$active_dims <- active_dims
                            self$latvars_tens <- torch_tensor(latvars, device = device)

                            centroid_dists <- mahalanobis_simple(centroids, centroids, vars, self$latvars_tens)$sqrt()$cpu() %>%
                              as.matrix()
                            diag(centroid_dists) <- 999999999999999
                            mins <- apply(centroid_dists, 1, min)
                            rho <- max(mins)
                            target_rho <- rho / rho_target_div
                            start_rho <- rho * rho_start_mult
                            self$target_rho <- torch_tensor(target_rho, device = device)
                            self$start_rho <- torch_tensor(start_rho, device = device)
                            self$rho <- torch_tensor(start_rho, device = device)

                            self$a <- nn_parameter(torch_randn(n_rates, n_dim) * 0.001)
                            self$b <- nn_parameter(torch_randn(n_rates, n_dim) * 0.001)
                            if(!is.null(init_rates)) {
                              self$rates <- nn_parameter(torch_tensor(init_rates))
                            } else {
                              self$rates <- nn_parameter(torch_randn(n_rates, n_dim) * 0.01)
                            }
                            self$root_values <- nn_parameter(torch_randn(n_dim) * 0.001)

                          },

                          forward = function(x) {

                            dbg <- Sys.getenv("DBG_SIZES", "off") == "on"
                            if (dbg) message(sprintf("[dbg] forward entered: x1 %s  rates %s  blens n=%d",
                                                     paste(dim(x[[1]]), collapse = "x"),
                                                     paste(dim(self$rates), collapse = "x"),
                                                     x[[5]]$size(1)))

                            z_ends <- call_torch_function("torch__sparse_mm", x[[1]], self$rates, quiet = TRUE) +
                              self$root_values$unsqueeze(1)
                            if (dbg) message("[dbg] sparse_mm #1 (ends) OK")
                            z_starts <- call_torch_function("torch__sparse_mm", x[[2]], self$rates, quiet = TRUE) +
                              self$root_values$unsqueeze(1)
                            z_tips <- call_torch_function("torch__sparse_mm", x[[3]], self$rates, quiet = TRUE) +
                              self$root_values$unsqueeze(1)

                            ## FAST_BATCH: when the dataset hands back the FULL resident matrices
                            ## (instead of host-sliced per-batch ones), z_ends/z_starts come out with
                            ## one row per EDGE IN THE TREE rather than one row per edge in the batch,
                            ## so select this batch's rows now.
                            ##   (row-slice ∘ multiply  ==  multiply ∘ row-slice)
                            ## Detected from the tensors themselves rather than from a global flag:
                            ## nn_module methods do NOT see script globals, and a silently-false flag
                            ## here means get_segments runs on all 4,040 edges instead of the batch's
                            ## 200 — a 20x blow-up that OOMs a 44 GB card (tried to allocate 24.3 GiB).
                            n_batch <- x[[5]]$size(1)
                            if (Sys.getenv("DBG_SIZES", "off") == "on") {
                              ## message() -> stderr, unbuffered: survives an abort. cat() to stdout
                              ## can be swallowed by R's buffer when the process dies (which is why
                              ## the first diagnostic run printed nothing at all).
                              message(sprintf("[dbg] z_ends rows=%d  z_starts rows=%d  n_batch(blens)=%d  idx len=%d",
                                              z_ends$size(1), z_starts$size(1), n_batch, x[[6]]$size(1)))
                            }
                            if (z_ends$size(1) != n_batch) {
                              idx <- x[[6]]$to(dtype = torch_long())
                              z_ends   <- torch_index_select(z_ends,   1, idx)
                              z_starts <- torch_index_select(z_starts, 1, idx)
                            }

                            zs <- get_segments(z_starts, z_ends, self$segs, self$a[x[[6]]], self$b[x[[6]]],
                                               self$bill_decoder, self$active_dims, device = self$device)

                            met <- get_metric_tensor(zs[[4]], centroids, vars, rho = self$rho, latvars_tens = self$latvars_tens)

                            blens_expanded <- x[[5]]$unsqueeze(-1)

                            manifold_energies <- get_manifold_energy(zs[[3]], met) / blens_expanded
                            code_energies <- get_manifold_energy(zs[[1]], 1.0) / blens_expanded
                            trophic_energies <- get_manifold_energy(zs[[2]], 1.0) / blens_expanded

                            manifold_energy_loss <- energy_loss(manifold_energies)
                            code_energy_loss <- energy_loss(code_energies)
                            trophic_energy_loss <- energy_loss(trophic_energies)

                            tip_loss_val <- tip_loss(x[[4]], z_tips)
                            rootval_loss <- root_loss(self$root_values)

                            list(manifold_energy_loss, code_energy_loss, trophic_energy_loss,
                                 0, tip_loss_val, rootval_loss,
                                 zs, manifold_energies, code_energies, trophic_energies)

                          })

mod <- mani_evo_mod(n_rates = nrow(bill_init_rates), n_dim = 16,
                    centroids, vars,
                    n_segs = n_segs, init_rate = bill_init_rates,
                    bill_decoder = bill_decoder, active_dims = active_dims,
                    latvars = lat_vars,
                    device = "cuda")
mod <- mod$cuda()

write_rds(as.numeric(mod$target_rho$cpu()), file.path(out_dir, "estimated_rho_16dim_rho_schedule_v3.rds"))

## Initialize mahalanobis function (needs a test forward pass for JIT)
## NB: this warm-up does its OWN sparse_mm and does NOT go through mod$forward(), so
## it needs the same row-selection. Under FAST_BATCH the dataset hands back the FULL
## tree matrices, so without this the JIT trace is built over all 4,040 edges and the
## Mahalanobis expand allocates 4040 x 2021 x 16 x 50 floats = 24.33 GiB -> OOM.
test <- dataloader_next(dataloader_make_iter(bill_dl))
init_rates_test <- torch_tensor(bill_init_rates, device = "cuda")
z_ends_test <- call_torch_function("torch__sparse_mm", test[[1]], init_rates_test, quiet = TRUE)
z_starts_test <- call_torch_function("torch__sparse_mm", test[[2]], init_rates_test, quiet = TRUE)
if (z_ends_test$size(1) != test[[5]]$size(1)) {
  idx_test <- test[[6]]$to(dtype = torch_long())
  z_ends_test   <- torch_index_select(z_ends_test,   1, idx_test)
  z_starts_test <- torch_index_select(z_starts_test, 1, idx_test)
}
a_test <- torch_randn(z_starts_test$size()[1], z_starts_test$size()[2], device = "cuda")
b_test <- torch_randn(z_starts_test$size()[1], z_starts_test$size()[2], device = "cuda")
segs_test <- torch_arange(0, 1, 1 / n_segs, device = "cuda")
zs_test <- get_segments(z_starts_test, z_ends_test, segs_test, a_test, b_test,
                        bill_decoder, active_dims, device = "cuda")
mahalanobis_squared <- mahalanobis_squared_fun(zs_test[[4]], centroids, vars, latvars_tens)

###############################################################################
## Rho annealing schedule (posterior tempering)
###############################################################################

cosine_schedule <- function(t, start=0, end=1, tau=1, clip_min=1e-9) {
  v_start <- cos(start * pi / 2) ^ (2 * tau)
  v_end <- cos(end * pi / 2) ^ (2 * tau)
  output <- cos((t * (end - start) + start) * pi / 2) ^ (2 * tau)
  output <- (v_end - output) / (v_end - v_start)
  output[output < clip_min] <- clip_min
  output[output > 1.0] <- 1.0
  output
}

###############################################################################
## Training loop
###############################################################################

## MAX_LR: peak LR of the one-cycle. Was hardcoded at 0.02 — a value that had NEVER been
## exercised (scheduler$step() was never called, so the LR sat at max_lr/25 = 8e-4 for
## the whole run). The first time the one-cycle actually ran, 0.02 spiked the loss to 637.
## Set this from the LR range test (LR_RANGE_TEST=on), don't trust the old default.
lr <- env_num("MAX_LR", 0.02)
save_every <- 50

optim1 <- optim_adam(mod$parameters, lr = lr)
## One-cycle over LR_EPOCHS (NOT the cap), and actually stepped — see the loop.
scheduler <- lr_one_cycle(optim1, max_lr = lr,
                          epochs = lr_epochs, steps_per_epoch = 1,
                          cycle_momentum = FALSE)

optim1$zero_grad()

## Initialize loss history CSV
loss_csv <- file.path(out_dir, "loss_history.csv")
cat("epoch,total_loss,tip_loss,manifold_energy,code_energy,trophic_energy,root_loss,rho,lr,epoch_secs\n",
    file = loss_csv)

checkpoint_every <- 100

test_dat <- list()
train_start <- Sys.time()

## Early-stopping state
best_loss  <- Inf
best_epoch <- 0L
stopped_early <- FALSE
best_path  <- file.path(out_dir, "checkpoint_best.to")

for(epoch in 1:n_epoch) {

  epoch_start <- Sys.time()
  optim1$zero_grad()
  total_loss <- 0
  total_recon_loss <- 0
  total_manifold_loss <- 0
  total_code_loss <- 0
  total_trophic_loss <- 0
  total_root_loss <- 0

  i <- 0
  prof <- Sys.getenv("PROFILE", "off") == "on"
  if (prof) { p_data <- p_fwd <- p_bwd <- p_sync <- 0 }
  if (prof) t_mark <- Sys.time()
  coro::loop(for (b in bill_dl) {
      i <- i + 1
      if (prof) { cuda_synchronize(); p_data <- p_data + as.numeric(Sys.time() - t_mark, units="secs"); t_mark <- Sys.time() }

      res <- mod(b)
      if (prof) { cuda_synchronize(); p_fwd <- p_fwd + as.numeric(Sys.time() - t_mark, units="secs"); t_mark <- Sys.time() }

      if(epoch %% save_every == 0) {
        with_no_grad({
          test_dat[[i]] <- purrr::map(purrr::list_flatten(res[-4]), ~.x$cpu()$detach())
        })
      }

      loss <- (manifold_weight * res[[1]] +
                 code_weight * res[[2]] +
                 trophic_weight * res[[3]] +
                 tip_weight * res[[5]] +
                 root_weight * res[[6]]) * 100

      loss$backward()
      if (prof) { cuda_synchronize(); p_bwd <- p_bwd + as.numeric(Sys.time() - t_mark, units="secs"); t_mark <- Sys.time() }

      ## NOTE: each of these six $cpu() pulls is a BLOCKING device->host sync — six
      ## per batch, ~21 batches per epoch = ~126 pipeline stalls per epoch, purely for
      ## per-epoch bookkeeping. Profiled below as p_sync.
      total_loss <- total_loss + as.numeric(loss$cpu())
      total_recon_loss <- total_recon_loss + as.numeric(res[[5]]$cpu())
      total_manifold_loss <- total_manifold_loss + as.numeric(res[[1]]$cpu())
      total_code_loss <- total_code_loss + as.numeric(res[[2]]$cpu())
      total_trophic_loss <- total_trophic_loss + as.numeric(res[[3]]$cpu())
      total_root_loss <- total_root_loss + as.numeric(res[[6]]$cpu())
      if (prof) { cuda_synchronize(); p_sync <- p_sync + as.numeric(Sys.time() - t_mark, units="secs"); t_mark <- Sys.time() }

  })

  if (prof) {
    tot <- p_data + p_fwd + p_bwd + p_sync
    message(sprintf("[prof] epoch %d  batches=%d  data %.2fs (%.0f%%)  fwd %.2fs (%.0f%%)  bwd %.2fs (%.0f%%)  d2h-sync %.2fs (%.0f%%)  total %.2fs",
                    epoch, i, p_data, 100*p_data/tot, p_fwd, 100*p_fwd/tot,
                    p_bwd, 100*p_bwd/tot, p_sync, 100*p_sync/tot, tot))
  }

  current_rho <- as.numeric(mod$rho$cpu())
  current_lr <- optim1$param_groups[[1]]$lr
  epoch_secs <- as.numeric(Sys.time() - epoch_start, units = "secs")

  cat("Epoch: ", epoch,
      "    loss: ", as.numeric(total_loss),
      "    tip recon loss: ", as.numeric(total_recon_loss),
      "    manifold energy: ", as.numeric(total_manifold_loss),
      "    code energy: ", as.numeric(total_code_loss),
      "    trophic energy: ", as.numeric(total_trophic_loss),
      "    root loss: ", as.numeric(total_root_loss),
      "    rho: ", current_rho,
      "    secs: ", round(epoch_secs, 3),
      "\n")

  ## Append to loss history CSV
  cat(paste(epoch, total_loss, total_recon_loss, total_manifold_loss,
            total_code_loss, total_trophic_loss, total_root_loss,
            current_rho, current_lr, epoch_secs, sep = ","), "\n",
      file = loss_csv, append = TRUE)

  ## Save checkpoint
  if(epoch %% checkpoint_every == 0) {
    torch_save(mod, file.path(out_dir, paste0("checkpoint_epoch_", epoch, ".to")))
  }
  if(save_latest_every > 0 && epoch %% save_latest_every == 0) {
    torch_save(mod, file.path(out_dir, "checkpoint_latest.to"))
  }

  if(epoch %% save_every == 0) {

    test_zs <- torch_cat(map(test_dat, 9)) %>%
      as.array()
    png(file.path(out_dir, paste0("paths_epoch_", epoch, ".png")),
        width = 800, height = 800)
    plot(cbind(as.vector(test_zs[ , 1, ]), as.vector(test_zs[, 2, ])), type = "n",
         xlab = "Latent dim 1", ylab = "Latent dim 2",
         main = paste("Epoch", epoch))
    points(as.matrix(centroids$cpu())[, 1:2], cex = 0.2, col = "green")
    for(j in seq_len(dim(test_zs)[1])) {
      points(cbind(as.vector(test_zs[j, 1, ]), as.vector(test_zs[j, 2, ])), type = "l", col = "red")
    }
    points(as.matrix(centroids$cpu())[, 1:2], cex = 0.2, col = "green")
    dev.off()

  }

  optim1$step()

  ## Step the LR schedule (one-cycle over LR_EPOCHS, then hold at final LR).
  ## LR_RANGE_TEST (Smith 2017): instead of the one-cycle, ramp the LR geometrically
  ## from LR_MIN to LR_MAX once the rho anneal is done (so the objective is fixed and
  ## loss-vs-LR is meaningful). Plot loss against lr afterwards and set max_lr below
  ## the divergence point. Needed because the one-cycle's nominal max_lr = 0.02 was
  ## NEVER exercised (scheduler$step() was never called) and, when finally run, spiked
  ## the loss to 637.
  if (lr_range_test) {
    if (epoch >= anneal_epochs) {
      frac <- (epoch - anneal_epochs) / max(1, (n_epoch - anneal_epochs))
      new_lr <- lr_min * (lr_max_test / lr_min)^frac
      for (g in seq_along(optim1$param_groups)) optim1$param_groups[[g]]$lr <- new_lr
    }
  } else if (epoch <= lr_epochs) {
    scheduler$step()
  }

  ## Update rho (posterior tempering schedule).
  ## prop_epoch runs over ANNEAL_EPOCHS, not the cap, and is clamped at 1 — so the
  ## anneal COMPLETES and rho is then held exactly at target for the rest of training.
  prop_epoch <- min(epoch / anneal_epochs, 1)
  noise <- cosine_schedule(prop_epoch, 0.1, 1, 6)
  new_rho <- mod$target_rho + (mod$start_rho - mod$target_rho) * noise
  if (epoch >= anneal_epochs) new_rho <- mod$target_rho   # exact, no residual
  mod$rho <- new_rho

  ## --- Early stopping bookkeeping (stationary phase only) ---
  ## Losses recorded while rho is still annealing are NOT comparable to each other
  ## (different objective each epoch), so we ignore them entirely: best-tracking
  ## starts once the anneal has completed. A "real" improvement is a new best
  ## beating the prior best by > TOL (relative).
  if (epoch >= anneal_epochs) {
    if (total_loss < best_loss * (1 - tol)) {
      best_loss  <- total_loss
      best_epoch <- epoch
      ## Save the BEST model, not just the last one. Previously only the final-epoch
      ## model was kept, so the ancestral estimates were extracted from a model past
      ## its own optimum (and, for a diverged run, from a badly broken one).
      torch_save(mod, best_path)
    }
    if (early_stop && epoch >= min_epochs && (epoch - best_epoch) >= patience) {
      cat(sprintf("\n[early stopping] epoch %d: best loss %.6f @ epoch %d not improved by >%.3g%% for %d epochs (rho fixed at target since epoch %d)\n",
                  epoch, best_loss, best_epoch, tol * 100, epoch - best_epoch, anneal_epochs))
      stopped_early <- TRUE
      break
    }
  }

}

train_elapsed <- as.numeric(Sys.time() - train_start, units = "secs")
epochs_run <- epoch          # last epoch executed (== n_epoch unless stopped early)
final_loss <- total_loss     # total loss at the final epoch
per_epoch  <- train_elapsed / epochs_run
cat(sprintf("\n=== TRAINING COMPLETE ===\n"))
cat(sprintf("epochs run:          %d  (cap N_EPOCHS=%d, stopped_early=%s)\n",
            epochs_run, n_epoch, stopped_early))
cat(sprintf("final loss:          %.6f  (epoch %d)\n", final_loss, epochs_run))
cat(sprintf("best loss:           %.6f  (epoch %d)\n", best_loss, best_epoch))
cat(sprintf("total training time: %.1f s (%.2f min)\n", train_elapsed, train_elapsed / 60))
cat(sprintf("per-epoch mean:      %.3f s\n", per_epoch))
cat(sprintf("extrapolated 2500ep: %.1f min (%.2f h)\n",
            per_epoch * 2500 / 60, per_epoch * 2500 / 3600))

## Machine-readable convergence record for the sweep to aggregate
conv <- data.frame(
  out_tag = out_tag, seed = seed_int, epochs_run = epochs_run, n_epochs_cap = n_epoch,
  stopped_early = stopped_early, final_loss = final_loss,
  best_loss = best_loss, best_epoch = best_epoch, restored_best = restored_best,
  anneal_epochs = anneal_epochs, lr_epochs = lr_epochs, stop_grace = stop_grace,
  patience = patience, min_epochs = min_epochs, tol = tol,
  per_epoch_secs = per_epoch, total_secs = train_elapsed,
  tip_w = tip_weight, manifold_w = manifold_weight, code_w = code_weight,
  trophic_w = trophic_weight, root_w = root_weight,
  rho_start_mult = rho_start_mult, rho_target_div = rho_target_div,
  lambda = lambda_val
)
write.csv(conv, file.path(out_dir, "convergence.csv"), row.names = FALSE)
tryCatch({
  if (cuda_is_available()) {
    cat(sprintf("peak GPU mem allocated: %.2f GB\n",
                as.numeric(cuda_max_memory_allocated()) / 1e9))
    cat(sprintf("peak GPU mem reserved:  %.2f GB\n",
                as.numeric(cuda_max_memory_reserved()) / 1e9))
  }
}, error = function(e) cat("GPU mem query unavailable:", conditionMessage(e), "\n"))

options(torch.serialization_version = 2)

## --- Restore the BEST model before extracting anything -----------------------
## Everything downstream (ancestral estimates, scoring) must come from the best
## model on the STATIONARY objective, not from whatever the last epoch happened to
## leave behind. Without this, a run that drifted (or diverged) past its optimum
## exported ancestors from a worse model than it had already found.
restored_best <- FALSE
if (file.exists(best_path) && best_epoch > 0L) {
  mod <- torch_load(best_path)
  mod$to(device = if (cuda_is_available()) "cuda" else "cpu")
  restored_best <- TRUE
  cat(sprintf("restored BEST model from epoch %d (loss %.6f) for extraction\n",
              best_epoch, best_loss))
} else {
  cat("WARNING: no best checkpoint found — extracting from the final-epoch model\n")
}
torch_save(mod, file.path(out_dir, "mani_evo_mod_v3_bayesian.to"))

###############################################################################
## Post-training: extract predictions
###############################################################################

bill_test_dl <- dataloader(bill_ds, 100, shuffle = FALSE, drop_last = FALSE)
test_dat <- list()
i <- 0
coro::loop(for (b in bill_test_dl) {
  i <- i + 1
  with_no_grad({
    test_dat[[i]] <- purrr::map(purrr::list_flatten(mod(b)[-4]), ~.x$cpu()$detach())
  })
  print(i)
})

## Linear baseline (a=0, b=0): straight lines in latent space
linear_mod <- mod$clone()
with_no_grad({
  linear_mod$a$zero_()
  linear_mod$b$zero_()
})

test_dat_linear <- list()
i <- 0
coro::loop(for (b in bill_test_dl) {
  i <- i + 1
  with_no_grad({
    test_dat_linear[[i]] <- purrr::map(purrr::list_flatten(linear_mod(b)[-4]), ~.x$cpu()$detach())
  })
  print(i)
})

## Save curved predictions
edge_ids <- rownames(bill_tree_mat)
test_zs <- torch_cat(map(test_dat, 9)) %>% as.array()
z_list <- array_branch(test_zs, 1) %>%
  purrr::map(~ as.data.frame(t(.x)) %>%
               setNames(paste0("latent_", 1:ncol(.))))

z_tree_df <- tibble(edge = edge_ids) %>%
  mutate(z_seqs = z_list) %>%
  left_join(bill_pf %>%
              select(edge = label, is_tip,
                     starts_with("latent_")))

write_rds(z_tree_df, file.path(out_dir, "bill_vae_aces_16dim_v3_bayesian.rds"))

## Save linear baseline
test_zs_linear <- torch_cat(map(test_dat_linear, 9)) %>% as.array()
z_list_linear <- array_branch(test_zs_linear, 1) %>%
  purrr::map(~ as.data.frame(t(.x)) %>%
               setNames(paste0("latent_", 1:ncol(.))))

z_tree_df_linear <- tibble(edge = edge_ids) %>%
  mutate(z_seqs = z_list_linear) %>%
  left_join(bill_pf %>%
              select(edge = label, is_tip,
                     starts_with("latent_")))

write_rds(z_tree_df_linear, file.path(out_dir, "bill_vae_aces_16dim_v3_bayesian_linear.rds"))

cat("\n=== DONE: outputs in", out_dir, "===\n")
