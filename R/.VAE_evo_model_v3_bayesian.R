###############################################################################
## Riemannian Manifold Evolutionary Model — v3 (Bayesian / Energy formulation)
##
## This script modifies v2 to use a proper Bayesian formulation where
## minimizing the loss is equivalent to finding the MAP estimate of:
##
##   Likelihood:
##     Y_tips | z_tips ~ N(z_tips, sigma_tip^2 I)
##
##   Prior on root:
##     z_root ~ N(0, sigma_root^2 I)
##
##   Prior on evolutionary paths (Onsager-Machlup action for Riemannian BM):
##     -log p(path_e) = (1 / (2 sigma_evo^2)) * E_e
##     where E_e = (1/T_e) * integral_0^1 ||dz/ds||^2_G ds
##     is the path energy divided by branch length T_e.
##
##   Smoothness priors in decoded (phenotype) space:
##     -log p_code(path_e) = alpha_code * (1/T_e) * integral ||d(dec(z))/ds||^2 ds
##     -log p_troph(path_e) = alpha_troph * (1/T_e) * integral ||d(troph(z))/ds||^2 ds
##
## Key changes from v2:
##   1. get_manifold_energy: squared Riemannian norm (no sqrt)
##      — corresponds to the Onsager-Machlup action for Riemannian BM
##   2. Division by branch length T_e (not T_e/n_segs)
##      — proper BM scaling: longer branches allow more variance
##   3. Same changes applied to decoder losses
##
## The loss weights map to inverse-variance (precision) parameters:
##   tip_weight      ~ 1 / (2 * sigma_tip^2)
##   (implicit 1.0)  ~ 1 / (2 * sigma_evo^2)    [manifold energy weight]
##   code_weight     ~ alpha_code
##   trophic_weight  ~ alpha_trophic
##   root_weight     ~ 1 / (2 * sigma_root^2)
##
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

## Create output directory (avoids overwriting existing data/)
out_dir <- "data/v3_bayesian"
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
                           },
                           .getbatch = function(i) {
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

n_segs <- 50

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
##
## G(z)_jj = 1 / [ sum_i (1/sigma_ij^2) * exp(-d_Mah(z, mu_i)^2 / rho^2) + lambda ]
##
## Interpretation: G defines a position-dependent diffusion coefficient.
## Near VAE centroids (data-dense), G is small -> diffusion is fast.
## Far from centroids (unsupported), G is large -> diffusion is slow.
###############################################################################

get_metric_tensor <- function(z, centroids, vars, lambda = 1e-2, rho, latvars_tens) {
  mh <- torch_exp(-(mahalanobis_squared[[1]](z, centroids, vars, latvars_tens) / (rho^2)))
  1 / (((1 / vars)$unsqueeze(1)$unsqueeze(-1) * mh$unsqueeze(3))$sum(dim = 2) + lambda)
}

###############################################################################
## KEY CHANGE: Energy formulation (squared Riemannian norm, no sqrt)
##
## v2 used path length:  L = integral ||dz/dt||_G dt     (sqrt of quadratic form)
## v3 uses path energy:  E = integral ||dz/dt||^2_G dt   (quadratic form itself)
##
## The energy is the Onsager-Machlup action for Brownian motion on the
## Riemannian manifold. Its MAP estimate corresponds to:
##   p(path) ~ exp(-E / (2 * sigma_evo^2 * T_e))
##
## The sqrt in v2 gave an L1-like penalty (Laplace prior on displacement);
## the squared form here gives an L2 penalty (Gaussian prior = true BM).
###############################################################################

get_manifold_energy <- function(vel, metric) {
  ## Squared Riemannian norm: vel^T G vel (summed over dimensions)
  ## Returns shape: (batch, n_segments)
  (vel * metric * vel)$sum(dim = 2)
}

## Keep the old function available for comparison
get_manifold_dist <- function(vel, metric) {
  (vel * metric * vel)$sum(dim = 2)$sqrt()
}

###############################################################################
## Loss functions
###############################################################################

## Mean energy across edges and segments
energy_loss <- function(energies) {
  torch_mean(energies)
}

tip_loss <- function(tip_data, tip_recon) {
  ## Gaussian likelihood: -log p(Y|z) ~ ||Y - z||^2 / (2 sigma_tip^2)
  tip_dists <- torch_square(tip_data - tip_recon)
  torch_mean(tip_dists)
}

root_loss <- function(root_values) {
  ## Gaussian prior: -log p(root) ~ ||root||^2 / (2 sigma_root^2)
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
                            target_rho <- rho / 3
                            start_rho <- rho * 3
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

                            z_ends <- call_torch_function("torch__sparse_mm", x[[1]], self$rates, quiet = TRUE) +
                              self$root_values$unsqueeze(1)
                            z_starts <- call_torch_function("torch__sparse_mm", x[[2]], self$rates, quiet = TRUE) +
                              self$root_values$unsqueeze(1)
                            z_tips <- call_torch_function("torch__sparse_mm", x[[3]], self$rates, quiet = TRUE) +
                              self$root_values$unsqueeze(1)

                            zs <- get_segments(z_starts, z_ends, self$segs, self$a[x[[6]]], self$b[x[[6]]],
                                               self$bill_decoder, self$active_dims, device = self$device)

                            met <- get_metric_tensor(zs[[4]], centroids, vars, rho = self$rho, latvars_tens = self$latvars_tens)

                            ## KEY CHANGE: energy (squared norm) divided by branch length
                            ## E_e = (1/T_e) * integral ||dz/ds||^2_G ds
                            ## Discretized: sum_k ||Delta_z_k||^2_G / T_e
                            ## The n_segs factor from the quadrature cancels with
                            ## using torch_mean over segments (which divides by n_segs-1).
                            blens_expanded <- x[[5]]$unsqueeze(-1)

                            ## Manifold energy: Onsager-Machlup action for Riemannian BM
                            manifold_energies <- get_manifold_energy(zs[[3]], met) / blens_expanded
                            ## Decoded beak shape energy: BM prior in phenotype space
                            code_energies <- get_manifold_energy(zs[[1]], 1.0) / blens_expanded
                            ## Decoded trophic niche energy: BM prior in ecological niche space
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

###############################################################################
## Loss weights — these are the Bayesian precision parameters
##
## The total loss is:  -log p(theta | Y) up to a constant, where:
##   manifold_weight  * manifold_energy  ~  1 / (2 sigma_evo^2)   per edge
##   code_weight      * code_energy      ~  alpha_code             per edge
##   trophic_weight   * trophic_energy   ~  alpha_trophic          per edge
##   tip_weight       * tip_MSE          ~  1 / (2 sigma_tip^2)   for observations
##   root_weight      * ||root||^2       ~  1 / (2 sigma_root^2)  for root prior
##
## NOTE: Because the energy formulation produces smaller values than the
## path length formulation (squared small numbers are smaller), you may
## need to increase manifold_weight, code_weight, or trophic_weight
## relative to v2 to get comparable behavior. Start with these defaults
## and adjust based on the relative magnitudes of each loss term.
###############################################################################

manifold_weight <- 1.0
code_weight <- 1/64
trophic_weight <- 1/10
tip_weight <- 10
root_weight <- 1/100

mod <- mani_evo_mod(n_rates = nrow(bill_init_rates), n_dim = 16,
                    centroids, vars,
                    n_segs = n_segs, init_rate = bill_init_rates,
                    bill_decoder = bill_decoder, active_dims = active_dims,
                    latvars = lat_vars,
                    device = "cuda")
mod <- mod$cuda()

write_rds(as.numeric(mod$target_rho$cpu()), file.path(out_dir, "estimated_rho_16dim_rho_schedule_v3.rds"))

## Initialize mahalanobis function (needs a test forward pass for JIT)
test <- dataloader_next(dataloader_make_iter(bill_dl))
init_rates_test <- torch_tensor(bill_init_rates, device = "cuda")
z_ends_test <- call_torch_function("torch__sparse_mm", test[[1]], init_rates_test, quiet = TRUE)
z_starts_test <- call_torch_function("torch__sparse_mm", test[[2]], init_rates_test, quiet = TRUE)
a_test <- torch_randn(z_starts_test$size()[1], z_starts_test$size()[2], device = "cuda")
b_test <- torch_randn(z_starts_test$size()[1], z_starts_test$size()[2], device = "cuda")
segs_test <- torch_arange(0, 1, 1 / n_segs, device = "cuda")
zs_test <- get_segments(z_starts_test, z_ends_test, segs_test, a_test, b_test,
                        bill_decoder, active_dims, device = "cuda")
mahalanobis_squared <- mahalanobis_squared_fun(zs_test[[4]], centroids, vars, latvars_tens)

###############################################################################
## Rho annealing schedule (posterior tempering)
##
## Large rho -> nearly flat (Euclidean) metric -> vague prior
## Small rho -> sharp manifold structure -> informative prior
## Cosine annealing avoids bad local optima early in training.
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

n_epoch <- 2500
lr <- 0.02
save_every <- 50

optim1 <- optim_adam(mod$parameters, lr = lr)
scheduler <- lr_one_cycle(optim1, max_lr = lr,
                          epochs = n_epoch, steps_per_epoch = 1,
                          cycle_momentum = FALSE)

optim1$zero_grad()

## Initialize loss history CSV
loss_csv <- file.path(out_dir, "loss_history.csv")
cat("epoch,total_loss,tip_loss,manifold_energy,code_energy,trophic_energy,root_loss,rho,lr\n",
    file = loss_csv)

checkpoint_every <- 100

test_dat <- list()

for(epoch in 1:n_epoch) {

  optim1$zero_grad()
  total_loss <- 0
  total_recon_loss <- 0
  total_manifold_loss <- 0
  total_code_loss <- 0
  total_trophic_loss <- 0
  total_root_loss <- 0

  i <- 0
  coro::loop(for (b in bill_dl) {
      i <- i + 1

      res <- mod(b)

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

      total_loss <- total_loss + as.numeric(loss$cpu())
      total_recon_loss <- total_recon_loss + as.numeric(res[[5]]$cpu())
      total_manifold_loss <- total_manifold_loss + as.numeric(res[[1]]$cpu())
      total_code_loss <- total_code_loss + as.numeric(res[[2]]$cpu())
      total_trophic_loss <- total_trophic_loss + as.numeric(res[[3]]$cpu())
      total_root_loss <- total_root_loss + as.numeric(res[[6]]$cpu())

  })

  current_rho <- as.numeric(mod$rho$cpu())
  current_lr <- optim1$param_groups[[1]]$lr

  cat("Epoch: ", epoch,
      "    loss: ", as.numeric(total_loss),
      "    tip recon loss: ", as.numeric(total_recon_loss),
      "    manifold energy: ", as.numeric(total_manifold_loss),
      "    code energy: ", as.numeric(total_code_loss),
      "    trophic energy: ", as.numeric(total_trophic_loss),
      "    root loss: ", as.numeric(total_root_loss),
      "    rho: ", current_rho,
      "\n")

  ## Append to loss history CSV
  cat(paste(epoch, total_loss, total_recon_loss, total_manifold_loss,
            total_code_loss, total_trophic_loss, total_root_loss,
            current_rho, current_lr, sep = ","), "\n",
      file = loss_csv, append = TRUE)

  ## Save checkpoint
  if(epoch %% checkpoint_every == 0) {
    torch_save(mod, file.path(out_dir, paste0("checkpoint_epoch_", epoch, ".to")))
  }
  torch_save(mod, file.path(out_dir, "checkpoint_latest.to"))

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

  ## Update rho (posterior tempering schedule)
  prop_epoch <- epoch / n_epoch
  noise <- cosine_schedule(prop_epoch, 0.1, 1, 6)
  new_rho <- mod$target_rho + (mod$start_rho - mod$target_rho) * noise
  mod$rho <- new_rho

}

options(torch.serialization_version = 2)
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
