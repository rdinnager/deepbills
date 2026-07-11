###############################################################################
## Post-training prediction extraction for v3 Bayesian model
##
## Loads the saved model from data/v3_bayesian/mani_evo_mod_v3_bayesian.to
## and extracts curved + linear baseline predictions.
##
## Run: Rscript R/.VAE_evo_model_v3_extract.R
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

out_dir <- "data/v3_bayesian"

## Load decoder (needed for model forward pass)
options(torch.serialization_version = 2)
bill_vae <- torch_load("data/bill_vae_w_trophic_v1.to")
bill_decoder <- bill_vae$decoder
bill_decoder$parameters %>% purrr::walk(function(param) param$requires_grad_(FALSE))
bill_decoder <- bill_decoder$cuda()

active_dims <- read_rds("data/active_dims_16dim.rds")

## Load phylogenetic data
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

## Dataset definition (needed for dataloader)
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
                            list(end_tree_mat, start_tree_mat, self$tip_tree_mat, self$tip_dat, self$blens[i], torch_tensor(i, device = "cuda"))
                         },
                         .length = function() {
                           nrow(self$end_tree_mat)
                         })

bill_ds <- bill_dataset(bill_tree_mat, start_mat, bill_tree_mat_tips, bill_tips, blens)

## Function definitions needed by the model
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

  mahalanobis_squared_tr <- jit_trace(mahalanobis_squared, z1, z2, v2, latvars_tens)

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

get_metric_tensor <- function(z, centroids, vars, lambda = 1e-2, rho, latvars_tens) {
  mh <- torch_exp(-(mahalanobis_squared[[1]](z, centroids, vars, latvars_tens) / (rho^2)))
  1 / (((1 / vars)$unsqueeze(1)$unsqueeze(-1) * mh$unsqueeze(3))$sum(dim = 2) + lambda)
}

get_manifold_energy <- function(vel, metric) {
  (vel * metric * vel)$sum(dim = 2)
}

get_manifold_dist <- function(vel, metric) {
  (vel * metric * vel)$sum(dim = 2)$sqrt()
}

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

## Centroids and vars (needed for metric tensor in forward pass)
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

## nn_module definition (needed for torch_load)
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

## Initialize mahalanobis JIT (needed before model forward pass)
bill_init_rates <- bill_pf %>%
  select(starts_with("rate_")) %>%
  as.matrix()

bill_dl_init <- dataloader(bill_ds, 200, shuffle = TRUE, drop_last = FALSE)
test <- dataloader_next(dataloader_make_iter(bill_dl_init))
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
## Load trained model
###############################################################################

cat("Loading trained model...\n")
mod <- torch_load(file.path(out_dir, "mani_evo_mod_v3_bayesian.to"))
mod <- mod$cuda()
cat("Model loaded.\n")

###############################################################################
## Extract curved predictions
###############################################################################

cat("Extracting curved predictions...\n")
bill_test_dl <- dataloader(bill_ds, 100, shuffle = FALSE, drop_last = FALSE)
test_dat <- list()
i <- 0
coro::loop(for (b in bill_test_dl) {
  i <- i + 1
  with_no_grad({
    test_dat[[i]] <- purrr::map(purrr::list_flatten(mod(b)[-4]), ~.x$cpu()$detach())
  })
  cat("Batch", i, "\n")
})

###############################################################################
## Extract linear baseline predictions (a=0, b=0)
###############################################################################

cat("Extracting linear baseline predictions...\n")
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
  cat("Linear batch", i, "\n")
})

###############################################################################
## Save curved predictions
###############################################################################

cat("Saving curved predictions...\n")
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

###############################################################################
## Save linear baseline predictions
###############################################################################

cat("Saving linear baseline predictions...\n")
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

cat("Done! Predictions saved to:\n")
cat("  ", file.path(out_dir, "bill_vae_aces_16dim_v3_bayesian.rds"), "\n")
cat("  ", file.path(out_dir, "bill_vae_aces_16dim_v3_bayesian_linear.rds"), "\n")
