library(tidyverse)
library(torch)
library(fibre)
library(phyf)
library(targets)

avonet <- tar_read(bird_beak_avonet)

options(torch.serialization_version = 2)
vae <- torch_load("data/bill_vae_w_trophic_v1.to")
vae <- vae$cuda()

beak_mod <- load_bird_beak_model()
beak_mod <- beak_mod$cuda()

z_tree_df <- read_rds("data/bill_vae_aces_16dim_noise_schedule_v2.rds")

trophic_dat <- read_rds("data/trophic_dat.rds")

z_tree_df <- z_tree_df %>%
  left_join(trophic_dat |> select(edge = label, true_trophic = trophic_niche))

active_dims <- read_rds("data/active_dims_16dim.rds")

rgl::open3d()
for(i in 1:nrow(z_tree_df)) {
  rgl::lines3d(z_tree_df$z_seqs[[i]]$latent_1,
               z_tree_df$z_seqs[[i]]$latent_2,
               z_tree_df$z_seqs[[i]]$latent_3,
               col = "black",
               lwd = 1)
}
rgl::points3d(z_tree_df$true_latent_1,
              z_tree_df$true_latent_2,
              z_tree_df$true_latent_3,
              col = "red",
              size = 5)
rgl::texts3d(z_tree_df$true_latent_1[z_tree_df$is_tip],
              z_tree_df$true_latent_2[z_tree_df$is_tip],
              z_tree_df$true_latent_3[z_tree_df$is_tip],
             z_tree_df$edge[z_tree_df$is_tip],
             adj = c(0, 0))

zseq <- z_tree_df$z_seqs[[1]]
trophic_levs <- levels(trophic_dat$trophic_niche)
decode_zseqs <- function(zseq, active_dims, trophic_levs) {
  z_mat <- matrix(0, nrow = nrow(zseq), ncol = 64)
  z_mat[ , active_dims] <- as.matrix(zseq)
  with_no_grad({
    z_tensor <- torch_tensor(z_mat, device = "cuda")
    x <- vae$decoder(z_tensor)
    pred_troph <- torch_argmax(nnf_softmax(x$out_trophic, 2L), dim = 2L)
    pred_codes <- x$out_codes
  })
  list(troph = trophic_levs[as.matrix(pred_troph$cpu())],
       codes = as.matrix(pred_codes$cpu()))
}

pred_seqs <- map(z_tree_df$z_seqs, ~ decode_zseqs(.x, active_dims, trophic_levs),
                 .progress = TRUE)

z_tree_df <- z_tree_df %>%
  mutate(prcodes = map(pred_seqs, ~ .x$codes),
         prtroph = map(pred_seqs, ~ .x$troph))

z_tree_df <- z_tree_df |>
  mutate(prcode_last = map(prcodes, ~ .x[nrow(.x), , drop = FALSE] |>
           as.data.frame()) |>
           list_rbind(),
         prtroph_last = map_chr(prtroph, ~ .x[length(.x)]))

z_tree_df <- z_tree_df |>
  left_join(avonet |>
              select(edge = label, starts_with("latent_code_")) |>
              rename_with(~ str_replace_all(.x, "latent_code_",
                                            "true_latent_code64_")))

bill_edge_trajs <- read_rds("data/bill_edge_trajs_16dim.rds")

z_tree_df <- z_tree_df |>
  left_join(bill_edge_trajs |>
              select(edge = end,
                     start_time,
                     end_time)) |>
  mutate(time_seqs = map2(start_time, end_time,
                          ~ seq(.x, .y, length.out = 50)))

write_rds(z_tree_df, "data/evo_predictions_all_edges_v2.rds")
