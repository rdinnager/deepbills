library(tidyverse)
library(phyf)
library(fibre)
library(Matrix)
library(MCMCglmm)
library(data.table)
library(INLA)
library(Morpho)

beak_folder <- "data/landmarks_oriented/all_landmarks"
get_aligned_beak_landmarks <- function(beak_folder) {
  beak_dats <- list.files(beak_folder, full.names = TRUE)
  specs <- word(fs::path_ext_remove(basename(beak_dats)),
                end = 2L, sep = fixed("_"))
  landmarks <- map2(beak_dats, specs,
                   ~ fread(.x, col.names = c("x", "y", "z")) |>
                     mutate(Species = .y),
                   .progress = TRUE) |>
    list_rbind()
  landmarks
}
landmarks <- get_aligned_beak_landmarks(beak_folder)
write_csv(landmarks, "data/landmarks_oriented/all_landmarks.csv")

landmark_df <- landmarks |>
  group_by(Species) |>
  mutate(point = 1:n()) |>
  ungroup() |>
  pivot_longer(c(-Species, -point), names_to = "coord", values_to = "value") |>
  mutate(point_name = paste0(coord, "_", point)) |>
  select(Species, value, point_name) |>
  pivot_wider(names_from = point_name, values_from = value)

bird_beak_coords <- bird_beak_codes |>
  select(-starts_with("latent_")) |>
  left_join(landmark_df, by = c("label" = "Species"))

bird_tree <- pf_as_phylo(bird_beak_coords)

bird_coord_mat <- landmark_df |>
  select(-Species) |>
  as.matrix() 

## First do PCA whitening on data to remove strong correlations
bird_coord_mat_pca <- prcompfast(bird_coord_mat, scale. = TRUE)
bird_coord_mat_rot <- bird_coord_mat_pca$x

bird_tree_prec <- inverseA(bird_tree, nodes = "TIPS")
tree_match <- match(landmark_df$Species, bird_tree_prec$node.names)
bird_tree_prec <- bird_tree_prec$Ainv[tree_match, tree_match]
bird_tree_chol <- t(chol(bird_tree_prec))

bird_coord_mat_tr <- bird_tree_chol %*% bird_coord_mat_rot

landmark_df_tr <- landmark_df |>
  select(Species) |>
  bind_cols(as.data.frame(as.matrix(bird_coord_mat_tr)))

# landmark_df_tr_long <- landmark_df_tr |>
#   pivot_longer(-Species, names_to = c("coord", "point"), names_sep = "_") |>
#   pivot_wider(names_from = coord, values_from = value)
# 
# plot((landmark_df_tr_long |> filter(Species == landmark_df_tr_long$Species[1]) |> select(y, z)))
# plot((landmarks |> filter(Species == landmark_df_tr_long$Species[1]) |> select(y, z)))
# 
# plot(landmark_df |> select(x_1, x_2))
# plot(landmark_df_tr |> select(PC1, PC2))
# 
# plot(landmark_df |> select(x_44, y_44))
# plot(landmark_df_tr |> select(x_44, y_44))
# 
# total_sd <- sd(c(landmark_df_tr_long$x, landmark_df_tr_long$y, landmark_df_tr_long$z))

landmark_tr_pf <- bird_beak_codes |>
  select(label, is_tip, phlo) |>
  left_join(landmark_df_tr, by = c("label" = "Species")) |>
  mutate(coord_mat = scale(as.matrix(across(starts_with("PC")))))
  
hyper1 = list(list(prec=list(initial = 0, fixed = TRUE)))

mod <- fibre(coord_mat ~ bre_brownian(phlo), data = landmark_tr_pf, 
             verbose = 2, engine_option = list(verbose = TRUE,
                                               num.threads = 8,
                                               control.family = list(hyper = list(hyper = list(prec = list(prior = "pc.prec", initial = 4, fixed = TRUE))))))
  
write_rds(mod, "output/landmark_phylo_model.rds")