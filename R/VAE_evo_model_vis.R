library(tidyverse)
library(torch)
library(fibre)
library(phyf)
library(targets)
library(rgl)
library(ggtree)

beak_mod <- load_bird_beak_model()
beak_mod <- beak_mod$cuda()

z_tree_df <- read_rds("data/evo_predictions_all_edges_v2.rds")

niche_pal <- tar_read(niche_pal)

scaling <- read_rds("data/code_dat_means_sds.rds")

tree_df <- phyf::bird_beak_codes |>
  pf_as_phylo() |>
  ggtree::fortify()

trophic_tree <- z_tree_df |>
  select(edge, prtroph, time_seqs, true_trophic) |>
  unnest_longer(-edge) |>
  left_join(tree_df |>
              select(edge = label, y)) 

tree_vert <- trophic_tree |>
  dplyr::filter(time_seqs != max(time_seqs)) |>
  mutate(timeslice = as.factor(time_seqs)) |>
  arrange(y)

ggplot(trophic_tree, aes(x = time_seqs, y = y, color = prtroph)) +
  geom_path(aes(group = edge)) +
  geom_path(aes(group = timeslice), data = tree_vert) +
  geom_point(aes(color = prtroph), data = trophic_tree |>
               group_by(edge) |>
               dplyr::filter(time_seqs == max(time_seqs)) |>
               dplyr::filter(!is.na(true_trophic))) +
  theme_tree2()

ggplot(trophic_tree, aes(x = time_seqs, y = y, color = prtroph)) +
  geom_path(aes(group = edge)) +
  geom_path(aes(group = timeslice), data = tree_vert) +
  geom_point(aes(color = prtroph), data = trophic_tree |>
               group_by(edge) |>
               dplyr::filter(time_seqs == max(time_seqs)) |>
               dplyr::filter(!is.na(true_trophic))) +
  coord_polar("y") +
  theme_tree2()

########## make meshes ###############
# how much distance is traversed by each edge?
zseqs <- z_tree_df$z_seqs[[1]]
time_seqs <- z_tree_df$time_seqs[[1]]
calc_rate_along <- function(zseqs, time_seqs) {
  zdiffs <- diff(as.matrix(zseqs))
  z_dist <- sum(sqrt(rowSums(zdiffs^2)))
  z_rate <- z_dist / (time_seqs[length(time_seqs)] - time_seqs[1])
  z_rate
}

calc_dist_along <- function(zseqs) {
  zdiffs <- diff(as.matrix(zseqs))
  z_dist <- sum(sqrt(rowSums(zdiffs^2)))
  z_dist
}

rates_along <- map2_dbl(z_tree_df$z_seqs, z_tree_df$time_seqs, 
                        ~ calc_rate_along(.x, .y),
                        .progress = TRUE)

sum(rates_along)

dist_along <- map_dbl(z_tree_df$z_seqs,
                       ~ calc_dist_along(.x),
                       .progress = TRUE)

num_meshes <- 10000
dist_split <- sum(dist_along) / num_meshes

num_splits <- dist_along %/% dist_split

indexes <- map(num_splits,
               ~ floor(100 * seq(0, 1, length.out = .x + 2))[-1])

prcodes <- map2(z_tree_df$prcodes, indexes, 
               ~ .x[.y, , drop = FALSE])

prcode <- prcodes[[1]][1, , drop = FALSE]
make_meshes <- function(prcode, scaling) {
  mesh <- try(beak_mod$get_mesh(torch_tensor((prcode * scaling$sds) + scaling$means, 
                                         device = "cuda"),
                            resolution = 150,
                            smooth = FALSE))
}
meshes <- map(array_branch(do.call(rbind, prcodes), 1), 
              ~ make_meshes(matrix(.x, nrow = 1), scaling = scaling),
              .progress = TRUE)

write_rds(meshes, "data/z_tree_meshes_v2.rds")

open3d()
shade3d(meshes[[1]], col = "gold", alpha = 0.25)
Sys.sleep(1)
#clear3d()
shade3d(meshes[[2]], col = "gold", alpha = 0.25)
Sys.sleep(1)
#clear3d()
shade3d(meshes[[3]], col = "gold", alpha = 0.25)
Sys.sleep(1)
#clear3d()
shade3d(meshes[[4]], col = "gold", alpha = 0.25)
Sys.sleep(1)
#clear3d()
shade3d(meshes[[5]], col = "gold", alpha = 0.25)

