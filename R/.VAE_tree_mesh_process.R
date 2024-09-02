library(tidyverse)
library(rgl)
library(uwot)
library(targets)

avonet <- tar_read(bird_beak_avonet)

trophic_dat <- as.factor(avonet$Trophic.Niche)
trophic_levs <- levels(trophic_dat)
trophic_nums <- seq_along(trophic_levs)
names(trophic_nums) <- trophic_levs

niche_pal <- tar_read(niche_pal)
names(niche_pal) <- trophic_levs

options(rstudio.help.showDataPreview = FALSE)

all_tree_df <- read_rds("data/all_tree_df.rds")

shade3d(all_tree_df$meshes[[100]][[1]])

shade3d(all_tree_df$meshes[[200]][[1]], col = "gold", alpha = 0.25)
shade3d(all_tree_df$meshes[[200]][[2]], col = "gold", alpha = 0.25)
shade3d(all_tree_df$meshes[[200]][[3]], col = "gold", alpha = 0.25)
shade3d(all_tree_df$meshes[[200]][[4]], col = "gold", alpha = 0.25)

big_z_mat <- do.call(rbind, all_tree_df$z_seqs) |>
  as.matrix()

embedding <- uwot::umap(big_z_mat, n_components = 3)
umap_df <- as.data.frame(embedding) |>
  mutate(rows = rep(seq_len(nrow(all_tree_df)), each = 50)) |>
  group_by(rows) |>
  summarise(umaps = list(cbind(V1, V2, V3)))

indexes <- map(all_tree_df$meshes,
               ~ floor(50 * seq(0, 1, length.out = length(.x) + 1))[-1])

all_tree_df2 <- all_tree_df |>
  mutate(indexes = indexes,
         umaps = umap_df$umaps) |>
  mutate(umaps = map2(umaps, indexes,
                      ~ .x[.y, , drop = FALSE]),
         times = map2(time_seqs, indexes,
                      ~ .x[.y]),
         prtroph = map2(prtroph, indexes,
                      ~ .x[.y])) |>
  select(edge, prtroph, times, umaps) |>
  mutate(col = map(prtroph, ~ niche_pal[.x])) |>
  unnest(cols = c(prtroph, times, umaps, col))

umaps_df <- all_tree_df |>
  mutate(umaps = map(umap_df$umaps,
                     ~ as.data.frame(.x) |>
                       bind_rows(tibble(V1 = NA, V2 = NA, V3 = NA))),
         time_seqs = map(time_seqs, ~ c(.x, NA))) |>
  select(edge, umaps, time_seqs) |>
  unnest(cols = c(umaps, time_seqs))

time_seq <- seq(0, max(all_tree_df2$times), length.out = 1000)
time_mat <- cbind(time_seq[-length(time_seq)], time_seq[-1])

for(i in 1:length(time_seq)) {
  new_df <- umaps_df |>
    filter(time_seqs <= time_seq[i])
  lines3d(new_df |>
            select(V1, V2, V3))
  new_df <- new_df |>
    group_by(edge) |>
    slice_tail(n = 1)
  new_df <- all_tree_df2 |>
    filter(times <= time_seq[i])
}
  
for(i in seq_along(umap_df$umaps)) {
  
  lines3d(umap_df$umaps[[i]])
  for(j in seq_along(indexes[[i]])) {
    mesh_place <- all_tree_df$meshes[[i]][[j]] |>
    scale3d(0.5, 0.5, 0.5) |>
    translate3d(umap_df$umaps[[i]][indexes[[i]][[j]], 1],
                umap_df$umaps[[i]][indexes[[i]][[j]], 2],
                umap_df$umaps[[i]][indexes[[i]][[j]], 3])
    shade3d(mesh_place, col = "gold")
  }
}






