library(tidyverse)
library(torch)
library(fibre)
library(phyf)
library(targets)
library(rgl)
library(ggtree)
library(ggrepel)
library(ggraph)
library(nat)

beak_mod <- load_bird_beak_model()
beak_mod <- beak_mod$cuda()

avonet <- tar_read(bird_beak_avonet)

z_tree_df <- read_rds("data/evo_predictions_all_edges_v2.rds")

niche_pal <- tar_read(niche_pal)
names(niche_pal) <- levels(as.factor(avonet$Trophic.Niche))

scaling <- read_rds("data/code_dat_means_sds.rds")

tree_df <- phyf::bird_beak_codes |>
  pf_as_phylo() |>
  ggtree::fortify()

trophic_tree <- z_tree_df |>
  select(edge, prtroph, time_seqs, true_trophic) |>
  unnest_longer(-edge) |>
  left_join(tree_df |>
              select(edge = label, y)) |>
  left_join(avonet |>
              select(edge = label, BLFamilyEnglish, Order))

fam_summ <- trophic_tree |>
  group_by(Order, BLFamilyEnglish) |>
  summarise(count = n() / 50) |>
  group_by(Order) |>
  mutate(ord_count = sum(count))

small_ords <- fam_summ |>
  filter(ord_count < 50) |>
  group_by(Order) |>
  mutate(prop = count / sum(count)) |>
  arrange(desc(prop)) |>
  summarise(lab = paste(na.omit(BLFamilyEnglish[prop > 0.05][1:3]), collapse = ", ")) 

big_fams <- fam_summ |>
  filter(ord_count >= 50) |>
  filter(count > 20) |>
  mutate(lab = BLFamilyEnglish)

choose_from <- avonet |>
  left_join(small_ords |>
              mutate(choose = 1)) |>
  filter(!is.na(choose)) |>
  bind_rows(avonet |>
              left_join(big_fams |>
                          select(BLFamilyEnglish, lab) |>
                          mutate(choose = 1)) |>
              filter(!is.na(choose))) 

choose_df <- choose_from |>
  group_by(lab, Trophic.Niche) |>
  sample_n(1) |>
  group_by(lab) |>
  sample_n(2, replace = TRUE) |>
  distinct(lab, .keep_all = TRUE) 

new_tree <- phyf::bird_beak_codes |>
  filter(label %in% choose_df$label) |>
  pf_as_phylo()



ord_summ <- trophic_tree |>
  group_by(Order, BLFamilyEnglish) |>
  summarize(count = n() / 50) |>
  mutate(prop = count / sum(count)) |>
  group_by(Order) |>
  arrange(desc(prop)) |>
  summarise(count = sum(count),
            fams = list(as.character(na.omit(BLFamilyEnglish[prop > 0.05][1:3])))) |>
  rowwise() |>
  mutate(order_label = paste(fams, collapse = ", "))

ord_tree_dat <- trophic_tree |>
  mutate(tip_close = abs(max(time_seqs) - time_seqs)) |>
  dplyr::filter(edge %in% phyf::bird_beak_codes$label[phyf::bird_beak_codes$is_tip],
                tip_close < 0.001) |>
  left_join(ord_summ |>
              select(Order, order_label)) |>
  mutate(time_seqs1 = time_seqs + 8,
         time_seqs2 = time_seqs1 + 2) |>
  group_by(Order) |>
  arrange(y) |>
  ungroup() |>
  drop_na() |>
  mutate(order_label = ifelse(Order == "PASSERIFORMES", "Passerines", order_label)) |>
  mutate(time_seqs2 = time_seqs2 + 0.1*as.numeric(as.factor(Order))) 

rl <- rle(ord_tree_dat$Order)
rl$values <- paste0("section_", seq_along(rl$values))
order_section <- inverse.rle(rl)
ord_tree_dat <- ord_tree_dat %>%
  mutate(order_section = order_section)

ord_tree_labs <- ord_tree_dat |>
  group_by(Order, order_label, order_section) |>
  summarise(y_mid = mean(y),
            time_mid = mean(time_seqs2)) |>
  group_by(Order) |>
  mutate(sec_num = as.numeric(as.factor(order_section))) |>
  mutate(order_label = ifelse(sec_num > 1, "", order_label),
         y_mid = ifelse(sec_num > 1, NA, y_mid))

# ord_tree_labs <- bind_rows(ord_tree_labs,
#                            ord_tree_dat |>
#                              select(Order, order_label,
#                                     y_mid = y,
#                                     time_mid = time_seqs2) |>
#                              mutate(order_label = ""))
  
tree_vert <- trophic_tree |>
  mutate(tip_close = abs(max(time_seqs) - time_seqs)) |>
  dplyr::filter(tip_close > 0.001) |>
  #dplyr::filter(!edge %in% phyf::bird_beak_codes$label[phyf::bird_beak_codes$is_tip]) |>
  mutate(timeslice = as.factor(time_seqs)) |>
  arrange(y)

# tt <- tree_vert |>
#   group_by(timeslice) |>
#   filter(n() == 2) |>
#   summarize(range = abs(y[2] - y[1]),
#             ys = list(y),
#             trophs = list(prtroph))

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
  scale_color_manual(name = "Trophic Niche", values = niche_pal) +
  guides(fill = "none") +
  theme_tree() +
  theme(legend.position = c(0, 0.5))

ggplot(trophic_tree, aes(x = time_seqs, y = y, color = prtroph)) +
  geom_path(aes(group = edge)) +
  geom_path(aes(group = timeslice), data = tree_vert) +
  geom_point(aes(color = prtroph), data = trophic_tree |>
               group_by(edge) |>
               dplyr::filter(time_seqs == max(time_seqs)) |>
               dplyr::filter(!is.na(true_trophic))) +
  geom_rect(aes(xmin = time_seqs1, xmax = time_seqs2, 
                ymin = y - 0.5, ymax = y + 0.5, fill = order_label), 
            data = ord_tree_dat, inherit.aes = FALSE) +
  scale_x_continuous(expand = c(0, 0), limits = c(0, 500)) +
  coord_polar("y") +
  geom_text_repel(aes(label = order_label, x = time_mid, y = y_mid), data = ord_tree_labs,
                  inherit.aes = FALSE,
                  min.segment.length = 0,
                  xlim = c(160, NA),
                  verbose = TRUE,
                  hjust = 0.5,
                  force = 1.5,
                  max.iter = 100000,
                  max.overlaps = 100,
                  max.time = 5) +
  scale_color_discrete(name = "Trophic Niche") +
  guides(fill = "none") +
  theme_tree() +
  theme(legend.position = c(0, 0.5))

trophic_rle <- map(z_tree_df$prtroph,
                   rle)
transitions <- map_if(trophic_rle,
                      ~ length(.x) > 1,
                      ~ cbind(.x$values[-length(.x$values)],
                              .x$values[-1]) |>
                        as.data.frame(),
                   .else = NULL) %>%
  list_rbind() 

transition_counts <- transitions |>
  group_by(V1, V2) |>
  summarise(n = n()) |>
  ungroup() 

transition_mat <- transition_counts |>
  mutate(n = n / sum(n)) |>
  pivot_wider(names_from = V2, values_from = n,
              values_fill = 0) 

trans_mat <- as.matrix(transition_mat |> select(-V1))
rownames(trans_mat) <- transition_mat$V1

trans_mat <- trans_mat[rownames(trans_mat), rownames(trans_mat)]

#rownames(trans_mat) <- colnames(trans_mat) <- NULL

diagram::plotmat(trans_mat, arr.lwd = floor(10 * (trans_mat / max(trans_mat)) + 1),
                 box.size = 0.035, arr.pos = 0.8, endhead = TRUE,
                 shadow.size = 0.005, dtext = 0.5)

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
               ~ floor(50 * seq(0, 1, length.out = .x + 2))[-1])

prcodes <- map2(z_tree_df$prcodes, indexes, 
               ~ .x[.y, , drop = FALSE])

# prcode <- prcodes[[1]][1, , drop = FALSE]
# make_meshes <- function(prcode, scaling) {
#   mesh <- try(beak_mod$get_mesh(torch_tensor((prcode * scaling$sds) + scaling$means, 
#                                          device = "cuda"),
#                             resolution = 150,
#                             smooth = FALSE))
# }

write_rds(meshes, "data/z_tree_meshes_v2.rds")

meshes <- read_rds("data/z_tree_meshes_v2.rds")

rows <- rep(seq_along(prcodes), times = map_int(prcodes, nrow))

all_tree_df <- z_tree_df |>
  mutate(meshes = tibble(row = rows,
                         mesh = meshes) |>
           group_by(row) |>
           summarise(mesh = list(mesh)) |>
           pull(mesh))

write_rds(all_tree_df, "data/all_tree_df.rds")

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


##### visualize first few dimensions ######
traj_df <- z_tree_df |>
  select(edge, z_seqs, prtroph) |>
  unnest(cols = c(z_seqs, prtroph))

just_pnts <- traj_df |>
               group_by(edge) |>
               slice_tail(n = 1) |>
               ungroup()

p1 <- ggplot(traj_df, aes(x = latent_1, y = latent_2, color = prtroph)) +
  geom_path(aes(group = edge), alpha = 0.6) +
  geom_point(aes(latent_1, latent_2, color = prtroph),
             data = just_pnts,
             size = 0.2,
             inherit.aes = FALSE,
             alpha = 0.4) +
  scale_colour_manual(values = niche_pal, name = "Trophic Niche") +
  theme_minimal()

p1

p2 <- ggplot(traj_df, aes(x = latent_3, y = latent_4, color = prtroph)) +
  geom_path(aes(group = edge), alpha = 0.6) +
  geom_point(aes(latent_3, latent_4, color = prtroph),
             data = just_pnts,
             size = 0.2,
             inherit.aes = FALSE,
             alpha = 0.4) +
  scale_colour_manual(values = niche_pal, name = "Trophic Niche") +
  theme_minimal()

p2

p3 <- ggplot(traj_df, aes(x = latent_15, y = latent_16, color = prtroph)) +
  geom_path(aes(group = edge), alpha = 0.6) +
  geom_point(aes(latent_15, latent_16, color = prtroph),
             data = just_pnts,
             size = 0.2,
             inherit.aes = FALSE,
             alpha = 0.4) +
  scale_colour_manual(values = niche_pal, name = "Trophic Niche") +
  theme_minimal()

p3

open3d()
for(i in 1:nrow(z_tree_df)) {
  lines3d(z_tree_df$z_seqs[[i]][ , 1:3], col = niche_pal[z_tree_df$prtroph[[i]]])
}

open3d()
for(i in 1:nrow(z_tree_df)) {
  lines3d(z_tree_df$z_seqs[[i]][ , 4:6], col = niche_pal[z_tree_df$prtroph[[i]]])
}


####### rgl trees #########

all_tree_df <- read_rds("data/all_tree_df.rds")

specs <- choose_from |>
  filter(Order == "PSITTACIFORMES")

bird_tree <- phyf::bird_beak_codes |>
  filter(label %in% specs$label) |>
  pf_as_phylo()

edges <- c(bird_tree$node.label, bird_tree$tip.label)
  
z_tree_df_red <- all_tree_df |>
  mutate(indexes = indexes) |>
  filter(edge %in% edges)

max_time <- max(z_tree_df_red$time_seqs |>
                  map_dbl(max))

rgl_tree <- z_tree_df_red |>
  #select(edge, prtroph, time_seqs, true_trophic) |>
  #unnest_longer(-edge) |>
  left_join(tree_df |>
              select(edge = label, y)) |>
  mutate(y = y / (max(y) - min(y)),
         time_seqs = map(time_seqs, ~.x / max_time)) 

tree_v <- rgl_tree |>
  select(edge, y, time_seqs, prtroph) |>
  mutate(first_time = map_dbl(time_seqs, ~.x[1]),
         troph = map_chr(prtroph, ~.x[1])) |>
  select(edge, y, first_time, troph) |>
  group_by(first_time) |>
  group_split() 

lens <- map_int(tree_v, nrow)


 pan3d <- function(button, dev = cur3d(), subscene = currentSubscene3d(dev)) {
   start <- list()
   
   begin <- function(x, y) {
     activeSubscene <- par3d("activeSubscene", dev = dev)
     start$listeners <<- par3d("listeners", dev = dev, subscene = activeSubscene)
     for (sub in start$listeners) {
       init <- par3d(c("userProjection","viewport"), dev = dev, subscene = sub)
       init$pos <- c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5)
       start[[as.character(sub)]] <<- init
     }
   }
   
   update <- function(x, y) {
     for (sub in start$listeners) {
       init <- start[[as.character(sub)]]
       xlat <- 2*(c(x/init$viewport[3], 1 - y/init$viewport[4], 0.5) - init$pos)
       mouseMatrix <- translationMatrix(xlat[1], xlat[2], xlat[3])
       par3d(userProjection = mouseMatrix %*% init$userProjection, dev = dev, subscene = sub )
      }
   }
   rgl.setMouseCallbacks(button, begin, update, dev = dev, subscene = subscene)
   cat("Callbacks set on button", button, "of RGL device", dev, "in subscene", subscene, "\n")
 }

open3d()

for(i in 1:nrow(rgl_tree)) {
  lines3d(rgl_tree$time_seqs[[i]], 0, rgl_tree$y[[i]], col = niche_pal[rgl_tree$prtroph[[i]]], lwd = 5)
}
for(i in 1:length(tree_v)) {
  lines3d(tree_v[[i]]$first_time, 0, tree_v[[i]]$y, col = niche_pal[tree_v[[i]]$troph], lwd = 5)
}
for(i in 1:nrow(rgl_tree)) {
  for(j in 1:length(rgl_tree$meshes[[i]])) {
    moved <- rgl_tree$meshes[[i]][[j]] |>
      scale3d(0.005, 0.005, 0.005) |>
      translate3d(rgl_tree$time_seqs[[i]][rgl_tree$indexes[[i]][[j]]], 0.001, rgl_tree$y[[i]])
    shade3d(moved, col = niche_pal[rgl_tree$prtroph[[i]]][rgl_tree$indexes[[i]][[j]]])
  }
}
pan3d(3)

######### Ancestral beak ###########

anc_z <- all_tree_df |>
  filter(edge == "Node2") |>
  pull(prcodes)

mesh <- beak_mod$get_mesh(torch_tensor((as.matrix(anc_z[[1]])[1, , drop = FALSE] * scaling$sds) + scaling$means,
                                       device = "cuda"),
                          resolution = 250,
                          smooth = TRUE)

shade3d(mesh, col = niche_pal["Invertivore"])

