
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

get_landmark_pf <- function(landmark_folder) {
  landmarks <- get_aligned_beak_landmarks(landmark_folder)
  
  ## rescale landmarks to same bounding box
  ## scale by the range into -1, 1
  landmarks <- landmarks |>
    group_by(Species) |>
    mutate(x_range = max(x) - min(x),
           y_range = max(y) - min(y),
           z_range = max(z) - min(z),
           range = max(x_range, y_range, z_range),
           x = ((x - min(x)) / x_range * 2 - 1) * (x_range / range), 
           y = ((y - min(y)) / y_range * 2 - 1) * (y_range / range),
           z = ((z - min(z)) / z_range * 2 - 1) * (z_range / range)) |>
    ungroup() |>
    select(Species, x, y, z)
  
  landmark_df <- landmarks |>
    group_by(Species) |>
    mutate(point = 1:n()) |>
    ungroup() |>
    pivot_longer(c(-Species, -point), names_to = "coord", values_to = "value") |>
    mutate(point_name = paste0(coord, "_", point)) |>
    select(Species, value, point_name) |>
    pivot_wider(names_from = point_name, values_from = value)
  
  bird_beak_codes |>
    select(label, is_tip, phlo) |>
    left_join(landmark_df, by = c("label" = "Species"))
  
}


## Code modified from https://github.com/mlcollyer/PACA_appendices/blob/main/SI.Appendix.2.Rmd
k.by.p <- function(x, tree) {
  
  #x <- as.matrix(x$x)
  
  Cov <- ape::vcv(tree)
  invC <- MCMCglmm::inverseA(tree, nodes = "TIPS", scale = FALSE)$Ainv
  colnames(invC) <- rownames(invC)
  invC <- invC[rownames(x), rownames(x)]
  
  #invC <- solve(Cov)
  D.mat <- RRPP:::Cov.proj(Cov, id = rownames(x))
  N <- nrow(x)
  ones <- matrix(1, N, 1) 
  a.adj <- ones %*% crossprod(ones, invC)/sum(invC)
  
  K <- function(x){
    x <- as.matrix(x)
    x.c <- x - a.adj%*%x
    MSEobs.d <- sum(x.c^2)  
    x.a <- D.mat%*%x.c
    MSE.d <- sum(x.a^2)  
    K.denom <- (sum(diag(Cov)) - N/sum(invC))/(N-1)
    (MSEobs.d/MSE.d) / K.denom
  }
  p <- NCOL(x)
  
  k.by.p <- apply(x, 2, K)
  k.to.p <- sapply(1:p, function(j) {
    xp <- as.matrix(x[, 1:j])
    K(x[,1:j])
  })
  
  list(k.by.p = k.by.p, k.to.p = k.to.p)
  
}

#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bird_beak_avonet
#' @param landmark_pf
#' @return
#' @author rdinnager
#' @export
compare_K <- function(bird_beak_avonet, landmark_pf) {

  codes <- bird_beak_avonet |>
    filter(is_tip) |>
    select(label, starts_with("latent_")) 
  
  x_lat <- as.matrix(codes[ , -1])
  rownames(x_lat) <- codes$label
  tree_lat <- pf_as_phylo(landmark_pf)
  x_lat <- scale(x_lat)
  
  lat_k <- k.by.p(x_lat, tree_lat)
  plot(lat_k$k.by.p, type = "l")
  plot(lat_k$k.to.p, type = "l")
  plot(lat_k$k.by.p, type = "l", ylim = c(0, 1), xlim = c(0, 64))
  plot(lat_k$k.to.p, type = "l", ylim = c(0, 1), xlim = c(0, 64))
  
  lat_pca <- gm.prcomp(x_lat, tree_lat, align.to.phy = FALSE)
  lat_pca_k <- k.by.p(lat_pca$x, tree_lat)
  plot(lat_pca_k$k.by.p, type = "l")
  plot(lat_pca_k$k.to.p, type = "l")
  
  lat_paca <- gm.prcomp(x_lat, tree_lat, align.to.phy = TRUE)
  lat_paca_k <- k.by.p(lat_paca$x, tree_lat)
  plot(lat_paca_k$k.by.p, type = "l")
  plot(lat_paca_k$k.to.p, type = "l")
  
  landmarks <- landmark_pf |>
    filter(is_tip) |>
    select(label, starts_with(c("x", "y", "z")))
  x_pca <- as.matrix(landmarks[ , -1])
  rownames(x_pca) <- landmarks$label
  tree_pca <- pf_as_phylo(landmark_pf)
  x_pca <- scale(x_pca)
  
  pca_k <- k.by.p(x_pca, tree_pca)
  plot(pca_k$k.by.p, type = "l")
  
  pca_pca <- gm.prcomp(x_pca, tree_pca, align.to.phy = FALSE)
  pca_paca <- gm.prcomp(x_pca, tree_pca, align.to.phy = TRUE)
  pca_phy <- gm.prcomp(x_pca, tree_pca, GLS = TRUE)
  
  pca_pca_k <- k.by.p(pca_pca$x, tree_pca)
  plot(pca_pca_k$k.by.p, type = "l")
  plot(pca_pca_k$k.to.p, type = "l")
  plot(pca_pca_k$k.by.p, type = "l", ylim = c(0, 1), xlim = c(0, 64))
  plot(pca_pca_k$k.to.p, type = "l", ylim = c(0, 1), xlim = c(0, 64))
  
  pca_paca_k <- k.by.p(pca_paca$x, tree_pca)
  plot(pca_paca_k$k.by.p, type = "l")
  plot(pca_paca_k$k.to.p, type = "l")
  plot(pca_paca_k$k.by.p, type = "l", ylim = c(0, 1), xlim = c(0, 64))
  plot(pca_paca_k$k.to.p, type = "l", ylim = c(0, 1), xlim = c(0, 64))
  
  pca_phy_k <- k.by.p(pca_phy$x, tree_pca)
  plot(pca_phy_k$k.by.p, type = "l")
  plot(pca_phy_k$k.to.p, type = "l")
  plot(pca_phy_k$k.by.p, type = "l", ylim = c(0, 1), xlim = c(0, 64))
  plot(pca_phy_k$k.to.p, type = "l", ylim = c(0, 1), xlim = c(0, 64))
  
  plot(pca_paca_k$k.by.p, type = "l", ylim = c(0, 1), xlim = c(0, 64), col = "red")
  points(lat_paca_k$k.by.p, type = "l", ylim = c(0, 1), xlim = c(0, 64), col = "green")
  points(pca_paca_k$k.to.p, type = "l", ylim = c(0, 1), xlim = c(0, 64), col = "red", lty = 2)
  points(lat_paca_k$k.to.p, type = "l", ylim = c(0, 1), xlim = c(0, 64), col = "green", lty = 2)
  
  list(x_lat = x_lat, lat_k = lat_k, 
       lat_pca = lat_pca, 
       lat_pca_k = lat_pca_k,
       lat_paca = lat_paca, 
       lat_paca_k = lat_paca_k,
       x_pca = x_pca, 
       pca_k = pca_k,
       pca_pca = pca_pca, 
       pca_paca = pca_paca, 
       pca_phy = pca_phy, 
       pca_pca_k = pca_pca_k, 
       pca_paca_k = pca_paca_k, 
       pca_phy_k = pca_phy_k)

}

#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param K_compare
#' @param landmark_pf
#' @return
#' @author rdinnager
#' @export
make_paca_phylo_plots <- function(K_compare, landmark_pf) {

  paca_lat_tree <- fortify(pf_as_phylo(landmark_pf))
  paca_lat_anc <- tibble(label = rownames(K_compare$lat_paca$anc.x),
                         `PACA DSDF 1` = K_compare$lat_paca$anc.x[, 1]) |>
    bind_rows(tibble(label = rownames(K_compare$lat_paca$x),
                     `PACA DSDF 1` = K_compare$lat_paca$x[, 1]))
  
  paca_lat_tree <- paca_lat_tree |>
    left_join(paca_lat_anc, by = "label") |>
    mutate(x0 = x[parent], x1 = x, y0 = `PACA DSDF 1`[parent], y1 = `PACA DSDF 1`)
  
  paca_lat_tree <- paca_lat_tree |>
    mutate(clade = "Clade0")
  
  lm_sp <- pf_as_sparse(pf_ones(landmark_pf))
  
  get_clade_sizes <- function(time = 100) {
    edges <- pf_epoch_info(landmark_pf$phlo, time)
    sizes <- map_int(edges$edge, ~sum(lm_sp[ , .x]))
  }
  
  pf_epoch_info(landmark_pf$phlo, 100)
  
  clade_1 <- rownames(lm_sp)[lm_sp[,"Node16"] == 1]
  
  paca_lat_tree <- paca_lat_tree |>
    mutate(clade = ifelse(label %in% clade_1, "Clade1", clade)) |>
    arrange(clade)
  
  ggplot(paca_lat_tree) +
    geom_segment(aes(x = x0, xend = x1, y = y0, yend = y1,
                     color = clade), 
                 linewidth = 0.1) +
    theme_tree2()
  
  pl_tr <- fortify(pf_as_phylo(landmark_pf))
  
  paca_lat_tree <- paca_lat_tree |>
    left_join(pl_tr |> select(label, x, y), by = "label")
  
  tr <- pf_as_phylo(landmark_pf)
  ancs <- K_compare$lat_paca$anc.x[, 1][tr$node.label]
  names(ancs) <- as.character(Ntip(tr) + 1:length(ancs))
  pacadsdf <- c(K_compare$lat_paca$x[, 1],  ancs)
  #names(pacadsdf) <- paca_lat_anc$label
  lm_tr <- pf_as_phylo(landmark_pf)
  labs <- c(lm_tr$tip.label, lm_tr$node.label)
  cols <- rep("black", length(labs))
  names(cols) <- labs
  phenogram(lm_tr, pacadsdf, ftype = "off", lwd = 0.1, colors = cols)
  
  lm_sp <- pf_as_sparse(pf_ones(landmark_pf))
  n_descs <- colSums(lm_sp)
  #plot(sort(n_descs, decreasing = TRUE)[1:100], type = "l")
  clade_1 <- rownames(lm_sp)[lm_sp[,"Node16"] == 1]
  cols[clade_1] <- "red"
  test <- phenogram(lm_tr, pacadsdf, ftype = "off", lwd = 0.1)
  
  btr <- matrix(pacadsdf, ncol = 1) %*% t(K_compare$lat_paca$rotation[ , 1])
  btr <- t((t(btr)*attr(K_compare$x_lat, "scaled:scale")) + attr(K_compare$x_lat, "scaled:center"))
  
  ggtree(paca_lat_tree |> mutate(yscale = `PACA DSDF 1`), 
         aes(colour = `PACA DSDF 1`), 
         continuous = 'colour',
         yscale = "yscale") +
    geom_tree() +
    scale_color_viridis_c() +
    labs(title = "PACA on Latent Space") +
    theme_minimal()
  
  ggtree(paca_lat_tree, aes(y = "PACA DSDF 1"), ladderize = FALSE) +
    theme_minimal()

}



