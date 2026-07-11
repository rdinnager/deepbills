#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param K_compare
#' @return
#' @author rdinnager
#' @export
make_K_plots <- function(K_compare) {

  lat_k_df <- tibble(latent_axis = 1:64, K = K_compare$lat_k$k.by.p)
  p1 <- ggplot(lat_k_df, aes(latent_axis, K)) +
    geom_smooth(colour = "black") +
    geom_bar(stat = "identity") +
    xlab("Latent Axis") +
    ylab("") +
    ylim(0, 0.5) +
    ggtitle("DeepSDF") +
    theme_minimal()
  
  pca_k_df <- tibble(pc_axis = 1:64, K = K_compare$pca_pca_k$k.by.p[1:64])
  p2 <- ggplot(pca_k_df, aes(pc_axis, K)) +
    geom_smooth(colour = "black") +
    geom_bar(stat = "identity") +
    xlab("PC Axis") +
    ylab("Blomberg's K") +
    ylim(0, 0.5) +
    ggtitle("PCA on Landmarks") +
    theme_minimal()
  
  pca_phy_k_df <- tibble(pc_axis = 1:64, K = K_compare$pca_phy_k$k.by.p[1:64])
  p3 <- ggplot(pca_phy_k_df, aes(pc_axis, K)) +
    geom_smooth(colour = "black") +
    geom_bar(stat = "identity") +
    xlab("PC Axis") +
    ylab("") +
    ylim(0, 0.5) +
    ggtitle("Phylogenetic PCA\non Landmarks") +
    theme_minimal()
  
  ragg::agg_png("figures/K_plots.png", width = 1380, height = 800,
                scaling = 2.5)
  
  pp <- p2 + p3 + p1 
  plot(pp)
  
  dev.off()
  
  lat_paca_df <- tibble(type = "DeepSDF-PACA", axis = 1:64, K = K_compare$lat_paca_k$k.by.p[1:64])
  pca_paca_df <- tibble(type = "Landmark-PACA", axis = 1:64, K = K_compare$pca_paca_k$k.by.p[1:64])
  paca_df <- bind_rows(lat_paca_df, pca_paca_df)
  
  pp2 <- ggplot(paca_df, aes(axis, K, fill = type)) +
    geom_col(position = "dodge") +
    #geom_smooth(aes(fill = type, colour = type)) +
    xlab("Axis") +
    ylab("Blomberg's K") +
    scale_wrap(n = 2) +
    theme_minimal() +
    ggtitle('Phylogenetically Aligned Components Analysis (PACA)') +
    theme(legend.position = "bottom",
          legend.title = element_blank())
  
  ragg::agg_png("figures/K_plot_PACA.png", width = 800, height = 1000,
                scaling = 2.5)
  
  print(pp2)
  
  dev.off()
  
  ragg::agg_png("figures/K_plots_combined.png", width = 1400, height = 1800,
                scaling = 3.25)
  pp + pp2 + plot_layout(nrow = 2, heights = c(0.5, 0.5)) + 
    plot_annotation(tag_levels = "a")
  
  dev.off()

  list(pp, pp2)
}
