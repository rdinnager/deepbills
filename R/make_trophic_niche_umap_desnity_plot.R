#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param pl_dat_unsupervised
#' @param niche_pal
make_trophic_niche_umap_density_plot <- function(pl_dat_unsupervised, niche_pal) {

  ggplot(pl_dat_unsupervised, aes(`UMAP Axis 1`, `UMAP Axis 2`)) +
    geom_point(aes(colour = `Trophic.Niche`), alpha = 0.1, show.legend = FALSE) +
    geom_density2d(aes(colour = `Trophic.Niche`), contour_var = "ndensity",
                   show.legend = FALSE) +
    coord_equal() +
    facet_wrap(vars(Trophic.Niche), nrow = 3) +
    scale_colour_manual(values = unname(niche_pal)) +
    theme_minimal()
    
  # ggplot(pl_dat_unsupervised, aes(`UMAP Axis 1`, `UMAP Axis 2`)) +
  #   geom_point(aes(colour = `Trophic.Niche`), alpha = 0.1) +
  #   geom_density2d(aes(colour = `Trophic.Niche`), contour_var = "ndensity") +
  #   theme_minimal()

}
