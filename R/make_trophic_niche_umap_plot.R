#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_umap_2
make_trophic_niche_umap_plot <- function(pl_dat,
                                         niche_pal) {

  ggplot(pl_dat, aes(`UMAP Axis 1`, `UMAP Axis 2`)) +
    geom_point(aes(colour = `Trophic.Niche`)) +
    coord_equal() +
    scale_colour_manual(values = unname(niche_pal)) +
    theme_minimal()
    

}
