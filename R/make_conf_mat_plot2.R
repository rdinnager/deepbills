#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_conf_mat_GLM_test_pca
#' @return
#' @author rdinnager
#' @export
make_conf_mat_plot2 <- function(trophic_niche_conf_mat_RF_test_ai64, file_name = "figures/test_ai64_conf_mat.png") {

  tt<-as.matrix(trophic_niche_conf_mat_RF_test_ai64$table)
  tt <- t(t(tt)/apply(tt, 2, sum))
  
  ragg::agg_png(file_name, width = 800, height = 640, scaling = 2)
  corrplot::corrplot(tt, "square", is.corr = FALSE)
  dev.off()
  
  file_name
  
}
