#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_conf_mat_GLM_test_pca
#' @return
#' @author rdinnager
#' @export
make_conf_mat_plot2 <- function(trophic_niche_conf_mat_RF_test_ai64) {

  tt<-as.matrix(trophic_niche_conf_mat_RF_test_pca$table)
  tt <- t(t(tt)/apply(tt, 2, sum))
  
  corrplot::corrplot(tt, "square", is.corr = FALSE)
  
}
