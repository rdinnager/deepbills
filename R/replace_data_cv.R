#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_cv_ai
#' @param trophic_niche_dat_train_pc
#' @return
#' @author rdinnager
#' @export
replace_data_cv <- function(trophic_niche_cv_ai, trophic_niche_dat_train_pc) {

  trophic_niche_cv_pc <- trophic_niche_cv_ai
  for(i in 1:length(trophic_niche_cv_ai$splits)) {
    trophic_niche_cv_pc$splits[[i]]$data <- trophic_niche_dat_train_pc
  }
  trophic_niche_cv_pc

}
