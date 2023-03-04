#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_RF
tune_RF_trophic <- function(trophic_niche_RF, trophic_niche_cv, trophic_niche_dat) {

  registerDoParallel(6)
  
  tuned <- trophic_niche_RF$wf %>%
               tune_bayes(resamples = trophic_niche_cv,
                          initial = trophic_niche_RF$tune,
                          param_info = extract_parameter_set_dials(trophic_niche_RF$wf) %>%
                            finalize(trophic_niche_dat))
  
  stopImplicitCluster()

  tuned
  
}
