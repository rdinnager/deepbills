#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_RF
tune_RF_trophic <- function(trophic_niche_RF_ai64_aug, 
                            trophic_niche_cv_ai64_aug_edited,
                            trophic_niche_dat_ai64_aug) {

  registerDoParallel(6)
  
  tuned <- trophic_niche_RF_ai64_aug$wf %>%
               tune_bayes(resamples = trophic_niche_cv_ai64_aug_edited,
                          initial = trophic_niche_RF_ai64_aug$tune,
                          param_info = extract_parameter_set_dials(trophic_niche_RF_ai64_aug$wf) %>%
                            finalize(trophic_niche_dat_ai64_aug),
                          metrics = metric_set(accuracy, roc_auc),
                          control = control_bayes(verbose = TRUE,
                                                  verbose_iter = TRUE))
  
  stopImplicitCluster()

  tuned
  
}
