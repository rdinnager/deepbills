#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_RF_preds_ai64
#' @param trophic_niche_RF_preds_pca
#' @param trophic_niche_RF_preds_test_ai64
#' @param trophic_niche_RF_preds_test_pca
#' @param trophic_niche_GLM_preds_ai64
#' @param trophic_niche_GLM_preds_pca
#' @param trophic_niche_GLM_preds_test_ai64
#' @param trophic_niche_GLM_preds_test_pca
#' @return
#' @author rdinnager
#' @export
get_niche_metrics <- function(metrics,
                              trophic_levs2,
                              trophic_niche_RF_preds_test_ai64,
                              trophic_niche_RF_preds_test_pca,
                              trophic_niche_GLM_preds_test_ai64,
                              trophic_niche_GLM_preds_test_pca) {

  
  RF_ai64 <- metrics(trophic_niche_RF_preds_test_ai64 |>
                       mutate(Trophic.Niche = as.factor(Trophic.Niche)) |>
                       group_by(Trophic.Niche),
                     starts_with(".pred"), -.pred_class,
                     truth = Trophic.Niche,
                     estimate = .pred_class) 
  
  RF_pca <- metrics(trophic_niche_RF_preds_test_pca |>
                       mutate(Trophic.Niche = as.factor(Trophic.Niche)) |>
                       group_by(Trophic.Niche),
                     starts_with(".pred"), -.pred_class,
                     truth = Trophic.Niche,
                     estimate = .pred_class) 
  
  mean(RF_ai64 |> filter(Trophic.Niche != "Omnivore") |> pull(.estimate))
  mean(RF_pca |> filter(Trophic.Niche != "Omnivore") |> pull(.estimate))
  
  GLM_ai64 <- metrics(trophic_niche_GLM_preds_test_ai64 |>
                       mutate(Trophic.Niche = as.factor(Trophic.Niche)) |>
                       group_by(Trophic.Niche),
                     starts_with(".pred"), -.pred_class,
                     truth = Trophic.Niche,
                     estimate = .pred_class)
  
  GLM_pca <- metrics(trophic_niche_GLM_preds_test_pca |>
                       mutate(Trophic.Niche = as.factor(Trophic.Niche)) |>
                       group_by(Trophic.Niche),
                     starts_with(".pred"), -.pred_class,
                     truth = Trophic.Niche,
                     estimate = .pred_class)
  
    
  GLM_train_ai64 <- metrics(trophic_niche_GLM_preds_ai64 |>
                       mutate(Trophic.Niche = as.factor(Trophic.Niche)) |>
                       group_by(Trophic.Niche),
                     starts_with(".pred"), -.pred_class,
                     truth = Trophic.Niche,
                     estimate = .pred_class)
  
  GLM_train_pca <- metrics(trophic_niche_GLM_preds_pca |>
                       mutate(Trophic.Niche = as.factor(Trophic.Niche)) |>
                       group_by(Trophic.Niche),
                     starts_with(".pred"), -.pred_class,
                     truth = Trophic.Niche,
                     estimate = .pred_class)
  
}
