#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bird_beak_avonet
fit_RF_trophic <- function(bird_beak_avonet, dat, cv) {

  trophic_recipe <- recipe(dat,
                           vars = colnames(bird_beak_avonet %>%
                                             select(Trophic.Niche,
                                                    starts_with("latent_"))),
                           roles = c("outcome", rep("predictor", 64))) %>%
    step_normalize(all_predictors()) %>%
    step_smote(Trophic.Niche, over_ratio = 0.5) %>%
    step_downsample(Trophic.Niche)

  trophic_mod <- rand_forest(mtry = tune(),
                             trees = tune(),
                             min_n = tune()) %>%
    set_engine("ranger") %>%
    set_mode("classification")

  trophic_wf <- workflow() %>%
    add_model(trophic_mod) %>%
    add_recipe(trophic_recipe)
  
  registerDoParallel(6)
  
  trophic_tune <- trophic_wf %>%
    tune_grid(resamples = cv,
              grid = 64,
              control = control_grid(verbose = TRUE,
                                     save_pred = TRUE))
  
  stopImplicitCluster()
  
  list(wf = trophic_wf, tune = trophic_tune)
    
}

#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bird_beak_avonet
fit_RF_trophic2 <- function(bird_beak_avonet, dat, cv) {

  trophic_recipe <- recipe(dat,
                           vars = colnames(bird_beak_avonet %>%
                                             select(Trophic.Niche,
                                                    starts_with("latent_"))),
                           roles = c("outcome", rep("predictor", 64))) %>%
    step_normalize(all_predictors()) %>%
    step_smote(Trophic.Niche)

  trophic_mod <- rand_forest(mtry = tune(),
                             trees = tune(),
                             min_n = tune()) %>%
    set_engine("ranger") %>%
    set_mode("classification")

  trophic_wf <- workflow() %>%
    add_model(trophic_mod) %>%
    add_recipe(trophic_recipe)
  
  registerDoParallel(6)
  
  trophic_tune <- trophic_wf %>%
    tune_grid(resamples = cv,
              grid = 64,
              control = control_grid(verbose = TRUE,
                                     save_pred = TRUE))
  
  stopImplicitCluster()
  
  list(wf = trophic_wf, tune = trophic_tune)
    
}