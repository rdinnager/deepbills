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


fit_RF_trophic3 <- function(trophic_niche_dat_train_ai,
                            trophic_niche_cv_ai,
                            ncodes = 15,
                            code_names = "latent_") {


  trophic_recipe <- recipe(trophic_niche_dat_train_ai,
                           vars = colnames(trophic_niche_dat_train_ai %>%
                                             select(Trophic.Niche,
                                                    starts_with(code_names),
                                                    weights)),
                           roles = c("outcome", rep("predictor", ncodes), "case_weights")) 

  trophic_mod <- rand_forest(mtry = tune(),
                            trees = tune(),
                            min_n = tune()) %>%
    set_mode("classification") %>%
    set_engine("ranger")

  trophic_wf <- workflow() %>%
    add_model(trophic_mod) %>%
    add_recipe(trophic_recipe) %>%
    add_case_weights(weights)
  
  registerDoParallel(6)
  
  trophic_tune <- trophic_wf %>%
    tune_grid(resamples = trophic_niche_cv_ai,
              grid = 64,
              control = control_grid(verbose = TRUE))
  
  stopImplicitCluster()
  
  list(wf = trophic_wf, tune = trophic_tune)
    
}

fit_RF_trophic3_aug <- function(trophic_niche_dat_train_ai,
                            trophic_niche_cv_ai,
                            ncodes = 15,
                            code_names = "latent_") {


  trophic_recipe <- recipe(trophic_niche_dat_train_ai,
                           vars = colnames(trophic_niche_dat_train_ai %>%
                                             select(Trophic.Niche,
                                                    starts_with(code_names))),
                           roles = c("outcome", rep("predictor", ncodes))) 

  trophic_mod <- rand_forest(mtry = tune(),
                            trees = tune(),
                            min_n = tune()) %>%
    set_mode("classification") %>%
    set_engine("ranger")

  trophic_wf <- workflow() %>%
    add_model(trophic_mod) %>%
    add_recipe(trophic_recipe) 
  
  registerDoParallel(6)
  
  trophic_tune <- trophic_wf %>%
    tune_grid(resamples = trophic_niche_cv_ai,
              grid = 64,
              control = control_grid(verbose = TRUE))
  
  stopImplicitCluster()
  
  list(wf = trophic_wf, tune = trophic_tune)
    
}

fit_GLM_trophic3 <- function(trophic_niche_dat_train_ai,
                            trophic_niche_cv_ai,
                            ncodes = 15,
                            code_names = "latent_") {


  trophic_recipe <- recipe(trophic_niche_dat_train_ai,
                           vars = colnames(trophic_niche_dat_train_ai %>%
                                             select(Trophic.Niche,
                                                    starts_with(code_names),
                                                    weights)),
                           roles = c("outcome", rep("predictor", ncodes), "case_weights")) 

  trophic_mod <- multinom_reg(penalty = tune(),
                              mixture = tune()) %>%
    set_mode("classification") %>%
    set_engine("glmnet")

  trophic_wf <- workflow() %>%
    add_model(trophic_mod) %>%
    add_recipe(trophic_recipe) %>%
    add_case_weights(weights)
  
  registerDoParallel(6)
  
  trophic_tune <- trophic_wf %>%
    tune_grid(resamples = trophic_niche_cv_ai,
              grid = 64,
              control = control_grid(verbose = TRUE))
  
  stopImplicitCluster()
  
  list(wf = trophic_wf, tune = trophic_tune)
    
}

fit_GLM_trophic3_aug <- function(trophic_niche_dat_train_ai,
                                 trophic_niche_cv_ai,
                                 ncodes = 15,
                                 code_names = "latent_") {


  trophic_recipe <- recipe(trophic_niche_dat_train_ai,
                           vars = colnames(trophic_niche_dat_train_ai %>%
                                             select(Trophic.Niche,
                                                    starts_with(code_names))),
                           roles = c("outcome", rep("predictor", ncodes))) 

  trophic_mod <- multinom_reg(penalty = tune(),
                              mixture = tune()) %>%
    set_mode("classification") %>%
    set_engine("glmnet")

  trophic_wf <- workflow() %>%
    add_model(trophic_mod) %>%
    add_recipe(trophic_recipe) 
  
  registerDoParallel(6)
  
  trophic_tune <- trophic_wf %>%
    tune_grid(resamples = trophic_niche_cv_ai,
              grid = 64,
              control = control_grid(verbose = TRUE))
  
  stopImplicitCluster()
  
  list(wf = trophic_wf, tune = trophic_tune)
    
}