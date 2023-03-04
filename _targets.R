## Load your packages, e.g. library(targets).
source("./packages.R")
conflict_prefer("filter", "dplyr")
conflict_prefer("save.image", "imager")
conflict_prefer("select", "dplyr")

library(future)
options(future.globals.onReference = "error")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

tar_option_set(error = "continue")

## tar_plan supports drake-style targets and also tar_target()
tar_plan(

# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style
  
  tar_target(avonet_birdtree, read_csv("data/AVONET3_BirdTree.csv")),
  
  tar_target(bird_beak_avonet, bird_beak_codes %>%
               filter(is_tip) %>%
               mutate(Species3 = Scientific) %>%
               left_join(avonet_birdtree, 
                         by = "Species3") %>%
               mutate(beak_narrowness = Beak.Length_Culmen / (Beak.Width + Beak.Depth / 2),
                      beak_flatness = Beak.Width / Beak.Depth)),
  
  tar_target(codes, bird_beak_avonet %>%
               select(starts_with("latent_")) %>%
               as.matrix()),
  
  tar_target(codes_anc, bird_beak_codes %>%
               select(starts_with("latent_")) %>%
               as.matrix()),
  
  tar_target(codes_anc_st_heavy, gaussianize(codes_anc, return.tau.mat = TRUE,
                                       return.u = TRUE)),
  
  tar_target(codes_anc_st, scale(codes_anc)),
  
  tar_target(bird_beak_avonet_codes, bird_beak_codes %>%
               mutate(codes = codes_anc_st)),
  
  tar_target(ancestral_codes_mod, fibre(codes ~ bre_brownian(phlo,
                                                             standardise = FALSE,
                                                             hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))),
                                        data = bird_beak_avonet_codes,
                                        engine = "inla",
                                        ncores = 6,
                                        verbose = 2,
                                        engine_options = list(inla.mode = "experimental",
                                                              control.compute = list(waic = TRUE)#,
                                                              #control.family = list(hyper = list(hyper = list(prec = list(initial = 10,
                                                                                                #fixed = TRUE))))))
                                                              ))),
  tar_target(ancestral_codes_mod_lasso, fibre(codes ~ bre_brownian(phlo,
                                                                   standardise = FALSE),
                                              family = "mgaussian",
                                              data = bird_beak_avonet_codes,
                                              engine = "glmnet",
                                              verbose = 2,
                                              ncores = 2,
                                              engine_options = list(alpha = 1,
                                                                    nlambda = 100))),
  
  tar_target(ancestral_codes_mod_elast_0.75, fibre(codes ~ bre_brownian(phlo,
                                                                   standardise = FALSE),
                                              family = "mgaussian",
                                              data = bird_beak_avonet_codes,
                                              engine = "glmnet",
                                              verbose = 2,
                                              ncores = 2,
                                              engine_options = list(alpha = 0.75,
                                                                    nlambda = 100))),
  
  tar_target(ancestral_codes_mod_elast_0.5, fibre(codes ~ bre_brownian(phlo,
                                                                   standardise = FALSE),
                                              family = "mgaussian",
                                              data = bird_beak_avonet_codes,
                                              engine = "glmnet",
                                              verbose = 2,
                                              ncores = 2,
                                              engine_options = list(alpha = 0.5,
                                                                    nlambda = 100))),
  
  tar_target(ancestral_codes_mod_elast_0.25, fibre(codes ~ bre_brownian(phlo,
                                                                   standardise = FALSE),
                                              family = "mgaussian",
                                              data = bird_beak_avonet_codes,
                                              engine = "glmnet",
                                              verbose = 2,
                                              ncores = 2,
                                              engine_options = list(alpha = 0.25,
                                                                    nlambda = 100))),
  
    tar_target(ancestral_codes_mod_ridge, fibre(codes ~ bre_brownian(phlo,
                                                                     standardise = FALSE),
                                              family = "mgaussian",
                                              data = bird_beak_avonet_codes,
                                              engine = "glmnet",
                                              verbose = 2,
                                              ncores = 2,
                                              engine_options = list(alpha = 0.01,
                                                                    nlambda = 100))),
  
  tar_target(ancestral_codes_mod_ridge_ungrouped, fibre(codes ~ bre_brownian(phlo,
                                                                     standardise = FALSE),
                                              family = "mgaussian",
                                              data = bird_beak_avonet_codes,
                                              engine = "glmnet",
                                              verbose = 2,
                                              ncores = 2,
                                              engine_options = list(alpha = 0.01,
                                                                    nlambda = 100))),
  
  tar_target(ancestral_codes_mod_2nd, fibre(codes ~ bre_second_order(phlo,
                                                                      hyper = list(prec = list(prior = "pc.prec", param = c(5, 0.1)))),
             data = bird_beak_avonet_codes,
             engine = "inla",
             engine_options = list(num.threads = 7,
                                   inla.mode = "experimental",
                                   verbose = TRUE,
                                   control.compute = list(waic = TRUE)))),
  
  tar_target(norms, sqrt(rowSums(codes ^ 2))),
  
  tar_target(beak_vec_dat, bird_beak_avonet %>%
               select(beak_narrowness, beak_flatness,
                      starts_with("latent_")) %>%
               mutate(beak_narrowness = (beak_narrowness - mean(beak_narrowness)) / sd(beak_narrowness),
                      beak_flatness = (beak_flatness - mean(beak_flatness)) / sd(beak_flatness))),
  
  tar_target(beak_vector_narrowness, lm.fit(x = beak_vec_dat %>% select(starts_with("latent_")) %>%
                                               as.matrix(),
                                             y = beak_vec_dat %>% pull(beak_narrowness))),
  
  tar_target(beak_vector_flatness, lm.fit(x = beak_vec_dat %>% select(starts_with("latent_")) %>%
                                               as.matrix(),
                                             y = beak_vec_dat %>% pull(beak_flatness))),
  
  tar_target(starting_codes, random_code(10)),
  
  tar_target(narrowness_codes, round(seq(-0.01/3, 0.01/5, length.out = 9), 5) %>%
               map(~ t(t(starting_codes) + .x * beak_vector_narrowness$coefficients))),
  
  tar_target(flatness_codes, seq(-0.01/5, 0.01/5, length.out = 9) %>%
               map(~ t(t(starting_codes) + .x * beak_vector_flatness$coefficients))),
  
  tar_target(sdfnet, load_model("bird_beaks"),
             format = "torch"),
  
  tar_target(pal_10, c(wes_palettes$Darjeeling1, wes_palettes$FantasticFox1) %>%
               col2rgb() %>%
               array_branch(2L)),
  
  tar_target(narrowness_images, narrowness_codes[[1]] %>%
               array_branch(1L) %>%
               map2(pal_10, ~ sdfnet$cuda()$render_image(latent_code = torch_tensor(matrix(.x, nrow = 1), 
                                                                                    requires_grad = FALSE), 
                                                         camera_position = get_camera_position(-3, 125, 200),
                                                         resolution = 600,
                                                         color = c(.y[1] / 255,
                                                                   .y[2] / 255,
                                                                   .y[3] / 255),
                                                         ssaa = 1,
                                                         cuda = TRUE)),
             pattern = map(narrowness_codes),
             iteration = "list"),
  
  tar_target(flatness_images, flatness_codes[[1]] %>%
               array_branch(1L) %>%
               map2(pal_10, ~ sdfnet$cuda()$render_image(latent_code = torch_tensor(matrix(.x, nrow = 1), 
                                                                                    requires_grad = FALSE), 
                                                         camera_position = get_camera_position(-3, 160, 200),
                                                         resolution = 600,
                                                         color = c(.y[1] / 255,
                                                                   .y[2] / 255,
                                                                   .y[3] / 255),
                                                         ssaa = 1,
                                                         cuda = TRUE)),
             pattern = map(flatness_codes),
             iteration = "list"),
  
  tar_target(starting_beak_images, starting_codes %>%
               array_branch(1L) %>%
               map(~ sdfnet$cuda()$render_image(latent_code = torch_tensor(matrix(.x, nrow = 1), 
                                                                    requires_grad = FALSE), 
                                         camera_position = get_camera_position(-3, 125, 200),
                                         resolution = 400,
                                         ssaa = 1,
                                         cuda = TRUE))),
  
  tar_target(narrowness_image, assemble_images(narrowness_images,
                                               real_col = 6)),
  
  tar_target(flatness_image, assemble_images(flatness_images,
                                             real_col = 5)),
  
  tar_target(narrowness_plot, {save.image(narrowness_image,
                                         ff <- "figures/narrowness_plot.png"); ff}),
  
  tar_target(flatness_plot, {save.image(flatness_image,
                                         ff <- "figures/flatness_plot.png"); ff}),
  
  tar_target(trophic_niche_dat, bird_beak_avonet %>%
               select(Trophic.Niche,
                      starts_with("latent_"))),
  
  tar_target(trophic_niche_split, initial_split(trophic_niche_dat,
                                                0.8,
                                                strata = Trophic.Niche)),
  
  tar_target(trophic_niche_dat_train, training(trophic_niche_split)),
  
  tar_target(trophic_niche_cv, mc_cv(trophic_niche_dat,
                                    0.8,
                                    5,
                                    strata = Trophic.Niche)),
  
  tar_target(trophic_niche_cv2, vfold_cv(trophic_niche_dat_train,
                                    5,
                                    strata = Trophic.Niche)),
  
  tar_target(trophic_niche_RF, fit_RF_trophic(bird_beak_avonet,
                                              trophic_niche_dat,
                                              trophic_niche_cv)),
  
  tar_target(trophic_niche_RF2, fit_RF_trophic2(bird_beak_avonet,
                                                trophic_niche_dat_train,
                                                trophic_niche_cv2)),
  
  tar_target(trophic_niche_RF_tuned, tune_RF_trophic(trophic_niche_RF, 
                                                     trophic_niche_cv, 
                                                     trophic_niche_dat)),
  
  tar_target(trophic_niche_RF_tuned2, tune_RF_trophic(trophic_niche_RF2, 
                                                      trophic_niche_cv2, 
                                                      trophic_niche_dat_train)),
  
  tar_target(trophic_niche_RF_final, trophic_niche_RF$wf %>%
               finalize_workflow(trophic_niche_RF_tuned %>% select_best()) %>%
               fit(trophic_niche_dat)),
  
  tar_target(trophic_niche_RF_final_fit, trophic_niche_RF$wf %>%
               finalize_workflow(trophic_niche_RF_tuned %>% select_best()) %>%
               last_fit(trophic_niche_split)),
  
  tar_target(trophic_niche_RF_final2, trophic_niche_RF2$wf %>%
               finalize_workflow(trophic_niche_RF_tuned2 %>% select_best()) %>%
               fit(trophic_niche_dat_train)),
  
  tar_target(trophic_niche_RF_final_fit2, trophic_niche_RF2$wf %>%
               finalize_workflow(trophic_niche_RF_tuned2 %>% select_best()) %>%
               last_fit(trophic_niche_split)),
  
  tar_target(trophic_niche_preds, augment(trophic_niche_RF_final,
                                          trophic_niche_dat)),
  
  tar_target(trophic_niche_preds2, augment(trophic_niche_RF_final2,
                                           trophic_niche_dat)),
  
  tar_target(trophic_niche_preds_test2, augment(trophic_niche_RF_final_fit2$.workflow[[1]],
                                                testing(trophic_niche_split))),
  
  tar_target(trophic_niche_conf_mat2, conf_mat(trophic_niche_preds2 %>%
                                                mutate(Trophic.Niche = as.factor(Trophic.Niche)),
                                              Trophic.Niche,
                                              .pred_class)),
  
  tar_target(trophic_niche_conf_mat_test2, conf_mat(trophic_niche_preds_test2 %>%
                                                mutate(Trophic.Niche = as.factor(Trophic.Niche)),
                                              Trophic.Niche,
                                              .pred_class)),
  
  tar_target(code_sample, get_code_sample(codes, 100000)),
  
  tar_target(code_sample_df, as.data.frame(unclass(code_sample)) %>%
               setNames(colnames(trophic_niche_dat_train)[-1])),
  
  tar_target(code_sample_trophic_pred, predict(trophic_niche_RF_final2,
                                               new_data = code_sample_df,
                                               type = "prob")),
  
  tar_target(code_sample_classify, code_sample_trophic_pred %>%
               rowwise() %>%
               mutate(max =  max(c_across())) %>%
               mutate(across(everything(), ~ .x / max)) %>%
               mutate(w = list(which(c_across(-max) > 0.75)))),
  
  tar_target(code_sample_classes, code_sample_classify %>%
               mutate(classes = list(colnames(code_sample_classify)[w])) %>%
               mutate(class_count = length(classes)) %>%
               select(classes, class_count)),
  
  tar_target(code_sample_w_classes, code_sample_df %>%
               bind_cols(code_sample_classes) %>%
               filter(class_count < 3)),
  
  tar_target(sure_classes, code_sample_w_classes %>%
               filter(class_count == 1) %>%
               mutate(classes = unlist(classes))),
  
  tar_target(two_classes, code_sample_w_classes %>%
               filter(class_count == 2) %>%
               unnest_wider(classes,
                            names_sep = "_",
                            simplify = TRUE)),
  
  tar_target(pred_rates, ancestral_codes_mod$random$phlo %>%
               separate(ID, c("var", "edge"), ":") %>%
               select(var, edge, mean) %>%
               pivot_wider(names_from = var, values_from = mean)),
  
  tar_target(rate_trajectories, pred_rates %>%
               mutate(across(-edge, ~ .x / sqrt(sum(.x^2))))),
  
  tar_target(try_ks, c(2:18)),
  
  tar_target(rate_cluster, movMF(rate_trajectories %>%
                                   select(-edge) %>%
                                   as.matrix(), 
                                 try_ks),
             pattern = map(try_ks),
             iteration = "list"),
  
  tar_target(best_k, which.min(sapply(rate_cluster, BIC))),
  
  tar_target(best_rate_cluster, rate_cluster[[best_k]]),
  
  tar_target(rate_traj_concs, sqrt(rowSums(best_rate_cluster$theta^2))),
  
  tar_target(rate_traj_means, t(t(best_rate_cluster$theta) / rate_traj_concs)),
  
  # tar_target(bird_samp_frac, sum(avonet$is_tip) / nrow(avonet_birdtree)),
  
  # tar_target(avonet_clads_model, fit_ClaDS(pf_as_phylo(avonet),
  #                                          bird_samp_frac)),
  
  #tar_target(narrowness_plots, ),
  
  tar_target(beak_dat, bird_beak_avonet %>%
               select(Species3, starts_with("Beak.")) %>%
               drop_na()),
  
  tar_target(beak_pca, princomp(as.matrix(beak_dat %>% select(-Species3)))),
  
  tar_quarto(report, "doc/report.qmd")

)
