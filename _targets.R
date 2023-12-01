## Load your packages, e.g. library(targets).
source("./packages.R")
conflict_prefer("filter", "dplyr")
conflict_prefer("save.image", "imager")
conflict_prefer("select", "dplyr")
conflict_prefer("%<-%", "zeallot")

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
  
  tar_target(clads_tree, read.tree("data/clads_tree.tre")),
  
  tar_target(bird_beak_avonet, bird_beak_codes %>%
               filter(is_tip) %>%
               mutate(Species3 = Scientific) %>%
               left_join(avonet_birdtree, 
                         by = "Species3") %>%
               mutate(beak_narrowness = Beak.Length_Culmen / (Beak.Width + Beak.Depth / 2),
                      beak_flatness = Beak.Width / Beak.Depth)),
  
  # tar_target(bird_beak_clads, pf_as_pf(clads_tree) %>%
  #              left_join(bird_beak_codes %>%
  #              filter(is_tip) %>%
  #              select(label, starts_with("latent_"))) %>%
  #              filter(!is.na(latent_code_1))),
  
  tar_target(clads_file, "data/clads_pf.rds", format = "file"),
  
  tar_target(div_rate_clads, read_rds(clads_file)),
  
  tar_target(bird_beak_clads, reconcile_clads(div_rate_clads,
                                              bird_beak_codes)),
  
  tar_target(codes, bird_beak_avonet %>%
               select(starts_with("latent_")) %>%
               as.matrix()),
  
  tar_target(two_stage_cvae, run_cvae(bird_beak_avonet),
             format = "torch"),
  
  tar_target(trophic_levs, bird_beak_avonet %>%
               select(Trophic.Niche) %>%
               mutate(Trophic.Niche = factor(Trophic.Niche)) %>%
               pull(Trophic.Niche) %>%
               unique()),
  
  tar_target(trophic_levs2, trophic_levs[order(trophic_levs)]),
  
  tar_target(cvae_latents, get_cvae_latents(two_stage_cvae,
                                            bird_beak_avonet)),
  
  tar_target(trophic_latent_samples, get_latent_samples(trophic_levs2,
                                                        two_stage_cvae,
                                                        niche_pal,
                                                        bird_beak_avonet),
             pattern = map(trophic_levs2, niche_pal)),
  
  tar_target(trophic_latent_samples_ims, trophic_latent_samples$image,
             pattern = map(trophic_latent_samples),
             iteration = "list"),
  
  tar_target(trophic_sample_pngs, save.image(trophic_latent_samples_ims,
                                             file.path("figures/trophic",
                                                       paste0(as.character(trophic_levs2),
                                                              ".png")),
                                             0.9),
             pattern = map(trophic_latent_samples_ims, trophic_levs2)),
  
  tar_target(trophic_sample_plot, make_trophic_sample_plot(trophic_latent_samples_ims,
                                                           trophic_levs2)),
  
  tar_target(codes_clads, bird_beak_clads %>%
               select(starts_with("latent_")) %>%
               as.matrix()),
  
  tar_target(codes_anc, bird_beak_codes %>%
               select(starts_with("latent_")) %>%
               as.matrix()),
  
  tar_target(codes_anc_clads, bird_beak_clads %>%
               select(starts_with("latent_")) %>%
               as.matrix()),
  
  tar_target(codes_anc_st_heavy, gaussianize(codes_anc, return.tau.mat = TRUE,
                                       return.u = TRUE)),
  
  tar_target(codes_anc_st_heavy_clads, gaussianize(codes_anc_clads, return.tau.mat = TRUE,
                                                   return.u = TRUE)),
  
  tar_target(codes_anc_st, scale(codes_anc)),
  
  tar_target(codes_anc_st_clads, scale(codes_anc_clads)),
  
  tar_target(bird_beak_avonet_codes, bird_beak_codes %>%
               mutate(codes = codes_anc_st)),
  
   tar_target(bird_beak_clads_codes, bird_beak_clads %>%
               mutate(codes = codes_anc_st_clads)),
  
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
  
   tar_target(ancestral_codes_mod_clads, fibre(codes ~ bre_brownian(phlo,
                                                                    standardise = FALSE,
                                                                    hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.1)))),
                                               data = bird_beak_clads_codes,
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
  
  tar_target(trophic_niche_codes_st, bird_beak_avonet %>%
               select(label, Trophic.Niche) %>%
               mutate(code = codes_anc_st)),
  
  tar_target(trophic_niche_dat, bird_beak_avonet %>%
               select(label,
                      Trophic.Niche,
                      starts_with("latent_")) %>%
               mutate(codes = scale(across(starts_with("latent_"), ~ .x))) %>%
               select(label, Trophic.Niche, codes) %>%
               group_split(Trophic.Niche),
             iteration = "list"),
  
  tar_target(trophic_niche_dat_all, bird_beak_avonet %>%
               select(label,
                      Trophic.Niche,
                      starts_with("latent_")) %>%
               mutate(codes = scale(across(starts_with("latent_"), ~ .x))) %>%
               select(label, Trophic.Niche, codes)),
  
  tar_target(trophic_niche_umap_2, umap(trophic_niche_dat_all %>% pull(codes),
                                        y = as.factor(trophic_niche_dat_all$Trophic.Niche),
                                        ret_model = TRUE)),
  
  tar_target(trophic_niche_umap_10, umap(trophic_niche_dat_all %>% pull(codes),
                                         n_components = 10,
                                         ret_model = TRUE)),
  
  tar_target(niche_dat_10, trophic_niche_umap_10$embedding %>%
               as.data.frame() %>%
               mutate(Trophic.Niche = as.factor(trophic_niche_dat_all$Trophic.Niche)) %>%
               group_by(Trophic.Niche) %>%
               group_split(),
             iteration = "list"),
  
  tar_target(trophic_hypervolumes_10, hypervolume_gaussian(niche_dat_10 %>%
                                                             select(starts_with("V")),
                                                           samples.per.point = 1000),
             pattern = map(niche_dat_10),
             iteration = "list"),
  
  tar_target(trophic_hypervolume_10, do.call(hypervolume_join, trophic_hypervolumes_10)),
  
  tar_target(trophic_hypervolume_10_df, hypervolume_to_data_frame(trophic_hypervolume_10) %>%
               slice_sample(n = 1000, weight_by = ValueAtRandomPoints)),
  
  tar_target(pl_dat, trophic_niche_umap_2$embedding %>%
               as.data.frame() %>%
               mutate(Trophic.Niche = as.factor(trophic_niche_dat_all$Trophic.Niche)) %>%
               rename(`UMAP Axis 1` = V1, `UMAP Axis 2` = V2)),
  
  tar_target(trophic_niche_umap_2_unsupervised, 
             umap(trophic_niche_dat_all %>% pull(codes),
                  n_neighbors = 25,
                  ret_model = TRUE)),
  
  tar_target(pl_dat_unsupervised, trophic_niche_umap_2_unsupervised$embedding %>%
               as.data.frame() %>%
               mutate(Trophic.Niche = as.factor(trophic_niche_dat_all$Trophic.Niche)) %>%
               rename(`UMAP Axis 1` = V1, `UMAP Axis 2` = V2)),
  
  tar_target(niche_pal, createPalette(nlevels(pl_dat$Trophic.Niche),
                                      wes_palettes$FantasticFox1[-1])),
  
  tar_target(trophic_niche_umap_plot, make_trophic_niche_umap_plot(pl_dat, niche_pal)),
  
  tar_target(trophic_niche_umap_unsupervised_plot, make_trophic_niche_umap_plot(pl_dat_unsupervised, niche_pal)),
  
  tar_target(trophic_niche_umap_unsupervised_density_plot, make_trophic_niche_umap_density_plot(pl_dat_unsupervised, niche_pal)),
  
  tar_target(trophic_niche_umap_unsupervised_plot_2, trophic_niche_umap_unsupervised_plot +
               trophic_niche_umap_unsupervised_density_plot +
               plot_layout(nrow = 2, heights = c(1, 1), guides = "collect") +
               plot_annotation(tag_levels = "a", tag_suffix = ")")),
  
  tar_target(trophic_hypervolumes, hypervolume_gaussian(trophic_niche_dat %>%
                                                          pull(codes)),
             pattern = map(trophic_niche_dat),
             iteration = "list"),
  
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
  
  tar_target(trophic_niche_preds_train2, augment(trophic_niche_RF_final_fit2$.workflow[[1]],
                                                 training(trophic_niche_split))),
  
  tar_target(metrics, metric_set(accuracy, bal_accuracy, roc_auc)),
  
  tar_target(trophic_niche_conf_mat2, conf_mat(trophic_niche_preds2 %>%
                                                mutate(Trophic.Niche = as.factor(Trophic.Niche)),
                                              Trophic.Niche,
                                              .pred_class)),
  
  tar_target(trophic_niche_conf_mat_test2, conf_mat(trophic_niche_preds_test2 %>%
                                                mutate(Trophic.Niche = as.factor(Trophic.Niche)),
                                              Trophic.Niche,
                                              .pred_class)),
  
  tar_target(trophic_niche_conf_mat_train2, conf_mat(trophic_niche_preds_train2 %>%
                                                mutate(Trophic.Niche = as.factor(Trophic.Niche)),
                                              Trophic.Niche,
                                              .pred_class)),
  
  tar_target(trophic_niche_metrics_train2, metrics(trophic_niche_preds_train2 %>%
                                                     mutate(Trophic.Niche = as.factor(Trophic.Niche)),
                                                   Trophic.Niche,
                                                   `.pred_Aquatic predator`:`.pred_Vertivore`,
                                                   estimate = .pred_class)),
  
  tar_target(trophic_niche_metrics_test2, metrics(trophic_niche_preds_test2 %>%
                                                     mutate(Trophic.Niche = as.factor(Trophic.Niche)),
                                                   Trophic.Niche,
                                                   `.pred_Aquatic predator`:`.pred_Vertivore`,
                                                   estimate = .pred_class)),
  
  tar_target(trophic_niche_conf_mat_train_plot, make_conf_mat_plot(trophic_niche_conf_mat2)),
  
  tar_target(trophic_niche_conf_mat_test_plot, make_conf_mat_plot(trophic_niche_conf_mat_test2)),
  
  tar_target(trophic_niche_conf_mat_train_plot_pdf, ggsave("figures/trophic_niche_confusion_train.pdf",
                                                           trophic_niche_conf_mat_train_plot,
                                                           height = 7, width = 7),
             format = "file"),
  
  tar_target(trophic_niche_conf_mat_test_plot_pdf, ggsave("figures/trophic_niche_confusion_test.pdf",
                                                           trophic_niche_conf_mat_test_plot,
                                                           height = 7, width = 7),
             format = "file"),
  
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
  
  tar_target(edge_lengths_clads, tibble(edge = pf_edge_names(bird_beak_clads_codes$phlo),
                                  len = pf_mean_edge_features(bird_beak_clads_codes$phlo))),
  
  tar_target(edge_lengths, tibble(edge = pf_edge_names(bird_beak_avonet_codes$phlo),
                                  len = pf_mean_edge_features(bird_beak_avonet_codes$phlo))),
  
  tar_target(pred_rates, ancestral_codes_mod$random$phlo %>%
               separate(ID, c("var", "edge"), ":") %>%
               left_join(edge_lengths) %>%
               mutate(mean = ifelse(mean < 0.00001, 0, mean)) %>%
               mutate(mean = mean * (1 / sqrt(len))) %>%
               select(var, edge, mean) %>%
               pivot_wider(names_from = var, values_from = mean)),
  
  tar_target(pred_rates_clads, ancestral_codes_mod_clads$random$phlo %>%
               separate(ID, c("var", "edge"), ":") %>%
               left_join(edge_lengths_clads) %>%
               mutate(mean = ifelse(mean < 0.00001, 0, mean)) %>%
               mutate(mean = mean * (1 / sqrt(len))) %>%
               select(var, edge, mean) %>%
               pivot_wider(names_from = var, values_from = mean)),
  
  tar_target(rate_w_norm, pred_rates %>%
               rowwise() %>%
               mutate(norm = sqrt(sum(c_across(-edge)^2)))),
  
  tar_target(rate_w_norm_clads, pred_rates_clads %>%
               rowwise() %>%
               mutate(norm = sqrt(sum(c_across(-edge)^2)))),
  
  tar_target(rate_trajectories, rate_w_norm %>%
               ungroup() %>%
               mutate(across(c(-edge, -norm), ~ .x / norm))),
  
  tar_target(rate_trajectories_clads, rate_w_norm_clads %>%
               ungroup() %>%
               mutate(across(c(-edge, -norm), ~ .x / norm))),
  
  tar_target(try_ks, c(2:18)),
  
  tar_target(rate_cluster, movMF(rate_trajectories %>%
                                   select(-edge, - norm) %>%
                                   as.matrix(), 
                                 try_ks),
             pattern = map(try_ks),
             iteration = "list"),
  
  tar_target(rate_cluster_clads, movMF(rate_trajectories_clads %>%
                                   select(-edge, - norm) %>%
                                   as.matrix(), 
                                 try_ks),
             pattern = map(try_ks),
             iteration = "list"),
  
  tar_target(best_k, which.min(sapply(rate_cluster, BIC))),
  
  tar_target(best_k_clads, which.min(sapply(rate_cluster_clads, BIC))),
  
  tar_target(best_rate_cluster, rate_cluster[[best_k]]),
  
  tar_target(best_rate_cluster_clads, rate_cluster_clads[[best_k_clads]]),
  
  tar_target(rate_traj_concs, sqrt(rowSums(best_rate_cluster$theta^2))),
  
  tar_target(rate_traj_concs_clads, sqrt(rowSums(best_rate_cluster_clads$theta^2))),
  
  tar_target(rate_traj_means_clads, best_rate_cluster_clads$theta / rate_traj_concs_clads),
  
  tar_target(rate_probs, dmovMF(rate_trajectories %>%
                                   select(-edge, - norm) %>%
                                   as.matrix(), 
                                best_rate_cluster$theta,
                                best_rate_cluster$alpha,
                                log = TRUE)),
  
  tar_target(rate_probs_clads, dmovMF(rate_trajectories_clads %>%
                                   select(-edge, - norm) %>%
                                   as.matrix(), 
                                best_rate_cluster_clads$theta,
                                best_rate_cluster_clads$alpha,
                                log = TRUE)),
  
  tar_target(rate_probs_clads_all, dmovMF(rate_trajectories_clads %>%
                                            select(-edge, - norm) %>%
                                            as.matrix(), 
                                          best_rate_cluster$theta,
                                          best_rate_cluster$alpha,
                                          log = TRUE)),
  
  tar_target(rates_w_summ, rate_w_norm %>%
               ungroup() %>%
               mutate(traj_dens = rate_probs,
                      traj_class = predict(best_rate_cluster))),
  
  tar_target(rates_w_summ_clads_all, rate_w_norm_clads %>%
               ungroup() %>%
               mutate(traj_dens = rate_probs_clads_all,
                      traj_class = predict(best_rate_cluster, rate_trajectories_clads %>%
                                             select(-edge, - norm) %>%
                                             as.matrix())) %>%
               left_join(bird_beak_clads %>%
                           select(edge = label, new_div_rate))),
  
  tar_target(rates_w_summ_clads, rate_w_norm_clads %>%
               ungroup() %>%
               mutate(traj_dens = rate_probs_clads,
                      traj_class = predict(best_rate_cluster_clads)) %>%
               left_join(bird_beak_clads %>%
                           select(edge = label, new_div_rate))),
  
  tar_target(beak_codes_aces, bind_cols(bird_beak_clads_codes %>%
                                          select(label, is_tip, phlo),
                                        map_dfc(predict(ancestral_codes_mod_clads),
                                                ~ .x$.pred_mean))),
  
  tar_target(phylo_w_rates_clads, bird_beak_clads_codes %>%
               left_join(rates_w_summ_clads %>%
                           select(label = edge, norm, traj_dens, traj_class, new_div_rate)) %>%
               left_join(edge_lengths_clads, by = c("label" = "edge")) %>%
               mutate(total_change = norm * len)),
  
  tar_target(phylo_w_rates_clads_all, bird_beak_clads_codes %>%
               left_join(rates_w_summ_clads_all %>%
                           select(label = edge, norm, traj_dens, traj_class, new_div_rate)) %>%
               left_join(edge_lengths_clads, by = c("label" = "edge")) %>%
               mutate(total_change = norm * len)),
  
  tar_target(phylo_w_rates, bird_beak_avonet_codes %>%
               left_join(rates_w_summ %>%
                           select(label = edge, norm, traj_dens, traj_class)) %>%
               left_join(edge_lengths, by = c("label" = "edge")) %>%
               mutate(total_change = norm * len)),
  
  tar_target(biggest_traj, phylo_w_rates %>%
               group_by(traj_class) %>%
               filter(total_change > 0.25) %>%
               slice_max(norm, n = 5)),
  
  tar_target(biggest_traj_clads, phylo_w_rates_clads %>%
               group_by(traj_class) %>%
               filter(total_change > 0.25) %>%
               slice_max(norm, n = 4)),
  
  tar_target(biggest_traj_clads_all, phylo_w_rates_clads_all %>%
               group_by(traj_class) %>%
               filter(total_change > 0.25) %>%
               slice_max(norm, n = 5)),
  
  tar_target(clust_pal, createPalette(length(best_rate_cluster_clads$alpha),
                                      wes_palettes$FantasticFox1)),
  
  tar_target(main_tree_plot_examples, make_main_tree_plot(phylo_w_rates_clads, biggest_traj_clads,
                                                          clust_pal)),
  
  tar_target(main_tree_plot_examples_pdf, ggsave("figures/main_tree_plot_examples.pdf", main_tree_plot,
                                                 width = 11, height = 8),
             format = "file"),
  
  tar_target(main_tree_plot, make_main_tree_plot(phylo_w_rates_clads, biggest_traj_clads,
                                                 clust_pal, add_examples = FALSE)),
  
  tar_target(main_tree_plot_pdf, ggsave("figures/main_tree_plot.pdf", main_tree_plot,
                                        width = 11, height = 8),
             format = "file"),
  
  tar_target(trajectory_codes, make_trajectory_codes(biggest_traj_clads, beak_codes_aces,
                                                     codes_anc_st_clads)),
  
  tar_target(trajectory_pal, clust_pal[as.numeric(as.character(biggest_traj_clads$traj_class))] %>%
               col2rgb() %>%
               array_branch(2)),
  
  tar_target(trajectory_images, trajectory_codes[[1]] %>%
               array_branch(1L) %>%
               map(~ sdfnet$cuda()$render_image(latent_code = torch_tensor(matrix(.x, nrow = 1),
                                                                           requires_grad = FALSE), 
                                                 camera_position = get_camera_position(-3, 150, 200),
                                                 resolution = 600,
                                                 color = c(trajectory_pal[[1]][1] / 255,
                                                           trajectory_pal[[1]][2] / 255,
                                                           trajectory_pal[[1]][3] / 255),
                                                 ssaa = 1,
                                                 cuda = TRUE)),
             pattern = map(trajectory_codes, trajectory_pal),
             iteration = "list"),
  
  tar_target(clust_mean_vects, array_branch(rate_traj_means_clads, 1),
             iteration = "list"),
  
  tar_target(cluster_codes, round(seq(0, 0.25, length.out = 9), 5) %>%
               map(~ t(t(starting_codes) + .x * clust_mean_vects)),
             pattern = map(clust_mean_vects),
             iteration = "list"),
  
  tar_target(cluster_images, cluster_codes %>%
               map(function(clust) {
                 clust %>%
                  array_branch(1L) %>%
                   map2(pal_10, ~ sdfnet$cuda()$render_image(latent_code = torch_tensor(matrix(.x, nrow = 1), 
                                                                                        requires_grad = FALSE), 
                                                             camera_position = get_camera_position(-3, 150, 200),
                                                             resolution = 600,
                                                             color = c(.y[1] / 255,
                                                                       .y[2] / 255,
                                                                       .y[3] / 255),
                                                             ssaa = 1,
                                                             cuda = TRUE))  
                 }),
             pattern = map(cluster_codes),
             iteration = "list"),
  
  tar_target(bird_anc_recon, sdfnet$cuda()$render_image(latent_code = 
                                                          torch_tensor(matrix((ancestral_codes_mod_clads$fixed$mean *
                                                                                attr(codes_anc_st_clads, "scaled:scale")) +
                                                                                attr(codes_anc_st_clads, "scaled:center"),
                                                                              nrow = 1),
                                                                       requires_grad = FALSE), 
                                                        camera_position = get_camera_position(-3, 125, 200),
                                                        resolution = 800,
                                                        ssaa = 1,
                                                        cuda = TRUE)),
  
  tar_target(bird_zero_recon, sdfnet$cuda()$render_image(latent_code = 
                                                          torch_tensor(matrix(0,
                                                                              nrow = 1,
                                                                              ncol = 64),
                                                                       requires_grad = FALSE), 
                                                        camera_position = get_camera_position(-3, 125, 200),
                                                        resolution = 800,
                                                        ssaa = 1,
                                                        cuda = TRUE)),
  
  tar_target(bird_anc_recon_joint, inla.posterior.sample(500, ancestral_codes_mod_clads$model,
                                                         selection = as.list(rep(1, 64)) %>% 
                                                           setNames(rownames(ancestral_codes_mod_clads$model$summary.fixed)),
                                                         num.threads = 6)),
  
  tar_target(bird_anc_recon_joint_mat, t(sapply(bird_anc_recon_joint,
                                                function(x) x$latent))),
  
  
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
