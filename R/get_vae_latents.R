#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param two_stage_vae_mod
#' @param bird_beak_avonet
#' @return
#' @author rdinnager
#' @export
get_vae_latents <- function(two_stage_vae_mod, bird_beak_avonet) {

  codes <- bird_beak_avonet %>%
    select(starts_with("latent_")) %>%
    as.matrix() %>%
    scale()
  
  trophic <- bird_beak_avonet %>%
    select(Trophic.Niche) %>%
    mutate(Trophic.Niche = factor(Trophic.Niche)) %>%
    pull(Trophic.Niche)
  
  bb_dataset <- dataset(name = "bb_ds",
                           initialize = function(codes, trophic) {
                             self$codes <- torch_tensor(codes)
                             self$troph <- torch_tensor(trophic) 
                           },
                           .getbatch = function(i) {
                             list(codes = self$codes[i, ], troph = self$troph[i])
                           },
                           .length = function() {
                             self$codes$size()[[1]]
                           })
  train_ds <- bb_dataset(codes, trophic)
  train_dl <- dataloader(train_ds, 1012, shuffle = TRUE)
  
  vae <- two_stage_vae_mod$cuda()
  
  post_dl <- dataloader(train_ds, 506, shuffle = FALSE)
  post_it <- as_iterator(post_dl)
  
  all_mean_var <- list()
  it <- 0
  loop(for(i in post_dl) {
    it <- it + 1
    with_no_grad({
      all_mean_var[[it]] <- vae(i$codes$cuda())
    })
    print(it)
  })
  
  all_means <- map(all_mean_var, 3) %>%
    map(~ as.matrix(.x$cpu())) %>%
    do.call(rbind, .)
  
  all_vars <- map(all_mean_var, 4) %>%
    map(~ as.matrix(torch_exp(.x)$cpu())) %>%
    do.call(rbind, .)
  
  vars <- apply(all_vars, 2, mean)
  
  vars_df <- tibble(`Mean Variance` = vars) %>%
    mutate(Status = ifelse(`Mean Variance` < 0.5, "Active", "Ignored"))
  
  p <- ggplot(vars_df, aes(`Mean Variance`)) +
    geom_histogram(aes(fill = Status)) +
    ylab("Count") +
    theme_minimal() +
    theme(legend.position = c(0.5, 0.5))
  
  plot(p)
  
  #active_dims <- which(vars < 0.5)
  #plot(all_means[ , active_dims][ , 1:2])
  #rgl::points3d(all_means)
  
 # write_rds(active_dims, "data/active_dims_16dim.rds")
  
  ######### Save latent codes ################
  
  active_dims <- which(vars < 0.5)
  latent_means <- all_means[ , active_dims]
  latent_vars <- all_vars[ , active_dims]
  
  colnames(latent_means) <- colnames(latent_vars) <- paste0("latent_", 1:ncol(latent_means))
  
  latent_df <- as.data.frame(latent_means) %>%
    mutate(Species = bird_beak_avonet$label) %>%
    pivot_longer(-Species, names_to = "latent_dim", values_to = "latent_mean") %>%
    left_join(as.data.frame(latent_vars) %>%
                mutate(Species = bird_beak_avonet$label) %>%
                pivot_longer(-Species, names_to = "latent_dim", values_to = "latent_var")) %>%
    mutate(latent_lower = latent_mean - 1.96 * sqrt(latent_var),
           latent_upper = latent_mean + 1.96 * sqrt(latent_var))
  
  #write_csv(latent_df, "data/bills_vae_latent_codes_16dim.csv")
  latent_df_means <- as.data.frame(latent_means) %>%
    mutate(Species = bird_beak_avonet$label) %>%
    select(Species, everything())
  
  #write_csv(latent_df_means, "data/bills_vae_latent_codes_means_only_16dim.csv")

  list(latent_df = latent_df, latent_df_means = latent_df_means,
       active_dims = active_dims, p = p)

}
