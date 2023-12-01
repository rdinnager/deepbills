#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bird_beak_avonet
#' @param two_stage_cvae
get_cvae_latents <- function(bird_beak_avonet, two_stage_cvae) {

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
  test_ds <- bb_dataset(codes, trophic)
  test_dl <- dataloader(test_ds, 1012, shuffle = FALSE)
  
  cvae <- two_stage_cvae$cuda()
  
  clatents <- list()
  uclatents <- list()
  
  i <- 0
  coro::loop(for (b in test_dl) {
    i <- i + 1
    with_no_grad({
    c(reconstruction, input, c_mean, c_log_var) %<-% cvae(b$codes$cuda(), b$troph$cuda())  
    clatents[[i]] <- list(means = as.matrix(c_mean$cpu()),
                          logvars = as.matrix(c_log_var$cpu()))
    
    c(reconstruction2, input2, uc_mean, uc_log_var) %<-% cvae(b$codes$cuda())
    uclatents[[i]] <- list(means = as.matrix(uc_mean$cpu()),
                           logvars = as.matrix(uc_log_var$cpu()))
    })
    print(i)
    
  })

  mean_df <- do.call(rbind, map(clatents, "means")) %>%
    as.data.frame() %>%
    mutate(type = "conditional",
           trophic = trophic) %>%
    bind_rows(do.call(rbind, map(uclatents, "means")) %>%
                as.data.frame() %>%
                mutate(type = "unconditional",
                       trophic = trophic))
  
  logvar_df <- do.call(rbind, map(clatents, "logvars")) %>%
    as.data.frame() %>%
    mutate(type = "conditional",
           trophic = trophic) %>%
    bind_rows(do.call(rbind, map(uclatents, "logvars")) %>%
                as.data.frame() %>%
                mutate(type = "unconditional",
                       trophic = trophic))
  
  uc_var_mean <- logvar_df %>%
    filter(type == "unconditional") %>%
    select(-type, -trophic) %>%
    map_dbl(mean) %>%
    tibble(mean_var = exp(.))
  
  manifold <- which(uc_var_mean$mean_var < 0.5)
  
  return(list(mean = mean_df, logvar = logvar_df, manifold = manifold))
  
  uc_var_mean <- logvar_df %>%
    filter(type == "unconditional") %>%
    select(-type, -trophic) %>%
    map_dbl(mean) %>%
    tibble(mean_var = exp(.))
  
  manifold <- which(uc_var_mean$mean_var < 0.5)
  
  ggplot(uc_var_mean, aes(mean_var)) +
    geom_histogram(bins = 50)
  
  c_var_mean <- logvar_df %>%
    filter(type == "conditional") %>%
    select(-type) %>%
    pivot_longer(-trophic, values_to = "val", names_to = "var") %>%
    mutate(val = exp(val)) %>%
    group_by(trophic, var) %>%
    summarise(means = mean(val))
  
  c_var_summ <- c_var_mean %>%
    group_by(trophic) %>%
    summarise(n_manifold = sum(means < 0.5))
  
  ggplot(mean_df %>% filter(type == "unconditional"), aes(V1, V2)) +
    geom_point(aes(colour = trophic)) +
    theme_minimal()
    
  
  ggplot(uc_var_mean, aes(mean_var)) +
    geom_histogram(bins = 50)
  
  qqnorm(mean_df[ , manifold][[1]])
  qqline(mean_df[ , manifold][[1]])
  
  qqnorm(mean_df[ , manifold][[2]])
  qqline(mean_df[ , manifold][[2]])
  
  qqnorm(mean_df[ , manifold][[3]])
  qqline(mean_df[ , manifold][[3]])
  
  qqnorm(mean_df[ , manifold][[4]])
  qqline(mean_df[ , manifold][[4]])
  
  qqnorm(mean_df[ , manifold][[5]])
  qqline(mean_df[ , manifold][[5]])
  
  qqnorm(mean_df[ , manifold][[6]])
  qqline(mean_df[ , manifold][[6]])
  
  qqnorm(mean_df[ , manifold][[7]])
  qqline(mean_df[ , manifold][[7]])
  
  qqnorm(mean_df[ , manifold][[8]])
  qqline(mean_df[ , manifold][[8]])
  
  qqnorm(mean_df[ , manifold][[9]])
  qqline(mean_df[ , manifold][[9]])
  
  qqnorm(mean_df[ , manifold][[10]])
  qqline(mean_df[ , manifold][[10]])
  
  qqnorm(mean_df[ , manifold][[11]])
  qqline(mean_df[ , manifold][[11]])

  qqnorm(mean_df[ , manifold][[12]]) ##
  qqline(mean_df[ , manifold][[12]])
  
  qqnorm(mean_df[ , manifold][[13]]) ##
  qqline(mean_df[ , manifold][[13]])
  
  qqnorm(mean_df[ , manifold][[14]])
  qqline(mean_df[ , manifold][[14]])
  
  qqnorm(mean_df[ , manifold][[15]])
  qqline(mean_df[ , manifold][[15]])
  
  qqnorm(mean_df[ , manifold][[16]])
  qqline(mean_df[ , manifold][[16]])
  
  qqnorm(mean_df[ , manifold][[17]])
  qqline(mean_df[ , manifold][[17]])
  
  qqnorm(mean_df[ , manifold][[18]])
  qqline(mean_df[ , manifold][[18]])
  
  qqnorm(mean_df[ , manifold][[19]])
  qqline(mean_df[ , manifold][[19]])
  
  qqnorm(mean_df[ , manifold][[20]])
  qqline(mean_df[ , manifold][[20]])
  
  qqnorm(mean_df[ , manifold][[21]])
  qqline(mean_df[ , manifold][[21]])
  
  qqnorm(mean_df[ , manifold][[22]])
  qqline(mean_df[ , manifold][[22]])
  
  qqnorm(mean_df[ , manifold][[23]])
  qqline(mean_df[ , manifold][[23]])
  
  qqnorm(mean_df[ , manifold][[24]]) ##
  qqline(mean_df[ , manifold][[24]])
  
  qqnorm(mean_df[ , manifold][[25]])
  qqline(mean_df[ , manifold][[25]])
  
  qqnorm(mean_df[ , manifold][[26]])
  qqline(mean_df[ , manifold][[26]])
  
  qqnorm(mean_df[ , manifold][[27]])
  qqline(mean_df[ , manifold][[27]])
  
  qqnorm(mean_df[ , manifold][[28]])
  qqline(mean_df[ , manifold][[28]])
  
  qqnorm(mean_df[ , manifold][[29]])
  qqline(mean_df[ , manifold][[29]])
}
