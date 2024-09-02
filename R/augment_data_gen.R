#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_dat_pca_all
#' @param two_stage_cvae_pc
#' @return
#' @author rdinnager
#' @export
augment_data_gen_pc <- function(trophic_niche_dat_pca_all, two_stage_cvae_pc,
                             trophic_levs2) {

  #code_means <- torch_tensor(trophic_niche_dat_pca_scaling[["scaled:center"]])
  #code_sds <- torch_tensor(trophic_niche_dat_pca_scaling[["scaled:scale"]])
  
  summ <- trophic_niche_dat_pca_all |>
    group_by(Trophic.Niche) |>
    summarise(count = n())
  # trophic <- bird_beak_avonet %>%
  #   select(Trophic.Niche) %>%
  #   mutate(Trophic.Niche = factor(Trophic.Niche)) %>%
  #   pull(Trophic.Niche)
  # 
  # trophic_levs <- torch_tensor(unique(trophic))
  
  num_samps <- max(summ$count)
  
  cvae <- two_stage_cvae_pc$cuda()
  lev <- torch_tensor(trophic_levs2)$`repeat`(c(num_samps))
  gen_codes <- cvae$decode(torch_randn(lev$size()[[1]], 64)$cuda(), lev$cuda())
  
  #gen_codes <- (gen_codes * code_sds$cuda()) + code_means$cuda()
  
  gen_codes <- as.matrix(gen_codes$cpu())
  colnames(gen_codes) <- paste0("PC", 1:8)
  code_df <- as.data.frame(gen_codes) |>
    mutate(Trophic.Niche = trophic_levs2[as.numeric(lev$cpu())]) |>
    left_join(summ) |>
    mutate(keep_num = num_samps - count, type = "augmented") |>
    group_by(Trophic.Niche) |>
    dplyr::filter(1:n() <= keep_num) |>
    ungroup()
  
  new_dat <- trophic_niche_dat_pca_all |>
    mutate(type = "original") |>
    bind_rows(code_df) |>
    mutate(strata = paste0(Trophic.Niche, " - ", type))
  
  #plot(gen_codes[as.logical(lev==5),1:2])
  #points(trophic_niche_dat_pca_all |> dplyr::filter(Trophic.Niche == trophic_levs2[5]) |> dplyr::select(PC1, PC2), col = "red")
  
  return(new_dat)

}


augment_data_gen_ai64 <- function(trophic_niche_dat_ai64_all, two_stage_cvae_ai64,
                             trophic_levs2) {

  #code_means <- torch_tensor(trophic_niche_dat_pca_scaling[["scaled:center"]])
  #code_sds <- torch_tensor(trophic_niche_dat_pca_scaling[["scaled:scale"]])
  
  summ <- trophic_niche_dat_ai64_all |>
    group_by(Trophic.Niche) |>
    summarise(count = n())
  # trophic <- bird_beak_avonet %>%
  #   select(Trophic.Niche) %>%
  #   mutate(Trophic.Niche = factor(Trophic.Niche)) %>%
  #   pull(Trophic.Niche)
  # 
  # trophic_levs <- torch_tensor(unique(trophic))
  
  num_samps <- max(summ$count)
  
  cvae <- two_stage_cvae_ai64$cuda()
  lev <- torch_tensor(trophic_levs2)$`repeat`(c(num_samps))
  gen_codes <- cvae$decode(torch_randn(lev$size()[[1]], 64)$cuda(), lev$cuda())
  
  #gen_codes <- (gen_codes * code_sds$cuda()) + code_means$cuda()
  
  gen_codes <- as.matrix(gen_codes$cpu())
  colnames(gen_codes) <- paste0("latent_code_", 1:64)
  code_df <- as.data.frame(gen_codes) |>
    mutate(Trophic.Niche = trophic_levs2[as.numeric(lev$cpu())]) |>
    left_join(summ) |>
    mutate(keep_num = num_samps - count, type = "augmented") |>
    group_by(Trophic.Niche) |>
    dplyr::filter(1:n() <= keep_num) |>
    ungroup()
  
  new_dat <- trophic_niche_dat_ai64_all |>
    mutate(type = "original") |>
    bind_rows(code_df) |>
    mutate(strata = paste0(Trophic.Niche, " - ", type))
  
  #plot(gen_codes[as.logical(lev==5),1:2])
  #points(trophic_niche_dat_pca_all |> dplyr::filter(Trophic.Niche == trophic_levs2[5]) |> dplyr::select(PC1, PC2), col = "red")
  
  return(new_dat)

}