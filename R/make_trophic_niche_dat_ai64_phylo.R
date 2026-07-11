#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bird_beak_avonet
#' @param bad_birds
#' @return
#' @author rdinnager
#' @export
make_trophic_niche_dat_ai64_phylo <- function(bird_beak_avonet, bad_birds) {
  
  tree_graph <- Phylo2DirectedGraph(pf_as_phylo(bird_beak_avonet))
  PEM <- PEM.build(tree_graph)
  
  PEM_df <- as.data.frame(PEM)
  PEM_df <- PEM_df / sd(unlist(PEM_df))
  
  PEM_df <- PEM_df |>
    mutate(Species = rownames(PEM_df))

  dat <- bird_beak_avonet %>%
    select(label,
           Trophic.Niche,
           starts_with("latent_")) %>%
    mutate(codes = scale(across(starts_with("latent_"), ~ .x))) %>%
    filter(!label %in% bad_birds) |>
    select(Species = label, Trophic.Niche, codes) |>
    mutate(codes = as.data.frame(codes)) |>
    unnest(codes) |>
    left_join(PEM_df, by = c("Species" = "Species")) |>
    group_by(Trophic.Niche) |>
    mutate(weights = 1 / n()) |>
    ungroup() |>
    mutate(weights = importance_weights(weights))
  
  dat

}
