#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param biggest_traj_clads
#' @param beak_codes_aces
make_trajectory_codes <- function(biggest_traj_clads, beak_codes_aces, codes_anc_st_clads) {

  m <- attr(codes_anc_st_clads, "scaled:center")
  s <- attr(codes_anc_st_clads, "scaled:scale")
  beak_codes_aces_coll <- beak_codes_aces %>%
    nest(codes = starts_with("codes.latent_code_")) %>%
    mutate(codes = map(codes, ~ (.x * s) + m))
  edges <- pf_ends(biggest_traj_clads$phlo)
  edge_codes <- edges %>%
    left_join(beak_codes_aces_coll %>%
                select(start = label, start_codes = codes)) %>%
    left_join(beak_codes_aces_coll %>%
                select(end = label, end_codes = codes)) %>%
    mutate(start_codes = map(start_codes, as.matrix),
           end_codes = map(end_codes, as.matrix))
  
  interp_codes <- map2(edge_codes$start_codes, edge_codes$end_codes,
                       ~ interpolate_codes(.x, .y))
  interp_codes

}

interpolate_codes <- function(start, end, n = 10) {
  steps <- seq(0, 1, length.out = n)
  new_codes <- map(steps,
                   ~ ((1 - .x) * start) + (.x * end))
  new_codes <- do.call(rbind, new_codes)
  new_codes
}
