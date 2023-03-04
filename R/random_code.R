#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param nameme1
#' @param nameme2
random_code <- function(n = 10, min_mag = 0.1, max_mag = 0.2) {

  sample_vec <- rmovMF(n, matrix(0, nrow = 1, ncol = 64))
  mag <- runif(n, min_mag, max_mag)
  sample_vec * mag
  
}
