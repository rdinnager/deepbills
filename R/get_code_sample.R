#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param codes
#' @param nameme1
get_code_sample <- function(codes, n = 50000) {

  norms <- sqrt(rowSums(codes^2)) 
  sample_vec <- rmovMF(n, matrix(0, nrow = 1, ncol = 64))
  res <- sample_vec * (sample(norms, n, replace = TRUE) + rnorm(n, sd = 0.025))
  res[res < 0] <- 0
  res
  
}
