#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param codes_anc
#' @param return.tau.mat
#' @param return.u
gaussianize <- function(codes_anc, return.tau.mat = TRUE, return.u = TRUE) {

  tt <- codes_anc
  ttg <- Gaussianize(tt[!is.na(tt[, 1]), ], return.tau.mat = return.tau.mat,
                     return.u = return.u)
  tt[!is.na(tt[, 1]), ] <- ttg$input

  list(input = tt, tau.mat = ttg$tau.mat)
  
}
