#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_latent_samples
make_trophic_sample_plot <- function(trophic_latent_samples_ims,
                                     trophic_levs2) {

  ims <- trophic_latent_samples_ims[-9]
  
  ims <- map(ims, ~ crop.borders(.x, nPix = 2))
  
  ims <- list(ims[1:3],
              ims[4:6],
              ims[7:9])
  
  im_rows <- map(ims, ~imappend(.x, "x"))
  
  im_final <- imappend(im_rows, "y")
  
  im_final

}
