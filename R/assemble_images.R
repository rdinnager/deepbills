#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param narrowness_images
#' @param real_col
assemble_images <- function(narrowness_images, real_col = 6) {

  ims <- narrowness_images %>%
    flatten()
  
  crops <- ims %>%
    map(~ bbox(!px.flood(.x, 1, 1)))
  
  mask <- grow(bbox(parany(crops)), 30, 0, 0)
  
  ims_cropped <- map(narrowness_images,
                     ~ map(.x,
                           function(x) crop.bbox(x, mask)))
  
  col_ims <- map(ims_cropped,
                 ~ imappend(.x, "y"))
  
  bords <- px.borders(col_ims[[real_col]], n = 10)
  col_ims[[real_col]] <- imager::colorise(col_ims[[real_col]], bords, "red")
  
  imag <- imappend(col_ims, "x")
  
  imag

}
