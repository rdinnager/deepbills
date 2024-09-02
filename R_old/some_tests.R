get_points_in_unit_sphere <- function(n) {

  x <- matrix(runif(as.integer(n * 2.5)*3) * 2 - 1, ncol = 3)
  mask <- which(sqrt(rowSums(x^2)) < 1)
  mask <- mask[1:n]
  x <- x[mask, ]

  return(x)
}

library(phyf)
latent <- bird_beak_codes %>%
  dplyr::select(dplyr::starts_with("latent_")) %>%
  dplyr::slice(500) %>%
  unlist()
 

points <- as.matrix(get_points_in_unit_sphere(100000))

sdf <- SDFNet(points, latent)

test <- get_meshes_from_latent(latent)

params_as_Robs <- purrr::map(SDFNet$state_dict(),
                             ~ .x$cpu()$numpy())

readr::write_rds(params_as_Robs, "data/sdf_net_params.rds")
