#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param two_stage_cvae
#' @param bird_beak_avonet
get_latent_samples <- function(trophic_lev, two_stage_cvae, pal_col, bird_beak_avonet) {

  codes <- bird_beak_avonet %>%
      select(starts_with("latent_")) %>%
      as.matrix() %>%
      scale()
  
  code_means <- torch_tensor(attr(codes, "scaled:center"))
  code_sds <- torch_tensor(attr(codes, "scaled:scale")  )
  
  # trophic <- bird_beak_avonet %>%
  #   select(Trophic.Niche) %>%
  #   mutate(Trophic.Niche = factor(Trophic.Niche)) %>%
  #   pull(Trophic.Niche)
  # 
  # trophic_levs <- torch_tensor(unique(trophic))
  
  num_samps <- 10000
  
  cvae <- two_stage_cvae$cuda()
  lev <- torch_tensor(trophic_lev)$`repeat`(c(num_samps))
  gen_codes <- cvae$decode(torch_randn(num_samps, 64)$cuda(), lev$cuda())
  
  gen_codes <- (gen_codes * code_sds$cuda()) + code_means$cuda()
  
  sdf_mod <- load_bird_beak_model()$cuda()

  
    # test_mesh <- sdf_mod$get_mesh(gen_codes[10, , drop = FALSE], resolution = 200,
  #                               cuda = TRUE, batch_size = 200000,
  #                               smooth = FALSE) %>%
  #   Rvcg::vcgSmooth(type = "HClaplace")
  # rgl::shade3d(test_mesh %>%
  #                rgl::rotate3d(pi, 1, 1, 0), col = "yellow")
  # 
  # rgl::shade3d(test_mesh %>%
  #                rgl::rotate3d(pi/2, 0, 0, 1) %>%
  #                rgl::rotate3d(pi/2, 0, 1, 0) %>%
  #                rgl::rotate3d(3*pi/4, 0, 0, 1), col = "yellow")
  
  do_transform <- function(transf, mesh) {
    switch(transf,
           mesh %>%
                 rgl::rotate3d(pi, 1, 0, 0),
           mesh %>%
                 rgl::rotate3d(pi/2, 0, 0, 1) %>%
                 rgl::rotate3d(pi/2, 0, 1, 0),
           mesh %>%
                 rgl::rotate3d(-pi/2, 0, 0, 1) %>%
                 rgl::rotate3d(-pi/2, 0, 1, 0),
           mesh %>%
             rgl::rotate3d(pi/2, 0, 0, 1) %>%
             rgl::rotate3d(pi/2, 0, 1, 1),
           mesh %>%
                 rgl::rotate3d(pi/2, 0, 0, 1) %>%
                 rgl::rotate3d(pi/2, 0, 1, 0) %>%
                 rgl::rotate3d(3*pi/4, 0, 0, 1),
           mesh %>%
                 rgl::rotate3d(pi/2, 0, 0, 1) %>%
                 rgl::rotate3d(pi/2, 0, 1, 0) %>%
                 rgl::rotate3d(pi/4, 0, 0, 1)
           )
  }
  
  generate_beak_im <- function(i) {
    test_mesh <- sdf_mod$get_mesh(gen_codes[sample.int(nrow(gen_codes), 1), , drop = FALSE], resolution = 200,
                                cuda = TRUE, batch_size = 200000,
                                smooth = FALSE) %>%
      Rvcg::vcgSmooth(type = "HClaplace")
    shade3d(do_transform(sample.int(6, 1), test_mesh), col =  pal_col)
    png_file <- tempfile(fileext = ".png")
    rgl::snapshot3d(filename = png_file, width = 256, height = 256,
                  webshot = FALSE)
    rgl::close3d()
  
    im2 <- imager::load.image(png_file)
    im <- imager::imfill(256, 256, val = c(0, 0, 0, 1))
    im[ , , , 1:3] <- im2 
    im[imager::R(im) == 1 & imager::G(im) == 1 & imager::B(im) == 1] <- 0
  
    im  
  }
  
  test <- impac(generate_beak_im, progress = TRUE, show_every = 10, bg = "white",
                scales = c(0.75, 0.75*0.95, 0.75*(0.95^2), 0.75*(0.95^3)), scale_fun = function(s, i, c) {     
                  if (c < (i * 0.5)) {         
                    c(s, rep(min(s) * 0.95, floor(1 / min(s))))     
                  } else 
                  { s } },
                width = 1024, height = 1024)
  
  # test <- generate_beak_im(2)
  # plot(test)
  # sdf_mod <- load_bird_beak_model()$cuda()
  # test_image <- sdf_mod$render_image(gen_codes[4, , drop = FALSE], verbose = TRUE, cuda = TRUE, resolution = 256, ssaa = 1, threshold = 0.001, max_ray_move = 0.005)
  # plot(test_image)
  
  test
  

}
