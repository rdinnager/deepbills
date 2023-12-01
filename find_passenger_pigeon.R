library(Morpho)
library(Rvcg)
library(readobj)
library(rgl)
library(geometry)
library(proxy)

library(tidyverse)
library(fibre)
library(torch)
library(conflicted)
library(targets)
library(zeallot)

conflict_prefer("dist", "proxy")
conflict_prefer("as.matrix", "base")

bird_beak_avonet <- tar_read(bird_beak_avonet)

codes <- bird_beak_avonet %>%
    select(starts_with("latent_")) %>%
    as.matrix() %>%
    scale()

code_means <- attr(codes, "scaled:center")
code_sds <- attr(codes, "scaled:scale")

clip_mesh_by_mesh <- function(mesh_to_clip, mesh_to_clip_with) {
  #rgl::shade3d(mesh_to_clip, alpha = 0.5, col = "red")
  #rgl::shade3d(mesh_to_clip_with, alpha = 0.5, col = "green")
  dists <- Rvcg::vcgClostKD(t(mesh_to_clip$vb), mesh_to_clip_with)
  clipped <- rgl::clipMesh3d(mesh_to_clip, dists$quality)
  clipped
  #rgl::shade3d(clipped)
}

test_mesh <- read.obj("D:/Projects/deepMorph/data/oriented_scaled_meshes_obj/Abeillia_abeillei_1.obj",
                      convert.rgl = TRUE)[[1]]

beak_mesh <- Rvcg::vcgPlyRead("data/pp_cropped2.ply")
beak_landmarks <- readr::read_csv("data/pp_landmarks2.csv")
shade3d(beak_mesh)
points3d(beak_landmarks[,3:5], col = "red")
main_landmarks <- beak_landmarks[c(3, 4, 1, 2), 3:5] %>%
  as.matrix()
      
mid_beak <- apply(main_landmarks[3:4, ], 2, mean) %>% matrix(nrow = 1)
line_1 <- rbind(mid_beak, main_landmarks[1, ])
      
line_1_length <- dist(line_1)
line_1_vec <- (line_1[1, , drop = FALSE] - line_1[2, , drop = FALSE]) / line_1_length
      
top_line <- main_landmarks[2, ] - mid_beak
      
tt <- top_line %*% t(line_1_vec)
      
point_on_line_1 <- mid_beak + as.vector(tt) * line_1_vec
      
landmarks_3 <- rbind(main_landmarks[1, ], point_on_line_1, main_landmarks[2, ])
landmark_dists <- dist(landmarks_3) %>% as.matrix()
target_3 <- matrix(c(landmark_dists[1, 2], 0, 0,
                     0, 0, 0,
                     0, 0, landmark_dists[2, 3]),
                   nrow = 3, byrow = TRUE)
      
rot <- Morpho::rotmesh.onto(beak_mesh, landmarks_3, target_3, adnormals = TRUE)
shade3d(rot$mesh)
points3d(main_landmarks)
      
# cut by sphere
dist_scale <- dist(target_3[c(1, 3), ])
      
sphere <- Rvcg::vcgSphere(5) %>%
  rgl::translate3d(target_3[1, 1], target_3[1, 2], target_3[1, 3]) %>%
  Morpho::scalemesh(dist_scale)

shade3d(rot$mesh)
shade3d(sphere, col = "blue", alpha = 0.5)
points3d(target_3, col = "red")
      
clipped <- clip_mesh_by_mesh(rot$mesh, sphere) %>%
  Rvcg::vcgClean(sel = 0:6, iterate = TRUE, silent = TRUE)
      
rgl::shade3d(clipped)
rgl::axes3d()
 
pnts <- Morpho::vert2points(clipped)
centroid <- apply(pnts, 2, mean) %>%
  matrix(nrow = 1)
dist_scale <- proxy::dist(centroid, pnts) %>%
  as.vector() %>%
  max()
      
scaled <- clipped %>%
  Morpho::scalemesh(1 / dist_scale, center = "none")

shade3d(scaled)
axes3d()
      
pnts <- Morpho::vert2points(scaled)
      
conv_hull <- pnts[unique(as.vector(convhulln(pnts))),]

mesh <- scaled
generate_sdf_sample <- function(mesh, sphere_rad = 1.5, sphere_centre = c(0.616886145, -0.001077347,  0.145010046),
                                close_sample = 0.005, very_close_sample = 0.0005,
                                n_pnts = 200000, save_folder = "sdf") {
  
  
  n_mesh <- n_pnts*0.4
  mesh_sample <- Rvcg::vcgSample(mesh, SampleNum = n_mesh, type = "mc")
  #points3d(mesh_sample)
  
  n_circ <- n_pnts * 0.2
  circ_sample <- matrix(runif((n_circ*2.5) * 3, -1, 1), ncol = 3)
  norms <- apply(circ_sample, 1, function(x) sqrt(sum(x^2)))
  circ_sample <- circ_sample[norms <= 1, ][1:n_circ, ]
  circ_sample <- t(t(circ_sample) + sphere_centre) * sphere_rad
  #points3d(mesh_sample)
  #points3d(circ_sample, col = "blue")
  
  samp <- rbind(mesh_sample + matrix(rnorm(n_mesh * 3, 0, close_sample), ncol = 3),
                mesh_sample + matrix(rnorm(n_mesh * 3, 0, very_close_sample), ncol = 3),
                circ_sample)
  
  #points3d(samp)
  
  sdf <- Rvcg::vcgClostKD(samp, mesh)
  
  pts <- Morpho::vert2points(mesh)
  max_x <- which.max(pts[ , 1])
  beak_tip <- pts[max_x, , drop = FALSE]
  max_dist <- proxy::dist(beak_tip, pts[-max_x, ]) %>%
    max()
  
  sphere <- Rvcg::vcgSphere(subdivision = 5) %>%
    Morpho::scalemesh(max_dist) %>%
    rgl::translate3d(beak_tip[1, 1], beak_tip[1, 2], beak_tip[1, 3])

  #sdf$quality[points_outside_sphere] <- ifelse(sdf$quality[points_outside_sphere] > 0, sdf_2$quality, sdf$quality[points_outside_sphere])
  
  ## edit sdfs inside mesh:
  ## First find any points inside mesh closer to sphere, reset those to sphere

  samp2 <- samp[sdf$quality >= 0, ]
  sdf_2 <- Rvcg::vcgClostKD(samp2, sphere)
  
  sdf_val <- ifelse(sdf_2$quality < 0, sdf_2$quality, 
                    ifelse(abs(sdf$quality[sdf$quality >= 0]) < abs(sdf_2$quality), sdf$quality[sdf$quality >= 0], sdf_2$quality))
  
  sdf$quality[sdf$quality >= 0] <- sdf_val
  
  ## Second if any point inside both mesh and sphere are closer to bounding box, reset sdf to bounding box
  
  samp3 <- samp[sdf$quality >= 0, ]
  
  bbox_x <- range(pts[ , 1]) 
  bbox_y <- range(pts[ , 2]) 
  bbox_z <- range(pts[ , 3])
  
  dist_to_lower_x <- samp3[ , 1] - bbox_x[1]
  dist_to_upper_x <- -samp3[ , 1] + bbox_x[2]
  dist_to_lower_y <- samp3[ , 2] - bbox_y[1]
  dist_to_upper_y <- -samp3[ , 2] + bbox_y[2]
  dist_to_lower_z <- samp3[ , 3] - bbox_z[1]
  dist_to_upper_z <- -samp3[ , 3] + bbox_z[2]
  
  closest_x <- ifelse(abs(dist_to_lower_x) < abs(dist_to_upper_x), dist_to_lower_x, dist_to_upper_x)
  closest_y <- ifelse(abs(dist_to_lower_y) < abs(dist_to_upper_y), dist_to_lower_y, dist_to_upper_y)
  closest_z <- ifelse(abs(dist_to_lower_z) < abs(dist_to_upper_z), dist_to_lower_z, dist_to_upper_z)
  
  closest_bbox <- ifelse(abs(closest_x) < abs(closest_y), abs(closest_x), abs(closest_y))
  closest_bbox <- ifelse(abs(closest_z) < abs(closest_bbox), abs(closest_z), abs(closest_bbox))
  
  closest_bbox <- ifelse(closest_x >= 0 & closest_y >= 0 & closest_z >= 0, closest_bbox, -closest_bbox)
  
  # rgl::points3d(samp3, color = colourvalues::color_values(as.numeric(closest_bbox >= 0)))
  # rgl::points3d(pts)
  # axes3d()
  # rgl::points3d(samp2, color = colourvalues::color_values(as.numeric(closest_bbox >= 0)))
  
  sdf_val <- ifelse(closest_bbox < 0, closest_bbox, 
                    ifelse(abs(sdf$quality[sdf$quality >= 0]) < abs(closest_bbox), sdf$quality[sdf$quality >= 0], closest_bbox))
  
  #sdf_val <- ifelse(abs(sdf$quality[sdf$quality >= 0]) < abs(closest_bbox), sdf$quality[sdf$quality >= 0], closest_bbox)
  sdf$quality[sdf$quality >= 0] <- sdf_val
  #sdf_val <- ifelse(abs(sdf_val[]) < abs(closest_bbox), sdf_val, closest_bbox) 
  
  #sdf$quality[sdf$quality >= 0] <- sdf_val
  
  # sdf$quality <- ifelse(samp[ , 1] > bbox_x[2] & sdf$quality > 0, samp[ , 1] - bbox_x[2],
  #        sdf$quality)
  
  # sdf$quality <- ifelse((samp[ , 1] < bbox_x[1] | samp[ , 1] > bbox_x[2]) & sdf$quality > 0, -pmin(abs(samp[ , 1] - bbox_x[1]), abs(samp[ , 1] - bbox_x[2])),
  #                       sdf$quality)
  
  # sdf$quality <- ifelse((samp[ , 2] < bbox_y[1] | samp[ , 2] > bbox_y[2]) & sdf$quality > 0, -pmin(abs(samp[ , 2] - bbox_y[1]), abs(samp[ , 2] - bbox_y[2])),
  #        sdf$quality)
  # 
  # sdf$quality <- ifelse((samp[ , 3] < bbox_z[1] | samp[ , 3] > bbox_z[2]) & sdf$quality > 0, -pmin(abs(samp[ , 3] - bbox_z[1]), abs(samp[ , 3] - bbox_z[2])),
  #                       sdf$quality)
  # 
  # rgl::points3d(samp, color = colourvalues::color_values(sdf$quality), size = 5)
  
  ## recenter and scale points
  samp_new <- t(t(samp) - sphere_centre) * (1 / sphere_rad)
  
  sdf_df <- tibble::tibble(x = samp_new[ , 1], y = samp_new[ , 2], z = samp_new[ , 3],
                           sdf = -sdf$quality)
  
  #readr::write_csv(sdf_df, file.path(save_folder, basename(mesh_file)))
  
  #rgl::points3d(sdf_df[,1:3], color = colourvalues::color_values(sdf_df$sdf), size = 5)
  
  # rgl::points3d(samp_new, color = colourvalues::color_values(as.numeric(sdf$quality >= 0)),size = 5)
  # axes3d()
  # 
  # rgl::points3d(samp, color = colourvalues::color_values(as.numeric((samp[ , 1] < bbox_x[1] | samp[ , 1] > bbox_x[2]) & sdf$quality > 0)), size = 5)
  # 
  # rgl::shade3d(sphere, alpha = 0.5, colour = "red")
  return(sdf_df)
}

pp_sdf <- generate_sdf_sample(scaled)

write_csv(pp_sdf, "data/pp_sdf.csv")

pp_ds_ds <- dataset("pp_ds",
                 initialize = function(pp_coord, pp_sdf) {
                   self$coords <- torch_tensor(pp_coord)
                   self$sdf <- torch_tensor(pp_sdf)
                 },
                 .getbatch = function(i) {
                   list(coords = self$coords[i, ],
                        sdf = self$sdf[i, ])
                 },
                 .length = function() {
                   self$coords$size()[[1]]
                 })

pp_sdf2 <- pp_sdf
pp_ds <- pp_ds_ds(pp_sdf %>% mutate(x = pp_sdf2$y, y = pp_sdf2$z, z = -pp_sdf2$x) %>% select(x, y, z) %>% as.matrix(),
                  pp_sdf %>% select(sdf) %>% as.matrix())

pp_dl <- dataloader(pp_ds, 2000, shuffle = TRUE)
pp_dl %>% dataloader_make_iter() %>% dataloader_next()
length(pp_dl)

cvae <- tar_read(two_stage_cvae)$cuda()
sdf_mod <- load_bird_beak_model()
beak_code_finder_mod <- nn_module("beak_finder",
                                  initialize = function(sdf_mod, cvae, latent_dim, sd = 0.01,
                                                        scale_means, scale_sds) {
                                    self$code <- nn_parameter(torch_randn(1, latent_dim) * sd^2)
                                    self$sdf_mod <- sdf_mod$clone()
                                    self$cvae <- cvae$clone()
                                    self$scale_means <- torch_tensor(scale_means, device = "cuda")
                                    self$scale_sds <- torch_tensor(scale_sds, device = "cuda")
                                    self$c_embedding <- torch_zeros(1, self$cvae$c_embed_dim, device = "cuda")
                                    walk(self$sdf_mod$parameters, ~ .x$requires_grad_(FALSE))
                                    walk(self$cvae$parameters, ~ .x$requires_grad_(FALSE))
                                  },
                                  forward = function(coords) {
                                    c(new_code, kl) %<-% self$get_sdf_code()
                                    sdf %<-% self$sdf_mod$forward(coords, new_code)
                                    list(sdf, kl)
                                  },
                                  get_sdf_code = function() {
                                    c(reconstruction, input, mean, log_var) %<-% self$cvae$forward(self$code)
                                    c(loss, reconstruction_loss, kl_loss) %<-% self$cvae$loss_function(reconstruction, input, mean, log_var)
                                    new_code <- (reconstruction * self$scale_sds) + self$scale_means
                                    #kl <- torch_sum(torch_exp(log_var) + torch_square(mean) - log_var, dim = 2L) - 64
                                    list(new_code, loss)
                                  })

latent_dim <- 64L

find_beak_code <- beak_code_finder_mod(sdf_mod$cuda(), cvae, 64L, scale_means = code_means,
                                       scale_sds = code_sds)
find_beak_code <- find_beak_code$cuda()

num_epochs <- 800

lr <- 0.005
optimizer <- optim_adam(find_beak_code$parameters, lr = lr)
scheduler <- lr_one_cycle(optimizer, max_lr = lr,
                          epochs = num_epochs, steps_per_epoch = length(pp_dl),
                          cycle_momentum = FALSE)
# scheduler <- lr_one_cycle(optimizer, max_lr = lr,
#                           epochs = num_epochs, steps_per_epoch = length(pp_dl),
#                           cycle_momentum = FALSE)

sdf_mod2 <- load_bird_beak_model()
# walk(sdf_mod2$parameters, ~ .x$detach()$cpu())
# sdf_mod2 <- sdf_mod2$cpu()
test <- sdf_mod2$get_mesh(find_beak_code$get_sdf_code()[[1]]$cpu())
close3d()
shade3d(test)
#shade3d(rotate3d(scaled), alpha = 0.5, col = "blue")
axes3d()
title3d('main', 'sub', 'xlab', 'ylab', 'zlab')
#find_beak_code <- find_beak_code$cuda()

intermediate_zs <- list()
this_epoch <- 1
for (epoch in 1:num_epochs) {

    batchnum <- 0
    this_epoch <- 0
    coro::loop(for (b in pp_dl) {

        batchnum <- batchnum + 1
        optimizer$zero_grad()

        c(recon, vae_loss) %<-% find_beak_code(b$coords$cuda())
        recon_loss <- torch_mean(abs(recon$clamp(-0.1, 0.1) - b$sdf$cuda()$clamp(-0.1, 0.1))) 
        loss <- 1000*recon_loss + 1e-06 * vae_loss

        if(epoch > this_epoch) {

          intermediate_zs[[epoch]] <- as.matrix(find_beak_code$code$cpu())

        }

        if(batchnum %% 50 == 0) {

            cat("Epoch: ", epoch,
                "  batch: ", batchnum,
                "  loss: ", as.numeric(loss$cpu()),
                "  recon loss: ", as.numeric(recon_loss$cpu()),
                "  vae loss: ", as.numeric(vae_loss$cpu()),
                "\n")
          
          if(epoch %% 25 == 0 && batchnum %% 100 == 0) {
            test <- sdf_mod2$get_mesh(find_beak_code$get_sdf_code()[[1]]$cpu())
            close3d()
            shade3d(test)
            #shade3d(scaled, alpha = 0.5, col = "blue")
            #find_beak_code <- find_beak_code$cuda()
          }

        }
        
        loss$backward()
        optimizer$step()
        #scheduler$step()
        this_epoch <- epoch
    })
}

test <- sdf_mod2$cpu()$get_mesh(find_beak_code$get_sdf_code()[[1]]$cpu(), 200)
close3d()
shade3d(test)

open3d()
shade3d(scaled)

im <- sdf_mod2$cuda()$render_image(find_beak_code$get_sdf_code()[[1]],
                            cuda = TRUE,
                            verbose = TRUE)



sdf_mod <- sdf_mod$cpu()
test <- find_beak_code$sdf_mod$cpu()$get_mesh(find_beak_code$code$cpu())
shade3d(test)
