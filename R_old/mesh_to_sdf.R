library(Morpho)
library(Rvcg)
library(readobj)
library(rgl)
library(geometry)
library(proxy)

## Step one is to iterate over each mesh, reorient the beak by the landmarks, and then
## rescale each beak to constant beak length. I will also save the convex hull of 
## each resulting mesh to use in the next step.

beak_files <- list.files("data/objs/", full.names = TRUE) 

beak_file <- beak_files[1]

clip_mesh_by_mesh <- function(mesh_to_clip, mesh_to_clip_with) {
  #rgl::shade3d(mesh_to_clip, alpha = 0.5, col = "red")
  #rgl::shade3d(mesh_to_clip_with, alpha = 0.5, col = "green")
  dists <- Rvcg::vcgClostKD(t(mesh_to_clip$vb), mesh_to_clip_with)
  clipped <- rgl::clipMesh3d(mesh_to_clip, dists$quality)
  clipped
  #rgl::shade3d(clipped)
}

#beak_file <- beak_files[2]
#plot_folder <- "beak_plots_trimmed"
reorient_beak <- function(beak_file, mesh_folder = "data/oriented_scaled_meshes_obj",
                          landmark_folder = "data/landmarks_oriented", hull_folder = "data/conv_hulls") {
  
  if(!dir.exists(mesh_folder)) {
    dir.create(mesh_folder)
  }
  if(!dir.exists(landmark_folder)) {
    dir.create(landmark_folder)
  }
  
  bird_file <- basename(beak_file)

  new_beak_file <- file.path(mesh_folder, bird_file)
  hull_file <- file.path(hull_folder, gsub(".obj", ".csv", bird_file))
  
  if(!dir.exists(hull_folder)) {
    dir.create(hull_folder)
  }
  
  
  
  if(!file.exists(new_beak_file)) {
    
    beak_mesh <- try(readobj::read.obj(beak_file, convert.rgl = TRUE)[[1]])
    suppressMessages(beak_landmarks <- try(readr::read_delim(file.path("data/landmarks/Markup_1", gsub(".obj", ".txt", bird_file, fixed = TRUE)), 
                                                             " ",
                                                             col_names = c("X", "Y", "Z"))))
    
    if(inherits(beak_mesh, 'try-error') | inherits(beak_landmarks, 'try-error')) {
      readr::write_lines(bird_file, "data/bad_models.txt", append = TRUE)
    } else {
      
      main_landmarks <- beak_landmarks[1:4, ] %>%
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
      
      L_landmark_folder <- file.path(landmark_folder, "L_landmarks")
      if(!dir.exists(L_landmark_folder)) {
        dir.create(L_landmark_folder)
      }
      
      all_landmark_folder <- file.path(landmark_folder, "all_landmarks")
      if(!dir.exists(all_landmark_folder)) {
        dir.create(all_landmark_folder)
      }
      
      readr::write_csv(as.data.frame(target_3), file.path(L_landmark_folder, gsub(".obj", ".csv", bird_file, fixed = TRUE)),
                       col_names = FALSE)
      
      rot_landmarks <- beak_landmarks %>% as.matrix() %>%
        Morpho::rotonmat(landmarks_3, target_3)
      
      readr::write_csv(as.data.frame(rot_landmarks), file.path(all_landmark_folder, gsub(".obj", ".csv", bird_file, fixed = TRUE)),
                       col_names = FALSE)
      
      rot <- Morpho::rotmesh.onto(beak_mesh, landmarks_3, target_3, adnormals = TRUE)
      
      # rgl::shade3d(rot$mesh, alpha = 0.1)
      # rgl::lines3d(target_3, col = "red")
      # rgl::points3d(rot_landmarks, size = 4, col = "blue")
      # rgl::shade3d(sphere, col = "red", alpha = 0.5)
      
      ## cut by sphere
      dist_scale <- dist(target_3[c(1, 3), ])
      
      sphere <- Rvcg::vcgSphere(5) %>%
        rgl::translate3d(target_3[1, 1], target_3[1, 2], target_3[1, 3]) %>%
        Morpho::scalemesh(dist_scale)
      
      clipped <- clip_mesh_by_mesh(rot$mesh, sphere) %>%
        Rvcg::vcgClean(sel = 0:6, iterate = TRUE, silent = TRUE)
      
      # rgl::shade3d(clipped)
      # rgl::axes3d()
      # 
      # rgl::shade3d(scaled)
      # rgl::axes3d()
      # rgl::points3d(conv_hull, col = "green")
      
      pnts <- Morpho::vert2points(clipped)
      
      centroid <- apply(pnts, 2, mean) %>%
        matrix(nrow = 1)
      dist_scale <- proxy::dist(centroid, pnts) %>%
        as.vector() %>%
        max()
      
      scaled <- clipped %>%
        Morpho::scalemesh(1 / dist_scale, center = "none")
      
      pnts <- Morpho::vert2points(scaled)
      
      conv_hull <- pnts[unique(as.vector(convhulln(pnts))),]
        
      Rvcg::vcgObjWrite(scaled, new_beak_file)
      readr::write_csv(conv_hull %>% as.data.frame(), hull_file)
        
    }
  }
  return(invisible(NULL))
}

library(furrr)

future::plan(future::multisession)

furrr::future_map(beak_files, ~reorient_beak(.x),
                  .progress = TRUE)


#### Find sphere that fits all meshes in it
library(pbapply)
hulls <- list.files("data/conv_hulls", full.names = TRUE)
all_pnts <- pbapply::pblapply(hulls, readr::read_csv) %>%
  dplyr::bind_rows() %>%
  as.matrix()

centroid <- apply(all_pnts, 2, mean)
all_pnts_cent <- t(t(all_pnts) - centroid)

norms <- apply(all_pnts_cent, 1, function(x) sqrt(sum(x^2)))
max_norm <- max(norms)

## sphere containing all meshes has radius 1.5 and centred at 0.616886145, -0.001077347,  0.145010046

##### Generate SDF samples for all meshes
library(Morpho)
library(Rvcg)
library(readobj)
library(rgl)
library(geometry)
library(proxy)

meshes <- list.files("data/oriented_scaled_meshes_obj/", full.names = TRUE)

mesh_file <- meshes[100]
generate_sdf_sample <- function(mesh_file, sphere_rad = 1.5, sphere_centre = c(0.616886145, -0.001077347,  0.145010046),
                                close_sample = 0.005, very_close_sample = 0.0005,
                                n_pnts = 200000, save_folder = "sdf") {
  
  
  mesh <- try(readobj::read.obj(mesh_file, convert.rgl = TRUE)[[1]])
  
  n_mesh <- n_pnts*0.4
  mesh_sample <- Rvcg::vcgSample(mesh, SampleNum = n_mesh, type = "mc")
  
  n_circ <- n_pnts * 0.2
  circ_sample <- matrix(runif((n_circ*2.5) * 3, -1, 1), ncol = 3)
  norms <- apply(circ_sample, 1, function(x) sqrt(sum(x^2)))
  circ_sample <- circ_sample[norms <= 1, ][1:n_circ, ]
  circ_sample <- t(t(circ_sample) + sphere_centre) * sphere_rad
  
  samp <- rbind(mesh_sample + matrix(rnorm(n_mesh * 3, 0, close_sample), ncol = 3),
                mesh_sample + matrix(rnorm(n_mesh * 3, 0, very_close_sample), ncol = 3),
                circ_sample)
  
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
  # rgl::points3d(samp, color = colourvalues::color_values(sdf$quality),
  #               size = 5)
  
  ## recenter and scale points
  samp_new <- t(t(samp) - sphere_centre) * (1 / sphere_rad)
  
  sdf_df <- tibble::tibble(x = samp_new[ , 1], y = samp_new[ , 2], z = samp_new[ , 3],
                           sdf = -sdf$quality)
  
  readr::write_csv(sdf_df, file.path(save_folder, basename(mesh_file)))
  
  # rgl::points3d(samp_new, color = colourvalues::color_values(as.numeric(sdf$quality >= 0)),
  #               size = 5)
  # axes3d()
  # 
  # rgl::points3d(samp, color = colourvalues::color_values(as.numeric((samp[ , 1] < bbox_x[1] | samp[ , 1] > bbox_x[2]) & sdf$quality > 0)),
  #               size = 5)
  # 
  # rgl::shade3d(sphere, alpha = 0.5, colour = "red")
  return(file.path(save_folder, basename(mesh_file)))
}

library(furrr)
future::plan(future::multiprocess)

sdf_files <- furrr::future_map(meshes, ~try(generate_sdf_sample(.x, save_folder = "data/sdf")))


## test some random ones out
sdf_1 <- readr::read_csv("data/sdf/Tricholestes_criniger_1.obj")
rgl::points3d(sdf_1[ , 1:3] %>% as.matrix(), color = colourvalues::color_values(sdf_1$sdf))
rgl::axes3d()

rgl::points3d(sdf_1[ , 1:3] %>% as.matrix(), color = colourvalues::color_values(as.numeric(sdf_1$sdf >= 0)))
rgl::axes3d()

########## combine sdf samples into one torch file #################
library(reticulate)
library(readr)
library(pbapply)
library(dplyr)
library(furrr)

csvs <- list.files("data/sdf", full.names = TRUE)

future::plan(future::multiprocess)

all_sdf <- furrr::future_map_dfr(csvs, ~readr::read_csv(.x, col_types = cols(
  x = col_double(),
  y = col_double(),
  z = col_double(),
  sdf = col_double()
)), .progress = TRUE)

torch <- reticulate::import("torch")
scipy <- reticulate::import("scipy.spatial.transform")

rot_1 <- scipy$Rotation$from_euler('y', 90, degrees = TRUE)$as_dcm()
rot_2 <- scipy$Rotation$from_euler('x', -90, degrees = TRUE)$as_dcm()
rot <- rot_1 %*% rot_2
#rot <- scipy$Rotation$from_euler('x', -90, degrees = TRUE)$as_dcm()

coord_mat <- t(rot %*% t(all_sdf[, 1:3] %>% as.matrix))
coords <- torch$tensor(reticulate::np_array(coord_mat, dtype = "float32"))
torch$save(coords, "data/torch/sdf_points.to")

values <- torch$tensor(reticulate::np_array(all_sdf$sdf, dtype = "float32"))
torch$save(values, "data/torch/sdf_values.to")

test <- torch$load("/data/dinnage/Projects/shapegan/data/sdf_points.to")
test2 <- torch$load("/data/dinnage/Projects/shapegan/data/sdf_values.to")