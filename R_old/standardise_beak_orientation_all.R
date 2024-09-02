library(Morpho)
library(readobj)
library(rgl)

beak_files <- list.files("data/objs/", full.names = TRUE) 

#beak_file <- beak_files[33]
#plot_folder <- "beak_plots_cut"
reorient_beak <- function(beak_file, cut_extraneous = TRUE, mesh_folder = "data/ply_oriented_cut",
                          landmark_folder = "data/landmarks_oriented", plot_folder = NULL) {
  
  if(!dir.exists(mesh_folder)) {
    dir.create(mesh_folder)
  }
  if(!dir.exists(landmark_folder)) {
    dir.create(landmark_folder)
  }
  
  bird_file <- basename(beak_file)
  if(!cut_extraneous) {
    new_beak_file <- file.path(mesh_folder, bird_file)
  } else {
    cut_1_folder <- file.path(mesh_folder, "cut_style_1")
    if(!dir.exists(cut_1_folder)) {
      dir.create(cut_1_folder)
    }
    
    cut_2_folder <- file.path(mesh_folder, "cut_style_2")
    if(!dir.exists(cut_2_folder)) {
      dir.create(cut_2_folder)
    }
    new_beak_file <- file.path(cut_1_folder, bird_file)
  }
  
  if(!is.null(plot_folder)) {
    if(!dir.exists(plot_folder)) {
      dir.create(plot_folder)
    }
  }
  
  
  if(!file.exists(gsub(".obj", ".ply", new_beak_file, fixed = TRUE))) {
    
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
      
      line_1_vec <- (line_1[1, , drop = FALSE] - line_1[2, , drop = FALSE]) / dist(line_1)
      
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
      
      if(cut_extraneous) {
      
        cut_plane_1 <- rot_landmarks[c(2, 4, 3), ] %>% as.matrix()
        cut_1 <- Morpho::cutMeshPlane(rot$mesh, v1 = cut_plane_1[1, ], v2 = cut_plane_1[2, ], v3 = cut_plane_1[3, ])
        
        hindmost <- rot_landmarks[which.min(rot_landmarks[ , 1]), , drop = FALSE]
        cut_plane_2 <- rbind(hindmost - c(0, 1, 0),
                                   hindmost - c(0, 0, 1),
                                   hindmost - c(0, -1, 0))
        
        cut_1 <- Morpho::cutMeshPlane(cut_1, v1 = cut_plane_2[1, ], v2 = cut_plane_2[2, ], v3 = cut_plane_2[3, ])
        
        cut_2 <- Morpho::cutMeshPlane(rot$mesh, v1 = cut_plane_2[1, ], v2 = cut_plane_2[2, ], v3 = cut_plane_2[3, ])
        
        cut_1 <- Rvcg::vcgIsolated(cut_1, silent = TRUE, facenum = 0.1*ncol(cut_1$it))
        cut_2 <- Rvcg::vcgIsolated(cut_2, silent = TRUE, facenum = 0.1*ncol(cut_2$it))
        
        # rgl::shade3d(cut_1, alpha = 0.1)
        # rgl::lines3d(target_3, col = "red")
        # rgl::points3d(rot_landmarks, size = 4, col = "blue")
        
        Morpho::mesh2ply(cut_1, gsub(".obj", "", new_beak_file, fixed = TRUE),
                         writeNormals = TRUE)
        
        Morpho::mesh2ply(cut_1, file.path(cut_2_folder, gsub(".obj", "", bird_file, fixed = TRUE)),
                         writeNormals = TRUE)
        
        if(!is.null(plot_folder)) {
          
          cut_1_plot_folder <- file.path(plot_folder, "cut_style_1")
          if(!dir.exists(cut_1_plot_folder)) {
            dir.create(cut_1_plot_folder)
          }
          
          cut_2_plot_folder <- file.path(plot_folder, "cut_style_2")
          if(!dir.exists(cut_2_plot_folder)) {
            dir.create(cut_2_plot_folder)
          }
          
          rgl::clear3d()
          rgl::shade3d(cut_1, col = "red")
          rgl::rgl.bringtotop()
          rgl::rgl.viewpoint(userMatrix = matrix(c(0.64, 0.71, -0.27, 0,
                                                   -0.32, 0.58, 0.74, 0,
                                                   0.69, -0.38, 0.61, 0,
                                                   0, 0, 0, 1), nrow = 4, byrow = TRUE),
                             zoom = 0.8)
          rgl::par3d(windowRect = c(20, 30, 800, 800))
          
          rgl::rgl.snapshot(file.path(cut_1_plot_folder, gsub(".obj", ".png", bird_file, fixed = TRUE)))
          
          rgl::clear3d()
          rgl::shade3d(cut_2, col = "red")
          rgl::rgl.bringtotop()
          rgl::rgl.viewpoint(userMatrix = matrix(c(0.64, 0.71, -0.27, 0,
                                                   -0.32, 0.58, 0.74, 0,
                                                   0.69, -0.38, 0.61, 0,
                                                   0, 0, 0, 1), nrow = 4, byrow = TRUE),
                             zoom = 0.8)
          rgl::par3d(windowRect = c(20, 30, 800, 800))
          
          rgl::rgl.snapshot(file.path(cut_2_plot_folder, gsub(".obj", ".png", bird_file, fixed = TRUE)))
          
        }
      
      
      } else {
        Morpho::mesh2ply(rot$mesh, gsub(".obj", "", new_beak_file, fixed = TRUE),
                         writeNormals = TRUE)
        
        if(!is.null(plot_folder)) {
          rgl::clear3d()
          rgl::shade3d(rot$mesh, col = "red")
          rgl::rgl.bringtotop()
          rgl::rgl.viewpoint(userMatrix = matrix(c(0.64, 0.71, -0.27, 0,
                                                   -0.32, 0.58, 0.74, 0,
                                                   0.69, -0.38, 0.61, 0,
                                                   0, 0, 0, 1), nrow = 4, byrow = TRUE),
                             zoom = 0.8)
          rgl::par3d(windowRect = c(20, 30, 800, 800))
          
          rgl::rgl.snapshot(file.path("beak_plots", gsub(".obj", ".png", bird_file, fixed = TRUE)))
        }
      }
    }
  }
  return(invisible(NULL))
}

for(i in seq_along(beak_files)) {
  reorient_beak(beak_files[i])
  cat("Finished reorienting beak", i, "of", length(beak_files), "->", basename(beak_files[i]), "\n")
}

########### do cuts and rotating ############

for(i in seq_along(beak_files)) {
  reorient_beak(beak_files[i], cut_extraneous = TRUE, mesh_folder = "data/ply_oriented_cut",
                landmark_folder = "data/landmarks_oriented", plot_folder = NULL)
  gc()
  cat("Finished reorienting beak", i, "of", length(beak_files), "->", basename(beak_files[i]), "\n")
}
