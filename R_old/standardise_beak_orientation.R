library(Morpho)
library(readobj)
library(rgl)

test_mesh <- readobj::read.obj("data/objs/Abeillia_abeillei_1.obj", convert.rgl = TRUE)[[1]]
test_landmarks <- readr::read_delim("data/landmarks/Markup_1/Abeillia_abeillei_1.txt", " ",
                                    col_names = c("X", "Y", "Z"))
main_landmarks <- test_landmarks[1:4, ] %>%
  as.matrix()

rgl::shade3d(test_mesh, alpha = 0.1)
rgl::points3d(main_landmarks, size = 5, col = "red")

mid_beak <- apply(main_landmarks[3:4, ], 2, mean) %>% matrix(nrow = 1)
line_1 <- rbind(mid_beak, main_landmarks[1, ])

rgl::shade3d(test_mesh, alpha = 0.1)
rgl::points3d(main_landmarks, size = 5, col = "red")
rgl::points3d(mid_beak, size = 5, col = "green")
rgl::lines3d(line_1, col = "green")

line_1_vec <- (line_1[1, , drop = FALSE] - line_1[2, , drop = FALSE]) / dist(line_1)

top_line <- main_landmarks[2, ] - mid_beak

tt <- top_line %*% t(line_1_vec)

point_on_line_1 <- mid_beak + as.vector(tt) * line_1_vec

rgl::shade3d(test_mesh, alpha = 0.1)
rgl::points3d(main_landmarks, size = 5, col = "red")
rgl::points3d(mid_beak, size = 5, col = "green")
rgl::lines3d(line_1, col = "green")
rgl::points3d(point_on_line_1, size = 5, col = "blue")

landmarks_3 <- rbind(main_landmarks[1, ], point_on_line_1, main_landmarks[2, ])
landmark_dists <- dist(landmarks_3) %>% as.matrix()
target_3 <- matrix(c(landmark_dists[1, 2], 0, 0,
                     0, 0, 0,
                     0, 0, landmark_dists[2, 3]),
                   nrow = 3, byrow = TRUE)

rgl::shade3d(test_mesh, alpha = 0.1)
rgl::points3d(main_landmarks, size = 5, col = "red")
rgl::points3d(mid_beak, size = 5, col = "green")
rgl::lines3d(line_1, col = "green")
rgl::points3d(point_on_line_1, size = 5, col = "blue")
rgl::lines3d(target_3, col = "red")

test_rot <- Morpho::rotmesh.onto(test_mesh, landmarks_3, target_3, adnormals = TRUE)

rgl::shade3d(test_rot$mesh, alpha = 0.1)
rgl::lines3d(target_3, col = "red")

rgl::shade3d(test_rot$mesh, col = "red")
rgl::rgl.bringtotop()
rgl::rgl.viewpoint(userMatrix = matrix(c(0.64, 0.71, -0.27, 0,
                                         -0.32, 0.58, 0.74, 0,
                                         0.69, -0.38, 0.61, 0,
                                         0, 0, 0, 1), nrow = 4, byrow = TRUE),
                   zoom = 0.8)
rgl::par3d(windowRect = c(20, 30, 800, 800))

rgl::rgl.snapshot("test.png")

beak_files <- list.files("data/objs/", full.names = TRUE) 

reorient_beak <- function(beak_file) {
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
  
  rot <- Morpho::rotmesh.onto(beak_mesh, landmarks_3, target_3, adnormals = TRUE)
  
  rgl::shade3d(rot$mesh, col = "red")
  rgl::rgl.bringtotop()
  rgl::rgl.viewpoint(userMatrix = matrix(c(0.64, 0.71, -0.27, 0,
                                           -0.32, 0.58, 0.74, 0,
                                           0.69, -0.38, 0.61, 0,
                                           0, 0, 0, 1), nrow = 4, byrow = TRUE),
                     zoom = 0.8)
  rgl::par3d(windowRect = c(20, 30, 800, 800))
  
  rgl::rgl.snapshot("test.png")
}

