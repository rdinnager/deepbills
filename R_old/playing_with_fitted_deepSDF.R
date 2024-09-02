library(reticulate)
library(rgl)

current_wd <- getwd()
setwd(file.path(current_wd, "sdf"))

torch <- import("torch")
np <- import("numpy")
trimesh <- import("trimesh")

sdf_net <- import_from_path("sdf_net", "model")

latent_codes <- torch$load(sdf_net$LATENT_CODES_FILENAME)$to()

latent_codes$requires_grad <- FALSE

latent_codes$size()

SDFNet <- sdf_net$SDFNet()

SDFNet$load()
SDFNet$eval()

test_voxels <- SDFNet$get_voxels(latent_code = latent_codes[1], voxel_resolution = 64L)
test_voxels <- SDFNet$get_voxels(latent_code = torch$from_numpy(np$array(rnorm(64, sd = 0.01)))$float()$to('cuda:0'), voxel_resolution = 64L)

test_voxels = np$pad(test_voxels, 1L, mode='constant', constant_values=1)

size <- 2
voxel_resolution <- 64L
level <- 0

test <- measure$marching_cubes_lewiner(test_voxels, level=level, spacing=tuple(size / voxel_resolution, size / voxel_resolution, size / voxel_resolution))
names(test) <- c("vertices", "faces", "normals", "_")
test$vertices <- test$vertices - (size / 2)

#scatter3Drgl(test$vertices[,1], test$vertices[,2], test$vertices[,3])

test_mesh <- trimesh$Trimesh(vertices=test$vertices, faces=test$faces, vertex_normals=test$normals)

test_mesh2 <- tmesh3d(t(test_mesh$vertices), indices = t(test_mesh$faces + 1L), homogeneous = FALSE)

#test_mesh2 <- tmesh3d(t(test$vertices), indices = test_mesh$faces + 1L, normals = t(test$normals), homogeneous = FALSE)

rgl::shade3d(test_mesh2, col = "grey")

rgl::par3d(userMatrix = matrix(c(-0.8517141, -0.07263903,  0.5189332,    0,
                          0.2424379,  0.82332683,  0.5131612,    0,
                          -0.4645315,  0.56288362, -0.6836360,    0,
                          0.0000000,  0.00000000,  0.0000000,    1
), nrow = 4, byrow = TRUE))

rgl::par3d(windowRect = c(164, 187, 284, 302))

rgl::rgl.snapshot("../figures/test_thumbnail.png")
