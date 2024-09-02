library(Morpho)
library(Rvcg)
library(readobj)
library(rgl)

beak_files <- list.files("data/oriented_trimmed_meshes_obj_circle/", full.names = TRUE) 

beak_file <- beak_files[16]

beak_mesh <- try(readobj::read.obj(beak_file, convert.rgl = TRUE)[[1]])
rgl::shade3d(beak_mesh)

beak_mesh <- Rvcg::vcgClean(beak_mesh, sel = c(6, 1), tol = 0.02)

borders <- Rvcg::vcgBorder(beak_mesh)

rgl::shade3d(beak_mesh)
points3d(t(beak_mesh$vb[1:3,])[which(borders$bordervb == 1),],col="green")

border_pts <- t(beak_mesh$vb[1:3,])[which(borders$bordervb == 1),]

border_roll <- Rvcg::vcgBallPivoting(border_pts)

new_mesh <- Morpho::mergeMeshes(beak_mesh, border_roll)
shade3d(new_mesh, col = "red")


shade3d(border_roll)

test_fill <- Rvcg::vcgBallPivoting(beak_mesh, radius = 0.2)

rgl::shade3d(beak_mesh)
rgl::shade3d(test_fill)

segments3d(t(humface$vb[1:3,])[c(rbind(edges$vert1[edges$border == 1],edges$vert2[edges$border == 1])),],col=2,lwd=3)
