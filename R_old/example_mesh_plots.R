library(rgl)
library(readobj)
library(Rvcg)

mesh <- read.obj("data/oriented_scaled_meshes_obj/Balaeniceps_rex_2.obj",
                 convert.rgl = TRUE)


mesh <- read.obj("data/objs/Balaeniceps_rex_2.obj",
                 convert.rgl = TRUE)

wire3d(mesh)

mesh2 <- Rvcg::vcgQEdecim(mesh[[1]], 2000L)

wire3d(mesh2)