library(Morpho)
library(readobj)
library(pbapply)
## test my_mesh2obj
source("R/my_mesh2obj.R")

ply <- Morpho::ply2mesh("data/plys_reoriented_cut/Abeillia_abeillei_1.ply")
ply2 <- Morpho::ply2mesh("data/plys_reoriented_cut/Abroscopus_schisticeps_1.ply")

my_mesh2obj(ply, "data/objs_reoriented_cut/Abeillia_abeillei_1")

test <- readobj::read.obj("data/objs_reoriented_cut/Abeillia_abeillei_1.obj", convert.rgl = TRUE)[[1]]

test2 <- readobj::read.obj("data/objs/Abeillia_abeillei_1.obj", convert.rgl = TRUE)[[1]]

rgl::shade3d(test, col = "red")


file_names <- list.files("data/plys_reoriented_cut", full.names = TRUE)

save_objs <- function(file_name) {
  ply <- try(Morpho::ply2mesh(file_name))
  if(!inherits(ply, "try-error")){
    Rvcg::vcgObjWrite(ply, file.path("data/objs_reoriented_cut", gsub(".ply", "", basename(file_name), fixed = TRUE)))
  } else {
    print(basename(file_name))
    return(basename(file_name))
  }
  return(invisible(NULL))
}

pblapply(file_names, save_objs)