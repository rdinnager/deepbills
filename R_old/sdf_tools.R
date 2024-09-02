setup_SDF <- function() {
  #current_wd <- getwd()
  #setwd(file.path(current_wd, "sdf"))
  
  reticulate::use_condaenv("tf-gpu", required = TRUE)
  
  torch <<- reticulate::import("torch")
  np <<- reticulate::import("numpy")
  trimesh <<- reticulate::import("trimesh")
  
  measure <<- reticulate::import("skimage.measure")
  
  sdf_net <<- reticulate::import_from_path("sdf_net", "model")
  
  SDFNet <<- sdf_net$SDFNet()
  
  SDFNet$load()
  SDFNet$eval()
  
  #setwd(current_wd)
}

get_meshes_from_latent <- function(latent_code, voxel_res = 64L, show = FALSE) {
  
  voxels <- SDFNet$get_voxels(latent_code = torch$from_numpy(np$array(latent_code))$float()$to('cuda:0'), voxel_resolution = voxel_res)
  
  voxels = np$pad(voxels, 1L, mode='constant', constant_values=1)
  
  size <- 2
  marched_cubes <- measure$marching_cubes_lewiner(voxels, level = 0, 
                                                  spacing = reticulate::tuple(size / voxel_res, 
                                                                              size / voxel_res, 
                                                                              size / voxel_res))
  
  names(marched_cubes) <- c("vertices", "faces", "normals", "_")
  marched_cubes$vertices <- marched_cubes$vertices - (size / 2)
  
  final_mesh <- rgl::tmesh3d(t(marched_cubes$vertices), 
                             indices = t(marched_cubes$faces + 1L), 
                             homogeneous = FALSE)
  
  if(show) {
    rgl::shade3d(final_mesh, col = "#d2b232")
  }
  
  return(final_mesh)
}