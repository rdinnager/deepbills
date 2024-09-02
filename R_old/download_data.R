library(readr)
library(dplyr)
library(tidyr)
library(purrr)

resources <- read_csv("D:/Data/Mark_My_Bird/resource.csv")

## fix errors in some filenames
resources <- resources %>%
  dplyr::mutate(new_download = paste0("http://lfs.nhm.ac.uk/mark-my-bird/",
                                      Binomal_Jetz, "/", ID, ".obj"))
  # dplyr::mutate(MMB_scan_download = 
  #                 purrr::map2_chr(ID, MMB_scan_download, ~gsub("[^/]*$", paste0(.x, ".obj"), .y)))
                  
long_data <- resources %>%
  dplyr::select(ID, Binomal_Jetz, starts_with("Markup")) %>%
  gather(markup_num, markup_label, -ID, -Binomal_Jetz) %>%
  mutate(url = paste0("http://lfs.nhm.ac.uk/mark-my-bird/", Binomal_Jetz, "/", 
                      ID, "_", markup_label, ".txt"),
         file_name = paste0("data/landmarks/", markup_num, "/", 
                            ID, ".txt"))

for(i in 1:length(resources$new_download)) {
  local_file <- file.path("data/objs", basename(resources$new_download[i]))
  if(!file.exists(local_file)) {
    try(download.file(resources$new_download[i], 
                      local_file, 
                      quiet = FALSE))
    message("Downloaded File # ", i, "of ", nrow(resources), ": ", basename(resources$new_download[i]))
  }
  Sys.sleep(1)
}


for(i in 1:length(long_data$url)) {
  if(!dir.exists(paste0("data/landmarks/", long_data$markup_num[i]))) {
    dir.create(paste0("data/landmarks/", long_data$markup_num[i]))
  }
  if(!file.exists(long_data$file_name[i])){
    try(download.file(long_data$url[i], long_data$file_name[i], quiet = TRUE))
    message("Downloaded File # ", i, "of ", nrow(long_data), ": ", basename(long_data$file_name[i]))
    Sys.sleep(0.2)
  }
}

