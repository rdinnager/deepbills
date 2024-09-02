library(reticulate)
library(readr)
library(taxize)

meta <- readr::read_csv("data/resource.csv")

np <- import("numpy")

latent_codes <- np$load("C:/Users/rdinn/Desktop/latent_codes.npy")

genera <- unique(meta$Museum_genus)

genera_tsn <- taxize::get_tsn(genera)

res <- taxize::gnr_resolve(genera, data_source_ids = c(1, 4),
                           with_context = TRUE, best_match_only = TRUE)

not_hit <- genera[!genera %in% res$user_supplied_name]


uids <- taxize::get_uid(gsub("_", " ", meta$Binomal_Jetz), ask = TRUE)
readr::write_lines(uids, "data/uids_temp.txt")

colids <- taxize::get_colid(gsub("_", " ", meta$Binomal_Jetz), ask = FALSE)
readr::write_lines(uids, "data/uids_temp.txt")

eolids <- taxize::get_eolid(gsub("_", " ", meta$Binomal_Jetz), ask = FALSE)
readr::write_lines(eolids, "data/eolids_temp.txt")

not_found <- which(is.na(uids))

res <- taxize::gnr_resolve(gsub("_", " ", meta$Binomal_Jetz)[not_found], data_source_ids = c(1, 4),
                           resolve_once = TRUE, with_context = TRUE, canonical = TRUE)
res <- res[[1]]

classify <- taxize::classification(gsub("_", " ", meta$Binomal_Jetz),
                                   db = "ncbi")
