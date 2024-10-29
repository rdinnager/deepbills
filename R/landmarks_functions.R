
get_aligned_beak_landmarks <- function(beak_folder) {
  beak_dats <- list.files(beak_folder, full.names = TRUE)
  specs <- word(fs::path_ext_remove(basename(beak_dats)),
                end = 2L, sep = fixed("_"))
  landmarks <- map2(beak_dats, specs,
                    ~ fread(.x, col.names = c("x", "y", "z")) |>
                      mutate(Species = .y),
                    .progress = TRUE) |>
    list_rbind()
  landmarks
}

get_landmark_pf <- function(beak_folder) {
  landmarks <- get_aligned_beak_landmarks(beak_folder)
  
  landmark_df <- landmarks |>
    group_by(Species) |>
    mutate(point = 1:n()) |>
    ungroup() |>
    pivot_longer(c(-Species, -point), names_to = "coord", values_to = "value") |>
    mutate(point_name = paste0(coord, "_", point)) |>
    select(Species, value, point_name) |>
    pivot_wider(names_from = point_name, values_from = value)
  
  bird_beak_codes |>
    select(label, is_tip, phlo) |>
    left_join(landmark_df, by = c("label" = "Species"))
  
}