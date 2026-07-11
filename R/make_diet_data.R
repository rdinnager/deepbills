#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bird_beak_codes
#' @param diet_data
#' @return
#' @author rdinnager
#' @export
make_diet_data <- function(bird_beak_avonet, diet_data) {

  diet_specs <- unique(diet_data$Scientific_Name)
  diet_genera <- unique(word(diet_specs, 1))
  sum(diet_specs %in% bird_beak_avonet$Scientific)
  sum(diet_specs %in% bird_beak_avonet$Species3)
  
  sum(diet_genera %in% word(bird_beak_avonet$Species3, 1))
  
  # genera_latent <- bird_beak_avonet |>
  #   left_join(trophic_niche_dat_pca_all, by = c("label" = "Species")) |>
  #   mutate(Genus = word(Species3, 1))
  
  diet <- diet_data |>
    #filter(Prey_Kingdom == "Animalia") |>
    mutate(Prey_Part = str_split(Prey_Part, ";")) |>
    mutate(Prey_Part_N = map_int(Prey_Part, ~length(.x))) |>
    unnest_longer(Prey_Part) |>
    mutate(Prey_Part = str_trim(Prey_Part)) |>
    mutate(Fraction_Diet = Fraction_Diet / Prey_Part_N) |>
    group_by(Scientific_Name, Prey_Kingdom, Prey_Class, Prey_Part) |>
    summarise(Fraction_Diet = sum(Fraction_Diet, na.rm = TRUE)) |>
    ungroup() |>
    drop_na(Fraction_Diet, Prey_Class) |>
    mutate(Prey_Class = ifelse(is.na(Prey_Part), Prey_Class,  Prey_Part)) |>
    mutate(Prey_Class = ifelse(Prey_Kingdom == "Plantae" & is.na(Prey_Part), NA_character_, Prey_Class)) |>
    mutate(Prey_Class = ifelse(Prey_Kingdom %in% c("Fungi", "Chromista", "Protozoa", "Bacteria"), "Other", Prey_Class)) |>
    group_by(Scientific_Name, Prey_Kingdom, Prey_Class, Prey_Part) |>
    summarise(Fraction_Diet = sum(Fraction_Diet, na.rm = TRUE)) |>
    group_by(Scientific_Name) |>
    mutate(Fraction_Diet = Fraction_Diet / sum(Fraction_Diet, na.rm = TRUE)) |>
    mutate(plant_items = sum(Prey_Kingdom == "Plantae"), 
           plant_total = sum(Fraction_Diet[Prey_Kingdom == "Plantae"])) |>
    filter(Prey_Kingdom != "Plantae" | (plant_items > 1 & !is.na(Prey_Class))) |>
    mutate(Fraction_Diet = ifelse(Prey_Kingdom == "Plantae", 
                                  Fraction_Diet/sum(Fraction_Diet[Prey_Kingdom == "Plantae"]), 
                                  Fraction_Diet)) |>
    mutate(Fraction_Diet = ifelse(Prey_Kingdom == "Plantae", 
                                  Fraction_Diet/plant_total, Fraction_Diet)) |>
    mutate(Fraction_Diet = Fraction_Diet / sum(Fraction_Diet, na.rm = TRUE))
 
  
  genera_diet <- diet |>
    mutate(Genus = word(Scientific_Name, 1)) |>
    group_by(Prey_Class) |>
    mutate(drop = sum(Fraction_Diet > 0) < 3) |>
    filter(!drop) |>
    group_by(Genus, Prey_Class) |>
    summarise(Fraction_Diet = sum(Fraction_Diet, na.rm = TRUE)) |>
    group_by(Genus) |>
    mutate(Fraction_Diet = Fraction_Diet / sum(Fraction_Diet, na.rm = TRUE)) |>
    pivot_wider(names_from = Prey_Class, values_from = Fraction_Diet, values_fill = 0)
  
  colnames(genera_diet)[-1] <- paste0("diet_", colnames(genera_diet)[-1])
  
  genera_diet
    

}
