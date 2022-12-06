## Load your packages, e.g. library(targets).
source("./packages.R")
conflict_prefer("filter", "dplyr")

## Load your R files
lapply(list.files("./R", full.names = TRUE), source)

## tar_plan supports drake-style targets and also tar_target()
tar_plan(

# target = function_to_make(arg), ## drake style

# tar_target(target2, function_to_make2(arg)) ## targets style
  
  tar_target(bird_beak_avonet, bird_beak_codes %>%
               filter(is_tip) %>%
               mutate(Species3 = gsub(" ", "_", Scientific)) %>%
               left_join(avonet %>%
                           filter(is_tip) %>%
                           select(-phlo, -is_tip, -label), 
                         by = "Species3")),
  
  tar_target(beak_dat, bird_beak_avonet %>%
               select(Species3, starts_with("Beak.")) %>%
               drop_na()),
  
  tar_target(beak_pca, princomp(as.matrix(beak_dat %>% select(-Species3)))),
  
  tar_quarto(report, "doc/report.qmd")

)
