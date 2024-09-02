library(JuliaCall)
library(phyf)
library(dplyr)
library(ape)
julia_setup()
julia_library("PANDA")

julia_command('test = load("data/clads_Birds.jld2", "output")')

julia_command('test = load("data/clads_Birds.jld2", "output")')
julia_command('save_ClaDS_in_R(test, "data/clad_Birds.rData")')

load("data/clad_Birds.rData")
ape::write.tree(CladsOutput$tree, "data/clads_tree.tre")

tree <- CladsOutput$tree

tree <- ape::makeNodeLabel(tree)
rate_df <- tibble(label = c(tree$tip.label, tree$node.label)[tree$edge[ , 2]],
                  node_num = tree$edge[ , 2]) %>%
  mutate(div_rate = CladsOutput$lambdai_map,
         tip_rate = CladsOutput$lambdatip_map[node_num])

clads_pf <- pf_as_pf(CladsOutput$tree) %>%
  left_join(rate_df)

#tree <- drop.tip(CladsOutput$tree, which(!CladsOutput$tree$tip.label %in% bird_beak_codes$label[bird_beak_codes$is_tip]))

readr::write_rds(clads_pf, "data/clads_pf.rds")
