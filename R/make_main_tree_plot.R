#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param phylo_w_rates_clads
make_main_tree_plot <- function(phylo_w_rates_clads, biggest_traj_clads, clust_pal, add_examples = TRUE) {
  
  examples <- biggest_traj_clads %>%
    ungroup() %>%
    mutate(id = seq_len(n()))
  
  size <- 0.5
  outline_size <- 1
  
  object <- phylo_w_rates_clads %>%
    mutate(traj_class = as.factor(traj_class))

  sel <- dplyr::select(object, traj_class, norm)
  
  ob <- phyf::pf_phyloflow(object)
  tree <- phyf::pf_as_phylo(ob)
  
  tree_df <- ggtree::fortify(tree)
  
  tree_df <- tree_df %>%
      dplyr::left_join(sel %>%
                         dplyr::mutate(label = phyf::pf_labels(ob)),
                       by = "label") %>%
    #drop_na(traj_class) %>%
    left_join(examples %>%
                select(label, id))
  
   p1 <- ggtree::ggtree(tree_df, size = 1.2, 
                        layout = "circular", ladderize = FALSE) + 
     ggtree::geom_tree(ggplot2::aes(color = traj_class,
                                    size = norm)) +
     scale_size_continuous(range = c(0.01, 0.8), trans = "exp") +
     scale_colour_manual(values = unname(clust_pal), na.translate = FALSE) +
     guides(color = guide_legend(ncol = 2))
   
   if(add_examples) {
     p1 <- p1 + 
     geom_point2(aes(x = branch, subset = !is.na(id)), fill = "black", size = 5,
                 shape = 21, colour = "white") +
     geom_text(aes(x = branch, label = id), color = "white", size = 3,
               fontface = "bold")
   }
  
   p1

}
