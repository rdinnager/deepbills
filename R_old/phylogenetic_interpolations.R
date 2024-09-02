library(readr)
library(plotly)
library(ape)
library(ggplot2)
library(ggtree)
library(reticulate)
library(rgl)
library(igraph)
library(future)
library(furrr)
library(stringr)
library(gifski)
library(imager)

big_latent_df <- readr::read_rds("data/big_latent_df_plotly.rds")
small_tree <- readr::read_rds("data/small_tree.rds")
dat <- readr::read_rds("data/latent_code_reconstructions.rds")

plot_mesh <- function(mesh) {
  
  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  
  fig <- plot_ly(
    x = mesh$x, y = mesh$y, z = mesh$z,
    i = mesh$i, j = mesh$j, k = mesh$k,
    facecolor = mesh$facecolor,
    type = "mesh3d"
  ) %>%
    layout(scene = list(aspectmode = "data",
                        xaxis = ax,
                        yaxis = ax,
                        zaxis = ax,
                        camera = list(eye = list(x = 1.5, y = 1.5, z = 1.5))))
  
  fig
  
}

small_tree_p <- ggtree::ggtree(small_tree)
small_tree_p

big_latent_df <- big_latent_df %>%
  dplyr::mutate(id = 1:dplyr::n())

metat <- small_tree_p$data %>%
  dplyr::left_join(dat %>%
                     dplyr::select(label = Binomal_Jetz,
                                   family = BLFamilyEnglish))

tips_to_interp <- sample(small_tree$tip.label, 2)
tips_to_interp <- c("Anhinga_anhinga", "Balaeniceps_rex")

make_phylo_interp <- function(tips_to_interp, big_latent_df, small_tree_p) {
  
  plot_name <- paste(tips_to_interp, collapse = "_to_")
  
  if(!dir.exists(file.path("plots", plot_name))) {
    dir.create(file.path("plots", plot_name))
  }
  
  latent_graph <- igraph::graph_from_data_frame(big_latent_df %>%
                                                  dplyr::select(node, parent) %>%
                                                  dplyr::distinct())
  
  nodes <- na.omit(c(metat$node[metat$label == tips_to_interp[1]],
             metat$node[metat$label == tips_to_interp[2]]))
  
  phylo_path <- igraph::shortest_paths(latent_graph, as.character(nodes[1]), 
                                       as.character(nodes[2]),
                                       mode = "all",
                                       output = "both")
  
  phylo_path_graph <- igraph::subgraph.edges(latent_graph, phylo_path$epath[[1]])
  root <- V(phylo_path_graph)[igraph::degree(phylo_path_graph, mode = "out") == 0]
  
  inward <- igraph::shortest_paths(phylo_path_graph, as.character(nodes[1]), 
                                   root,
                                   mode = "all",
                                   output = "both")
  
  outward <- igraph::shortest_paths(phylo_path_graph, root, 
                                    as.character(nodes[2]),
                                   mode = "all",
                                   output = "both")
  
  # edge <- inward$epath[[1]][1]
  get_edge_seq <- function(edge, inward = TRUE) {
    the_nodes <- igraph::ends(phylo_path_graph, edge) %>%
      as.vector()
    latent_seq <- big_latent_df %>%
      dplyr::filter(node == as.numeric(the_nodes[1]),
                    parent == as.numeric(the_nodes[2])) 
    if(inward){
      latent_seq <- latent_seq %>%
        dplyr::arrange(dplyr::desc(xs))
    } else {
      latent_seq <- latent_seq %>%
        dplyr::arrange(xs)
    }
  }
  
  latent_seq_df <- purrr::map_dfr(inward$epath[[1]],
                                  ~get_edge_seq(.x, inward = TRUE)) %>%
    dplyr::bind_rows(purrr::map_dfr(outward$epath[[1]],
                                    ~get_edge_seq(.x, inward = FALSE)))
  
  
  p2 <- small_tree_p +
    geom_point(data = metat %>%
               dplyr::filter(label %in% tips_to_interp),
             aes(x = x,
                 y = y,
                 colour = family)) +
    geom_text(data = metat %>%
                 dplyr::filter(label %in% tips_to_interp) %>%
                  dplyr::mutate(label = paste0("  ", label)),
               aes(x = x,
                   y = y,
                   label = label),
              hjust = 0) +
    xlim(0, 160) +
    theme(legend.position = "none")
  
  p2
  
  # row_num <- 1196
  make_phylo_interp_frame <- function(row_num) {
    mesh_dat <- latent_seq_df[row_num, ]
    mesh_plot <- plot_mesh(mesh_dat$recon_mesh[[1]])
    tree_plot <- p2 + geom_vline(xintercept = mesh_dat$xs, colour = "grey") +
      annotate("point", x = mesh_dat$xs, y = mesh_dat$ys,
               colour = "blue") 
    
    l1 <- plotly::ggplotly(tree_plot) %>%
      plotly::style(textposition = "right")
    l2 <- subplot(l1, mesh_plot, widths = c(0.5, 0.5)) %>%
      layout(showlegend = FALSE)
    l2
    
    orca(l2, file = file.path("plots", plot_name, paste0("plot_", 
                                                         stringr::str_pad(row_num,
                                                                          5,
                                                                          pad = "0"), 
                                                         ".png")),
         width = 1024, height = 728)
  }
  
  plan(multisession(workers = 7L))
  
  furrr::future_walk(seq_len(nrow(latent_seq_df)),
                     ~make_phylo_interp_frame(.x),
                     .progress = TRUE)
  
  
}

plan(multisession(workers = 7L))
make_phylo_interp(tips_to_interp, big_latent_df, small_tree_p)

make_gif <- function(png_folder, file_name, desired_len = 15, width = 1024, height = 728) {
  pngs <- list.files(png_folder, full.names = TRUE, pattern = ".png")
  delay <- desired_len / length(pngs)
  gifski::gifski(pngs, file_name, delay = delay, width = width, height = height)
}

make_gif("plots/Anhinga_anhinga_to_Balaeniceps_rex",
         file_name = "media/Anhinga_anhinga_to_Balaeniceps_rex.gif")


#### diversification animations

n_frames <- 30*10
rotate_num <- 3

angles <- seq(0, (2 * pi) * rotate_num, length.out = n_frames)
circle_coords <- purrr::map(angles,
                            ~ c(1.5 * cos(.x) + 1.5 * sin(.x),
                                -1.5 * sin(.x) + 1.5 * cos(.x))) %>%
  do.call(rbind, .)

all_xs <- unique(big_latent_df$xs)
all_xs <- all_xs[order(all_xs)]
x_seq <- c(all_xs[c(TRUE, FALSE)], all_xs[length(all_xs)])

plan(multisession(workers = 7L))

i <- 1
for(i in seq_along(x_seq)) {
  
  i_padded <- stringr::str_pad(i, 3, pad = "0")
  dir.create(file.path("plots/diversification", i_padded))
  
  dat <- big_latent_df %>%
    dplyr::filter(xs == x_seq[i]) %>%
    dplyr::arrange(dplyr::desc(ys))
  
  p2 <- small_tree_p +
    geom_vline(xintercept = dat$xs, colour = "grey") +
    geom_point(data = dat,
               aes(x = xs,
                   y = ys),
               colour = "darkgreen",
               size = 2) +
    theme(legend.position = "none")
  
  p2
  
  ggsave(file.path("plots/diversification", i_padded, "tree_plot.png"), p2,
         height = 6, width = 4)
  
  # j <- 1
  beak_pngs <- purrr::map_chr(seq_len(nrow(dat)), 
                              ~file.path("plots/diversification", i_padded, paste0("plot_", 
                                                               stringr::str_pad(.x,
                                                                                3,
                                                                                pad = "0"), 
                                                               ".png")))
  make_meshes <- function(j) {
    pm <- plot_mesh(dat$recon_mesh[[j]]) %>%
      layout(scene = list(camera = list(eye = list(x = circle_coords[i, 1], 
                                                   y = circle_coords[i, 2], 
                                                   z = 1.5))))
    
    
    orca(pm, file = beak_pngs[j],
         width = 512, height = 512)
    
      
  }
  
  furrr::future_walk(seq_len(nrow(dat)),
                     ~make_meshes(.x))
  
  beak_imgs <- purrr::map(beak_pngs,
                          ~imager::load.image(.x))
  
  img_app <- imager::imappend(beak_imgs, 'y')
  imager::save.image(img_app, file.path("plots/diversification", i_padded, "combined.png"))
  
  print(i)
  
}

test <- big_latent_df %>% filter(xs == big_latent_df$xs[30])

i <- 8
for(i in 1:300) {
  i_padded <- stringr::str_pad(i, 3, pad = "0")
  tree_p <- imager::load.image(file.path("plots/diversification", i_padded, "tree_plot.png"))
  beak_p <- imager::load.image(file.path("plots/diversification", i_padded, "combined.png")) %>%
    imager::flatten.alpha()
  
  if(height(beak_p) > height(tree_p)) {
    h <- height(tree_p)
    w <- ceiling(width(beak_p) * (height(tree_p) / height(beak_p)))
    beak_p <- imager::resize(beak_p, w, h, interpolation_type = 6)
  } else {
    beak_p <- imager::pad(beak_p, height(tree_p) - height(beak_p), 'y', val = c(1, 1, 1))
  }
  
  if(height(beak_p) < height(tree_p)) {
    beak_p <- imager::pad(beak_p, height(tree_p) - height(beak_p), 'y', val = c(1, 1, 1))
  }
  
  big_im <- imager::imappend(list(tree_p, beak_p), 'x')
  
  big_im <- imager::pad(big_im, 3000 - width(big_im), 'x', val = c(1, 1, 1))
  
  imager::save.image(big_im, file.path("plots/diversification/animation", 
                                       paste0("frame_", stringr::str_pad(i,
                                                                         3,
                                                                         pad = "0"),
                                              ".png")))
  
}

gifski(list.files("plots/diversification/animation", full.names = TRUE, pattern = ".png"), 
       "media/diversification_animation.gif",
       delay = 1/30,
       width = 3000, height = 1800)
