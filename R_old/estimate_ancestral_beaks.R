library(readr)
library(plotly)
library(ape)
library(ggplot2)
library(ggtree)
library(reticulate)
library(rgl)

dat <- readr::read_rds("data/latent_code_reconstructions.rds")
tree <- ape::read.tree("data/phylogenies/Stage2_MayrParSho_Ericson_set1_decisive.tre")
tree <- tree[[1]]
tree$tip.label

sum(dat$Binomal_Jetz %in% tree$tip.label)
tree_beaks <- ape::drop.tip(tree, which(!tree$tip.label %in% dat$Binomal_Jetz))

latent_codes <- do.call(rbind, dat$latent_code)
latent_codes_st <- scale(latent_codes)
rownames(latent_codes_st) <- dat$Binomal_Jetz

latent_codes_st <- latent_codes_st[tree_beaks$tip.label, ]

#mod <- RRphylo::RRphylo(tree = tree_beaks, latent_codes_st, clus = 0.8)

#readr::write_rds(mod, "data/latent_code_aces.rds")

p1 <- ggtree(tree_beaks, layout = "slanted")
p1

metat <- p1$data %>%
  dplyr::inner_join(dat %>%
                      dplyr::select(label = Binomal_Jetz,
                                    family = BLFamilyEnglish))

p2 <- p1 +
  geom_point(data = metat,
             aes(x = x,
                 y = y,
                 colour = family))

p2

l <- plotly::ggplotly(p2)
l

y_range <- 182:264
specs_filt <- metat$label[metat$y %in% y_range]
small_tree <- ape::keep.tip(tree_beaks, which(tree_beaks$tip.label %in% specs_filt))
plot(small_tree)

#readr::write_rds(small_tree, "data/small_tree.rds")

latent_small <- latent_codes_st[small_tree$tip.label, ]

mod_small <- RRphylo::RRphylo(tree = small_tree, y = latent_small, clus = 0.8)

library(RPANDA)

rates <- mod_small$rates
end_nodes <- rownames(rates) %>%
  as.numeric()
end_nodes[is.na(end_nodes)] <- sapply(rownames(rates)[is.na(end_nodes)],
                                      function(x) which(small_tree$tip.label == x))

rates_reorder <- rates[match(small_tree$edge[ , 2], end_nodes), ]

plot_ClaDS_phylo(small_tree, abs(rates_reorder), log = FALSE)

scales <- attr(scale(latent_codes), "scaled:scale")
centres <- attr(scale(latent_codes), "scaled:center")

aces_latent <- t((t(mod_small$aces) * scales) + centres)

aces_latent_list <- purrr::array_branch(aces_latent, 1)
latent_aces_df <- dplyr::tibble(node = as.integer(rownames(aces_latent))) %>%
  dplyr::mutate(latent_code = aces_latent_list)

latent_aces_df <- latent_aces_df %>%
  dplyr::bind_rows(dplyr::tibble(node = 1:length(small_tree$tip.label),
                                 Binomal_Jetz = small_tree$tip.label) %>%
                     dplyr::left_join(dat %>%
                                        dplyr::select(Binomal_Jetz,
                                                      latent_code)) %>%
                     dplyr::select(-Binomal_Jetz))


small_tree_p <- ggtree::ggtree(small_tree)
small_tree_p

tree_w_latent <- small_tree_p$data %>%
  dplyr::left_join(latent_aces_df)

x_range <- range(tree_w_latent$x)
depth <- max(x_range)
traversal_frames <- 600
#frames_per_len <- floor(traversal_frames / depth)
x_seq <- seq(x_range[1], x_range[2], length.out = traversal_frames)
end_node <- 1
make_edge_sequence <- function(end_node) {
  node_row <- tree_w_latent %>%
    dplyr::filter(node == end_node)
  parent_row <- tree_w_latent %>%
    dplyr::filter(node == node_row$parent[1])
  
  x_start <- parent_row$x[1]
  x_end <- node_row$x[1]
  edge_xs <- x_seq[x_seq > x_start & x_seq < x_end]
  n_frames <- length(edge_xs)
  latent_vector <- node_row$latent_code[[1]] - parent_row$latent_code[[1]]
  latent_codes_along <- purrr::map(seq(0, 1, length.out = n_frames),
                                   ~parent_row$latent_code[[1]] + 
                                     latent_vector * .x)
  
  dplyr::tibble(node = node_row$node[1], parent = node_row$parent[1],
                xs = edge_xs, ys = node_row$y,
                latent_codes = latent_codes_along)
}

big_latent_df <- purrr::map_dfr(tree_w_latent$node,
                                ~make_edge_sequence(.x))


rgl::shade3d(big_latent_df$recon_mesh[[200]])

source("R/sdf_tools.R")
setup_SDF()

big_latent_df <- big_latent_df %>%
  dplyr::mutate(recon_mesh = pbapply::pblapply(latent_codes, get_meshes_from_latent, voxel_res = 128L))

readr::write_rds(big_latent_df, "data/big_latent_df.rds")

make_plotly_mesh <- function(mesh) {
  
  deg2rad <- function(deg) {(deg * pi) / (180)}
  
  mesh <- Morpho::rotaxis3d(mesh, c(0, 1, 0), theta = deg2rad(90))
  mesh <- Morpho::rotaxis3d(mesh, c(1, 0, 0), theta = deg2rad(-90))
  
  x <- mesh$vb[1,]
  y <- mesh$vb[2,]
  z <- mesh$vb[3,]
  m <- matrix(c(x,y,z), ncol=3, dimnames=list(NULL,c("x","y","z")))
  
  zmean <- apply(t(mesh$it),MARGIN=1,function(row){mean(m[row,3])})
  
  facecolor <- colourvalues::color_values(zmean)
  
  i <- mesh$it[1,]-1 
  j <- mesh$it[2,]-1
  k <- mesh$it[3,]-1
  
  list(x = x, y = y, z = z, facecolor = facecolor,
       i = i, j = j, k= k)
}



library(readr)
library(plotly)
library(ape)
library(ggplot2)
library(ggtree)
library(reticulate)
library(rgl)

library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dashTable)
library(pbapply)

# big_latent_df <- readr::read_rds("data/big_latent_df.rds")
# big_latent_df <- big_latent_df %>%
#   dplyr::mutate(recon_mesh = pbapply::pblapply(recon_mesh,
#                                                make_plotly_mesh))
# 
# 
# 
# readr::write_rds(big_latent_df, "data/big_latent_df_plotly.rds")
big_latent_df <- readr::read_rds("data/big_latent_df_plotly.rds")
small_tree <- readr::read_rds("data/small_tree.rds")
dat <- readr::read_rds("data/latent_code_reconstructions.rds")

small_tree_p <- ggtree::ggtree(small_tree)
small_tree_p

big_latent_df <- big_latent_df %>%
  dplyr::mutate(id = 1:dplyr::n())

metat <- small_tree_p$data %>%
  dplyr::left_join(dat %>%
                      dplyr::select(label = Binomal_Jetz,
                                    family = BLFamilyEnglish))

p2 <- small_tree_p +
  geom_point(data = metat,
             aes(x = x,
                 y = y,
                 colour = family))

p2

l <- plotly::ggplotly(p2)
l

phy_plot <- l %>%
  add_markers(data = big_latent_df,
              x = ~xs,
              y = ~ys,
              size = 0.01,
              customdata = ~id,
              inherit = FALSE)

phy_plot

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
                        zaxis = ax))
  
  fig
  
}

default <- plot_mesh(big_latent_df$recon_mesh[[1]])


app <- Dash$new()

app$layout(
  htmlDiv(
    list(
      htmlDiv(list(
        dccGraph(
          id = 'beak-mesh-plot',
          figure = default
        )), style = list(
          width ='24%',
          display = 'inline-block',
          `padding-bottom` = '300px')
      ),
      htmlDiv(list(
        dccGraph(
          id = 'phy-plot',
          figure = phy_plot
        )), style = list(
          width = '60%',
          height = '800px',
          display = 'inline-block')
      )
    )
  )
)

app$callback(
  output = list(id='beak-mesh-plot', property='figure'),
  params = list(input(id='phy-plot', property='clickData')),
  function(hoverData) {
    id <- hoverData$points[[1]]$customdata
    mesh <- big_latent_df$recon_mesh[[which(big_latent_df[["id"]] == id)]]
    
    return(plot_mesh(mesh))
  }
)

app$run_server()


















vid_files <- list.files("E:/UC_comp_lab_backup/shapegan/images/",
                        full.names = TRUE)
library(gifski)

gifski::gifski(vid_files, width = 2160, height = 1080,
               delay = 1/30, gif_file = "animation.gif")

gifski::gifski(vid_files, width = 1200, height = 600,
               delay = 1/30, gif_file = "animation2")