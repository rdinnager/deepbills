library(readr)
library(plotly)
library(ape)
library(ggplot2)
library(ggtree)
library(reticulate)
library(rgl)

#library(dash)
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
