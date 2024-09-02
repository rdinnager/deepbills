library(readr)
library(plotly)
library(colourvalues)

library(dash)
library(dashHtmlComponents)
library(dashCoreComponents)
library(dashTable)

dat <- readr::read_rds("data/latent_code_reconstructions.rds")

deg2rad <- function(deg) {(deg * pi) / (180)}

mesh <- dat$recon_mesh[[2]]
species <- dat$Binomal_Jetz[2]

plot_mesh <- function(mesh, species) {
  
  ax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    showticklabels = FALSE,
    showgrid = FALSE
  )
  
  mesh <- Morpho::rotaxis3d(mesh, c(0, 1, 0), theta = deg2rad(90))
  mesh <- Morpho::rotaxis3d(mesh, c(1, 0, 0), theta = deg2rad(-90))
  
  x <- mesh$vb[1,]
  y <- mesh$vb[2,]
  z <- mesh$vb[3,]
  m <- matrix(c(x,y,z), ncol=3, dimnames=list(NULL,c("x","y","z")))
  
  zmean <- apply(t(mesh$it),MARGIN=1,function(row){mean(m[row,3])})
  
  facecolor <- colourvalues::color_values(zmean)
  
  fig <- plot_ly(
    x = x, y = y, z = z,
    i = mesh$it[1,]-1, j = mesh$it[2,]-1, k = mesh$it[3,]-1,
    facecolor = facecolor,
    type = "mesh3d"
  ) %>%
    layout(scene = list(aspectmode = "data",
                        xaxis = ax,
                        yaxis = ax,
                        zaxis = ax),
           title = list(text = gsub("_", " ", species)))
  
  fig
  
}


# fig <- plot_ly(
#   x = x, y = y, z = z,
#   i = mesh$it[1,]-1, j = mesh$it[2,]-1, k = mesh$it[3,]-1,
#   facecolor = rep("#d2b232", length(zmean)),
#   type = "mesh3d"
# ) %>%
#   layout(scene = list(aspectmode = "data"))
# 
# 
# fig

defaultPlot <- plot_mesh(mesh, species)
defaultPlot


beak_umap <- plot_ly(dat, 
                     x = ~X0, 
                     y = ~Y0, 
                     z = ~Z0,
                     text = ~Binomal_Jetz,
                     color = ~BLFamilyEnglish, 
                     colors = hues::iwanthue(n = dplyr::n_distinct(dat$BLFamilyEnglish))) %>%
  add_markers() %>%
  layout(scene = list(aspectmode = "data"),
         height = 800)

beak_umap
  

app <- Dash$new()

app$layout(
  htmlDiv(
    list(
      htmlDiv(list(
        dccGraph(
          id = 'beak-mesh-plot',
          figure = defaultPlot
        )), style = list(
          width ='24%',
          display = 'inline-block',
          `padding-bottom` = '300px')
      ),
      htmlDiv(list(
        dccGraph(
          id = 'beak-umap',
          figure = beak_umap
        )), style = list(
          width = '60%',
          height = '800px',
          display = 'inline-block')
      )
    )
  )
)

app$callback(
  output = list(output(id='beak-mesh-plot', property='figure')),
  params = list(input(id='beak-umap', property='clickData')),
  function(hoverData) {
    spec <- hoverData$points[[1]]$text
    print(spec)
    mesh <- dat$recon_mesh[[which(dat[["Binomal_Jetz"]] == spec)]]
    print(mesh)
    
    return(list(plot_mesh(mesh, spec)))
  }
)

app$run_server()


