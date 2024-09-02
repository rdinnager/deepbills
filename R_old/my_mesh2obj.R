my_mesh2obj <- function (x, filename = dataname, writeNormals = TRUE) 
{
  ismatrix <- FALSE
  x.it <- NULL
  dataname <- deparse(substitute(x))
  if (is.matrix(x)) {
    ismatrix <- TRUE
    dimsx <- dim(x)
    if (dimsx[2] == 3 && dimsx[1] != 3) 
      x <- t(x)
    x <- list(vb = x)
  }
  x.vb <- cbind("v", t(x$vb[1:3, ]))
  if(writeNormals) {
    if (!ismatrix) 
      x.normals <- cbind("vn", t(x$normals))
    if (!ismatrix) 
      n_col <- nrow(x$it)
      n_row <- ncol(x$it)
      x.it <- cbind("f", matrix(paste(t(x$it), t(x$it), sep = "//"), nrow = n_row, ncol = n_col))
  } else {
    if (!ismatrix) 
      x.it <- cbind("f", t(x$it))
  }
  if(writeNormals) {
    obj <- rbind(x.vb, x.normals, x.it)
  } else {
    obj <- rbind(x.vb, x.it)
  }
  Morpho:::write.obj(obj, filename = filename)
}