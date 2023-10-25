library(rcdd)
library(tessellation)
library(SimplicialCubature)

#' @importFrom rcdd makeH validcdd scdd q2d
#' @importFrom tessellation delaunay getDelaunaySimplicies
#' @importFrom SimplicialCubature adaptIntegrateSimplex

A <- rbind(
  c(-1, 0, 0), # -x
  c( 1, 0, 0), # x
  c( 0,-1, 0), # -y
  c( 1, 1, 0), # x+y
  c( 0, 0,-1), # -z
  c( 1, 1, 1)  # x+y+z
)
b <- c(5, 4, 5, 3, 10, 6)

A <- d2q(A)
b <- d2q(b)

H <- makeH(A, b)
. <- validcdd(H)

V <- q2d(scdd(H)[["output"]])
if(any(V[, 1L] != 0) || any(V[, 2L] != 1)) {
  stop("Invalid input.")
}
vertices <- V[, -c(1L, 2L)]

dlny <- delaunay(vertices)
simplices <- getDelaunaySimplicies(dlny)
nsimplices <- length(simplices)

# union of simplicies
d <- ncol(A)
S <- array(NA_real_, dim=c(d, d+1L, nsimplices))
for(i in seq_len(nsimplices)) {
  S[, , i] <- t(simplices[[i]])
}

# integration
f <- function(x, y, z) {
  x + y*z
}
stopifnot(length(formals(f)) == d)
g <- function(v) NULL
bdy <- sprintf(
  "f(%s)", paste0(sprintf("v[%d]", 1:d), collapse = ", ")
)
body(g) <- parse(text = bdy)

adaptIntegrateSimplex(g, S)


