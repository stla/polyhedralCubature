library(rcdd)
library(tessellation)
library(SimplicialCubature)
library(spray)

#' @importFrom rcdd makeH validcdd scdd q2d
#' @importFrom tessellation delaunay getDelaunaySimplicies
#' @importFrom SimplicialCubature adaptIntegrateSimplex definePoly integrateSimplexPolynomial
#' @importFrom spray is.spray index coeffs


integrateOverPolyhedron <- function(f, A, b) {
  stopifnot(is.matrix(A))
  stopifnot(nrow(A) == length(b))
  # make H-representation
  H <- makeH(A, b)
  . <- validcdd(H)
  # make V-representation
  V <- scdd(H)[["output"]]
  if(is.character(V)) {
    V <- q2d(V)
  }
  if(any(V[, 1L] != 0) || any(V[, 2L] != 1)) {
    stop("Invalid input.")
  }
  # Delaunay
  vertices <- V[, -c(1L, 2L)]
  dlny <- delaunay(vertices)
  simplices <- getDelaunaySimplicies(dlny)
  nsimplices <- length(simplices)
  # union of the simplices
  d <- ncol(A)
  S <- array(NA_real_, dim=c(d, d+1L, nsimplices))
  for(i in seq_len(nsimplices)) {
    S[, , i] <- t(simplices[[i]])
  }
  # integrate
  if(is.function(f)) {
    if(length(formals(f)) != d) {
      stop("")
    }
    g <- function(v) NULL
    bdy <- sprintf(
      "f(%s)", paste0(sprintf("v[%d]", 1:d), collapse = ", ")
    )
    body(g) <- parse(text = bdy)
    adaptIntegrateSimplex(g, S)
  } else if(is.spray(f)) {
    P <- definePoly(coeffs(f), index(f))
    integrateSimplexPolynomial(P, S)
  }
}


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
f <- function(x, y, z) {
  x + y*z
}

integrateOverPolyhedron(f, A, b)

P <- f(lone(1, 3), lone(2, 3), lone(3, 3))
integrateOverPolyhedron(P, A, b)
