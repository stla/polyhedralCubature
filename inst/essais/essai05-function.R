library(rcdd)
library(tessellation)
library(SimplicialCubature)
library(spray)
library(qspray)
library(gmp) # to load in .onLoad function

#' @importFrom rcdd makeH validcdd scdd q2d d2q qsum
#' @importFrom tessellation delaunay getDelaunaySimplicies
#' @importFrom SimplicialCubature adaptIntegrateSimplex definePoly integrateSimplexPolynomial
#' @importFrom spray is.spray index coeffs
#' @importFrom qspray integratePolynomialOnSimplex


integrateOverPolyhedron <- function(f, A, b) {
  stopifnot(is.matrix(A))
  stopifnot(nrow(A) == length(b))
  stopifnot(is.numeric(A) || is.character(A))
  mode_A <- mode(A)
  if(mode_A != mode(b)) {
    stop(
      "The matrix `A` and the vector `b` must have the same mode, ",
      "numeric or character."
    )
  }
  if(!is.element(mode_A, c("numeric", "character"))) {
    stop(
      "The matrix `A` and the vector `b` must be of mode ",
      "numeric or character."
    )
  }
  is_qspray <- inherits(f, "qspray")
  isnot_qspray <- is.function(f) || is.spray(f)
  if(!(is_qspray || isnot_qspray)) {
    stop("Invalid argument `f`.")
  }
  # make H-representation
  if(is_qspray) {
    if(!is.character(A)) {
      A <- d2q(A)
    }
    if(!is.character(b)) {
      b <- d2q(b)
    }
  }
  H <- makeH(A, b)
  . <- validcdd(H)
  # make V-representation
  V <- scdd(H)[["output"]]
  if(isnot_qspray && is.character(V)) {
    V <- q2d(V)
  }
  if(any(V[, 1L] != 0) || any(V[, 2L] != 1)) {
    stop("Invalid arguments `A` and/or `b`.")
  }
  # Delaunay
  vertices <- V[, -c(1L, 2L)]
  if(is_qspray) {
    dlny <- delaunay(q2d(vertices))
  } else {
    dlny <- delaunay(vertices)
  }
  simplices <- getDelaunaySimplicies(dlny)
  nsimplices <- length(simplices)
  if(isnot_qspray) {
    # union of the simplices
    d <- ncol(A)
    U <- array(NA_real_, dim=c(d, d+1L, nsimplices))
    for(i in seq_len(nsimplices)) {
      U[, , i] <- t(simplices[[i]])
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
      adaptIntegrateSimplex(g, U)
    } else if(is.spray(f)) {
      P <- definePoly(coeffs(f), index(f))
      integrateSimplexPolynomial(P, U)
    }
  } else { # qspray
    results <- character(nsimplices)
    for(i in seq_len(nsimplices)) {
      S <- simplices[[i]]
      Svs <- as.integer(rownames(S))
      qS <- vertices[Svs, ]
      results[i] <- as.character(integratePolynomialOnSimplex(f, qS))
    }
    qsum(results)
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

# function
integrateOverPolyhedron(f, A, b)

# spray
P <- f(lone(1, 3), lone(2, 3), lone(3, 3))
integrateOverPolyhedron(P, A, b)

# qspray
Q <- f(qlone(1), qlone(2), qlone(3))
integrateOverPolyhedron(Q, A, b)
