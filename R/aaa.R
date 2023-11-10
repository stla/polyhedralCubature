.onLoad <- function(lib, pkg) {
  loadNamespace("gmp")    # to use as.character.bigq
  loadNamespace("Matrix") # to use as.matrix
  # Work around bug in code checking in R 4.2.2 for use of packages
  dummy_gmp <- function(x) gmp::as.bigq(x)
  dummy_matrix <- function(x) Matrix::as.matrix(x)
}

#' @importFrom magrittr %>%
#' @export %>%
NULL
