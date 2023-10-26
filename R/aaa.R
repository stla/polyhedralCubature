.onLoad <- function(lib, pkg) {
  loadNamespace("gmp") # to use as.character.bigq
  # Work around bug in code checking in R 4.2.2 for use of packages
  dummy <- function(x) gmp::as.bigq(x)
}
