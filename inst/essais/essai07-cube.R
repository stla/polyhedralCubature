library(polyhedralCubature)
library(ompr)

model <- MIPModel() %>%
  add_variable(x) %>% add_variable(y) %>% add_variable(z) %>%
  add_constraint(0 <= x) %>% add_constraint(x <= 1) %>%
  add_constraint(0 <= y) %>% add_constraint(y <= 1) %>%
  add_constraint(0 <= z) %>% add_constraint(z <= 1)
Ab <- getAb(model)
A <- Ab[["A"]]
b <- Ab[["b"]]

f <- function(x, y, z) {
  exp(x + y + z)
}
integrateOverPolyhedron(f, A, b)
