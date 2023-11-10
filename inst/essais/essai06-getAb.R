library(ompr)
library(magrittr)

model <- MIPModel() %>%
  add_variable(x) %>% add_variable(y) %>% add_variable(z) %>%
  add_constraint(-5 <= x) %>% add_constraint(x <= 4) %>%
  add_constraint(-5 <= y) %>% add_constraint(y <= 3 - x) %>%
  add_constraint(-10 <= z) %>% add_constraint(z <= 6 - x - y)

ineqs <- extract_constraints(model)
signs <- ifelse(ineqs[["sense"]] == "<=", 1, -1)
A <- as.matrix(signs * ineqs[["matrix"]])
b <- signs * ineqs[["rhs"]]
