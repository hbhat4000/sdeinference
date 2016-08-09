library('nloptr')
eval_f <- function(x) {
  return(100 * (x[2] - x[1] * x[1])^2 + (1 - x[1])^2)
}

eval_grad_f <- function(x) {
  return(c(-400 * x[1] * (x[2] - x[1] * x[1]) - 2 * (1 - x[1]), 200 * (x[2] - x[1] * x[1]) ))
}

x0 <- c(-1.2, 1)

opts <- list("algorithm" = "NLOPT_LD_LBFGS", xtol_rel = "1.0e-6", "print_level" = 3, "check_derivatives" = TRUE, "check_derivatives_print" = "all")

res <- nloptr(x0 = x0, eval_f = eval_f, eval_grad_f = eval_grad_f, opts = opts)