
# LogUtility

<!-- badges: start -->
<!-- badges: end -->

The goal of LogUtility is to provide the exact solution and numerical solvers for the HJB PDE for maximizing log utility of
a portfolio comprised of a risky asset (Ito diffusion) and money-market account.

## Installation

You can install the latest version from GitHub via devtools

``` r
devtools::install_github("shill1729/LogUtility")
```

## Example

A basic example demonstrating numerical convergence for certain inputs (stability has not been analyzed yet for either scheme).

Note exact_sol will produce a warning about checking a vector in an if condition--it's harmless as long as you do not enter a negative wealth-level and I will fix it eventually.

```r
library(LogUtility)
mu <- 0.10
rate <- 0.05
volat <- 0.55
bankroll <- 1000
space_epsilon <- 0.001
time_length <- 0.5
time_res <- 500
space_res <- 500
space_length <- 2*bankroll-space_epsilon
x <- seq(space_epsilon, space_length, length.out = space_res+1)

num_sol <- explicit_integrator(bankroll, time_length, mu, rate, volat, "value", space_epsilon, time_res, space_res)
num_sol2 <- imex_integrator(bankroll, time_length, mu, rate, volat, "value", space_epsilon, time_res, space_res)
ex_sol <- exact_sol(bankroll, 0, time_length, mu, rate, volat)

num_control <- explicit_control(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res)
num_control2 <- imex_control(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res)
ex_control <- exact_control(mu, rate, volat)
sol <- data.frame(rbind(num_sol, num_sol2, ex_sol))
names(sol) <- c("value", "dv1", "dv2", "x")
print(sol)
print(data.frame(num_control, num_control2,  ex_control))

# Graphing solutions
u_ex <- explicit_integrator(bankroll, time_length, mu, rate, volat, "grid", space_epsilon, time_res, space_res)
u_imex <- imex_integrator(bankroll, time_length, mu, rate, volat, "grid", space_epsilon, time_res, space_res)
u_exact <- exact_sol(x, 0, time_length, mu, rate, volat)$v
plot(x, u_exact, type = "l", main = "Optimal Value")
lines(x, u_imex[time_res+1, ], col = "blue", lty = "dashed")
lines(x, u_ex[time_res+1, ], col = "red", lty = "dashed")
result <- data.frame(x, explicit = u_ex[time_res+1, ], imex = u_imex[time_res+1,], exact = u_exact)
print(head(result))
```

