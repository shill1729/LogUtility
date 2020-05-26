#' The exact solution to HJB-PDE problem for optimal log-utility
#' 
#' @param x the current wealth level
#' @param t the current time 
#' @param time_length the investment time-horizon
#' @param mu the mean rate of return of the stock
#' @param rate the risk-free rate of the money-market account
#' @param volat the annual volatility of the stock
#' 
#' @description {The exact solution to the Hamilton-Jacobi-Bellman PDE for the maximal log utility value of a portfolio on a
#' geometric Brownian motion and risk-free money-market account.}
#' @return vector containing the value, the first derivative, second derivative and initial bankroll
#' @export exact_sol
exact_sol <- function(x, t, time_length, mu, rate, volat)
{
  
  if(x == 0)
  {
    stop("argument 'x' cannot be zero")
  }
  lambda <- (mu - rate) / volat
  sol <- log(x) + (time_length - t)*(rate + 0.5*lambda^2)
  return(data.frame(v = sol, vx = 1/x, vxx = -1/(x^2), x = x))
}

#' The exact optimal control for the HJB PDE problem
#'
#' @param mu the mean return of the risky asset
#' @param rate the risk-free rate of return of the risk-free asset
#' @param volat the volatility of the risky asset
#' 
#' @description {The optimal control is the so-called Kelly fraction \eqn{(\mu-r)/\sigma}. }
#' @return numeric
#' @export exact_control
exact_control <- function(mu, rate, volat)
{
  return((mu-rate)/volat^2)
}
