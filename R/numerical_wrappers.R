#' The numerical solution to HJB-PDE problem for optimal log-utility
#' 
#' @param bankroll the initial bankroll to invest with
#' @param time_length the time-horizon of the investment
#' @param mu the mean rate of return of the stock
#' @param rate the risk-free rate of the money-market account
#' @param volat the annual volatility of the stock
#' @param output string to choose output: 'value' or 'grid'
#' @param space_epsilon the truncation towards zero of the space (wealth) variable
#' @param time_res the time resolution for the time-grid
#' @param space_res the space resolution for the space-grid
#' 
#' @description {Numerically solve the Hamilton-Jacobi-Bellman PDE for the maximal log utility value of a portfolio on a
#' geometric Brownian motion and risk-free money-market account. An implicit-explicit scheme is used, where the
#' linear part of the differential equation is handled implicitly and the non-linear term is handled explicitly in time.}
#' @return vector
#' @export imex_integrator
imex_integrator <- function(bankroll, time_length,mu, rate, volat, output = "value", space_epsilon = 0.001, time_res = 1000, space_res = 400)
{
  if(output == "value")
  {
    return(imex_integrator1(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res))
  } else if(output == "grid")
  {
    return(imex_integrator_grid(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res))
  }
}

#' The numerical solution to HJB-PDE problem for optimal log-utility
#' 
#' @param bankroll the initial bankroll to invest with
#' @param time_length the time-horizon of the investment
#' @param mu the mean rate of return of the stock
#' @param rate the risk-free rate of the money-market account
#' @param volat the annual volatility of the stock
#' @param output string to choose output: 'value' or 'grid'
#' @param space_epsilon the truncation towards zero of the space (wealth) variable
#' @param time_res the time resolution for the time-grid
#' @param space_res the space resolution for the space-grid
#' 
#' @description {Numerically solve the Hamilton-Jacobi-Bellman PDE for the maximal log utility value of a portfolio on a
#' geometric Brownian motion and risk-free money-market account. An explicit scheme is used in the time variable for the 
#' entire discretized PDE, thus it stability is only available for small time-steps i.e. high time resolutions.}
#' @return vector
#' @export explicit_integrator
explicit_integrator <- function(bankroll, time_length,mu, rate, volat, output = "value", space_epsilon = 0.001,time_res = 1000, space_res = 400)
{
  if(output == "value")
  {
    return(explicit_integrator1(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res))
  } else if(output == "grid")
  {
    return(explicit_integrator_grid(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res))
  }
}
