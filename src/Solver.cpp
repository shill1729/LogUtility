// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]


std::vector<double> discretize_space(double space_length, double space_epsilon, int space_res)
{
  // The price space grid is x_0, ..., x_M, so M+1 points, partitioning [epsilon,B] into a grid.
  double space_step = (space_length-space_epsilon) / space_res;
  std::vector<double> space_grid(space_res + 1);
  for (int i = 0; i < space_res + 1; ++i)
  {
    space_grid[i] = space_epsilon + i * space_step;
  }
  return space_grid;
}

// Thomas algorithm for solving Ax=d when A = triDiag(a, b, c). Stability only for diagonally dominant or symmetric positive definite.
std::vector<double> thomas_algorithm(std::vector<double> a, std::vector<double> b, std::vector<double> c, std::vector<double> d)
{
  // On paper, b and d are from 1 to n, the vector a, the lower diag is from 2 to n, and c, the upper diag, is from 1 to n-1.
  int n = d.size();
  // To preserve the coefficients, we make new ones for the forward sweep and x for the solution to be returned.
  std::vector<double> cc(n - 1);
  std::vector<double> dd(n);
  std::vector<double> x(n);
  cc[0] = c[0] / b[0];
  // Forward sweep for c vector
  for (int i = 1; i < n - 1; ++i)
  {
    cc[i] = c[i] / (b[i] - a[i] * cc[i - 1]);
  }
  // Repeating nearly the same for d
  dd[0] = d[0] / b[0];
  for (int i = 1; i < n; ++i)
  {
    // On paper a_i is written indexed from 2 to n, but since the values in the computer are indexed starting from 0, to n-1, and we are looping from 1 to n-1, we must decrement a_i to a_{i-1}
    // It actually doesn't matter for a_i=a constants
    dd[i] = (d[i] - a[i-1] * dd[i - 1]) / (b[i] - a[i-1] * cc[i - 1]);
  }
  // Back substitution for the solution
  x[n - 1] = dd[n - 1];
  for (int i = n - 2; i >= 0; --i)
  {
    x[i] = dd[i] - cc[i] * x[i + 1];
  }
  return x;
}

//' The numerical solution to HJB-PDE problem for optimal log-utility
//' 
//' @param bankroll the initial bankroll to invest with
//' @param time_length the time-horizon of the investment
//' @param mu the mean rate of return of the stock
//' @param rate the risk-free rate of the money-market account
//' @param volat the annual volatility of the stock
//' @param space_epsilon the truncation towards zero of the space (wealth) variable
//' @param time_res the time resolution for the time-grid
//' @param space_res the space resolution for the space-grid
//' 
//' @description {Numerically solve to the Hamilton-Jacobi-Bellman PDE for the maximal log utility value of a portfolio on a
//' geometric Brownian motion and risk-free money-market account.}
//' @return matrix
// [[Rcpp::export]]
Eigen::MatrixXd explicit_integrator_grid(double bankroll, double time_length, double mu, double rate, double volat, double space_epsilon = 0.001, int time_res = 1000, int space_res = 400)
{
  // Set length of space to twice the bankroll
  double space_length = 2*bankroll-space_epsilon;
  // Get space-step and time-step
  double space_step = (space_length - space_epsilon) / space_res;
  double time_step = (time_length) / time_res;
  // Create space-grid
  std::vector<double> x = discretize_space(space_length, space_epsilon, space_res);
  // Compute market price of risk
  double lambda = (mu - rate) / volat;
  // Initialize solution matrix
  Eigen::MatrixXd solution(time_res + 1, space_res + 1);
  solution = solution.Zero(time_res + 1, space_res + 1);
  // Set initial conditions
  for (int j = 0; j < space_res + 1; ++j)
  {
    solution(0, j) = std::log(x[j]);
  }
  // Set boundary conditions
  for (int i = 0; i < time_res + 1; ++i)
  {
    solution(i, 0) = std::log(x[0]);
    solution(i, space_res) = std::log(x[space_res]);
  }
  // Time-stepping explicit integrator
  for (int i = 0; i < time_res; ++i)
  {
    for (int j = 1; j < space_res; ++j)
    {
      // Central difference approx. to first derivative of value function
      double fd1 = (solution(i, j + 1) - solution(i, j - 1)) / (2 * space_step);
      // Non-linear term: square of first derivative divided by second derivative
      double nlq = std::pow((solution(i, j + 1) - solution(i, j - 1)), 2.0) / (4 * (solution(i, j + 1) - 2 * solution(i, j) + solution(i, j - 1)));
      // Time-stepping iteration for the explicit integrator
      solution(i + 1, j) = solution(i, j) + time_step * (rate*x[j] * fd1 - 0.5*std::pow(lambda, 2.0)*nlq);
    }
  }
  return solution;
  
}



//' The numerical solution to HJB-PDE problem for optimal log-utility
//' 
//' @param bankroll the initial bankroll to invest with
//' @param time_length the time-horizon of the investment
//' @param mu the mean rate of return of the stock
//' @param rate the risk-free rate of the money-market account
//' @param volat the annual volatility of the stock
//' @param space_epsilon the truncation towards zero of the space (wealth) variable
//' @param time_res the time resolution for the time-grid
//' @param space_res the space resolution for the space-grid
//' 
//' @description {Numerically solve to the Hamilton-Jacobi-Bellman PDE for the maximal log utility value of a portfolio on a
//' geometric Brownian motion and risk-free money-market account.}
//' @return matrix
// [[Rcpp::export]]
Eigen::MatrixXd imex_integrator_grid(double bankroll, double time_length, double mu, double rate, double volat, double space_epsilon = 0.001, int time_res = 1000, int space_res = 400)
{
  // Set length of space to twice the bankroll
  double space_length = 2*bankroll-space_epsilon;
  // Get space-step and time-step
  double space_step = (space_length - space_epsilon) / space_res;
  double time_step = (time_length) / time_res;
  // Create space-grid
  std::vector<double> x = discretize_space(space_length, space_epsilon, space_res);
  // Compute market price of risk
  double lambda = (mu - rate) / volat;
  
  
  
  // Initialize solution matrix
  Eigen::MatrixXd solution(time_res + 1, space_res + 1);
  solution = solution.Zero(time_res + 1, space_res + 1);
  // Set initial conditions
  for (int j = 0; j < space_res + 1; ++j)
  {
    solution(0, j) = std::log(x[j]);
  }
  // Set boundary conditions
  for (int i = 0; i < time_res + 1; ++i)
  {
    solution(i, 0) = std::log(x[0]);
    solution(i, space_res) = std::log(x[space_res]);
  }
  // Auxillary vector
  std::vector<double> d(space_res - 1);
  // First entry is alpha*LB, last entry is delta*UB, all else zero
  d[0] = time_step * rate * x[1] /(2 * space_step) * std::log(x[0]);
  d[space_res - 2] = - time_step * rate * x[space_res - 1] / (2 * space_step) * std::log(x[space_res]);
  for (int i = 1; i < space_res - 2; ++i)
  {
    d[i] = 0;
  }

  // First we set up the matrix coefficients
  std::vector<double> alphaV(space_res - 2);
  std::vector<double> betaV(space_res - 1);
  std::vector<double> deltaV(space_res - 2);
  for (int j = 0; j < space_res - 2; ++j)
  {
    alphaV[j] = time_step*rate*x[j + 2]/(2*space_step);
    deltaV[j] = -time_step*rate*x[j + 1]/(2*space_step);
  }
  for (int j = 0; j < space_res - 1; ++j)
  {
    betaV[j] = 1.0;
  }
  
  
  // Time-stepping implicit-explicit integrator
  for (int i = 0; i < time_res; ++i)
  {
    // The recursive RHS vector of the solution at the previous time step
    std::vector<double> bb(space_res - 1);
    for (int j = 1; j < space_res; ++j)
    {
      // Non-linear term: square of first derivative divided by second derivative
      double nlq = std::pow((solution(i, j + 1) - solution(i, j - 1)), 2.0) / (4 * (solution(i, j + 1) - 2 * solution(i, j) + solution(i, j - 1)));
      bb[j - 1] = solution(i, j)-time_step*0.5*std::pow(lambda, 2.0)*nlq-d[j - 1];
    }
    // Solve the linear equation using any method; Thomas algorithm works best for tridiagonal systems
    std::vector<double> sol = thomas_algorithm(alphaV, betaV, deltaV, bb);
    for (int j = 1; j < space_res; ++j)
    {
    
      // Time-stepping iteration for the explicit integrator
      solution(i + 1, j) = sol[j-1];
    }
  }
  return solution;
}

//' The numerical solution to HJB-PDE problem for optimal log-utility
//' 
//' @param bankroll the initial bankroll to invest with
//' @param time_length the time-horizon of the investment
//' @param mu the mean rate of return of the stock
//' @param rate the risk-free rate of the money-market account
//' @param volat the annual volatility of the stock
//' @param space_epsilon the truncation towards zero of the space (wealth) variable
//' @param time_res the time resolution for the time-grid
//' @param space_res the space resolution for the space-grid
//' 
//' @description {Numerically solve to the Hamilton-Jacobi-Bellman PDE for the maximal log utility value of a portfolio on a
//' geometric Brownian motion and risk-free money-market account.}
//' @return vector
// [[Rcpp::export]]
std::vector<double> explicit_integrator1(double bankroll, double time_length, double mu, double rate, double volat, double space_epsilon = 0.001, int time_res = 1000, int space_res = 400)
{
  // Set length of space to twice the bankroll
  double space_length = 2*bankroll-space_epsilon;
  // Get space-step and time-step
  double space_step = (space_length - space_epsilon) / space_res;
  Eigen::MatrixXd solution = explicit_integrator_grid(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res);
  
  std::vector<double> output(4);
  int index = space_res/2;
  double value = solution(time_res, index);
  double vx = (solution(time_res, index + 1) - solution(time_res, index - 1)) / (2 * space_step);
  double vxx = (solution(time_res, index + 1) - 2 * solution(time_res, index) + solution(time_res, index - 1))/(std::pow(space_step, 2.0));
  output[0] = value;
  output[1] = vx;
  output[2] = vxx;
  output[3] = bankroll;
  return output;
}

//' The numerical solution to HJB-PDE problem for optimal log-utility
//' 
//' @param bankroll the initial bankroll to invest with
//' @param time_length the time-horizon of the investment
//' @param mu the mean rate of return of the stock
//' @param rate the risk-free rate of the money-market account
//' @param volat the annual volatility of the stock
//' @param space_epsilon the truncation towards zero of the space (wealth) variable
//' @param time_res the time resolution for the time-grid
//' @param space_res the space resolution for the space-grid
//' 
//' @description {Numerically solve to the Hamilton-Jacobi-Bellman PDE for the maximal log utility value of a portfolio on a
//' geometric Brownian motion and risk-free money-market account.}
//' @return vector
// [[Rcpp::export]]
std::vector<double> imex_integrator1(double bankroll, double time_length, double mu, double rate, double volat, double space_epsilon = 0.001, int time_res = 1000, int space_res = 400)
{
  // Set length of space to twice the bankroll
  double space_length = 2*bankroll-space_epsilon;
  // Get space-step and time-step
  double space_step = (space_length - space_epsilon) / space_res;
  Eigen::MatrixXd solution = imex_integrator_grid(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res);
  
  std::vector<double> output(4);
  int index = space_res/2;
  double value = solution(time_res, index);
  double vx = (solution(time_res, index + 1) - solution(time_res, index - 1)) / (2 * space_step);
  double vxx = (solution(time_res, index + 1) - 2 * solution(time_res, index) + solution(time_res, index - 1))/(std::pow(space_step, 2.0));
  output[0] = value;
  output[1] = vx;
  output[2] = vxx;
  output[3] = bankroll;
  return output;
}


//' The numerical optimal control value
//' 
//' @param bankroll the initial wealth level
//' @param time_length the time-horizon of the investment
//' @param mu the mean rate of return of the stock
//' @param rate the risk-free rate of the money-market account
//' @param volat the annual volatility of the stock
//' @param space_epsilon the truncation towards zero of the space (wealth) variable
//' @param time_res the time resolution for the time-grid
//' @param space_res the space resolution for the space-grid
//' 
//' @description {Numerically compute the optimal control function/value of the HJB PDE for optimal log utility.}
//' @return matrix
//' @export explicit_control
// [[Rcpp::export]]
double explicit_control(double bankroll, double time_length, double mu, double rate, double volat, double space_epsilon = 0.001, int time_res = 1000, int space_res = 400)
{
  double control = 0.0;
  std::vector<double> v = explicit_integrator1(bankroll,time_length, mu, rate, volat, space_epsilon, time_res, space_res);
  // Return the first derivative
  double vx = v[1];
  // Second derivative
  double vxx = v[2];
  double x = v[3];
  control = -((mu-rate)/std::pow(volat, 2.0))*(vx/(x*vxx));
  return control;
}  

//' The numerical optimal control value
//' 
//' @param bankroll the initial wealth level
//' @param time_length the time-horizon of the investment
//' @param mu the mean rate of return of the stock
//' @param rate the risk-free rate of the money-market account
//' @param volat the annual volatility of the stock
//' @param space_epsilon the truncation towards zero of the space (wealth) variable
//' @param time_res the time resolution for the time-grid
//' @param space_res the space resolution for the space-grid
//' 
//' @description {Numerically compute the optimal control function/value of the HJB PDE for optimal log utility.}
//' @return matrix
//' @export imex_control
// [[Rcpp::export]]
double imex_control(double bankroll, double time_length, double mu, double rate, double volat, double space_epsilon = 0.001, int time_res = 1000, int space_res = 400)
{
  double control = 0.0;
  std::vector<double> v = imex_integrator1(bankroll, time_length, mu, rate, volat, space_epsilon, time_res, space_res);
  // Return the first derivative
  double vx = v[1];
  // Second derivative
  double vxx = v[2];
  double x = v[3];
  control = -((mu-rate)/std::pow(volat, 2.0))*(vx/(x*vxx));
  return control;
}
  
  
