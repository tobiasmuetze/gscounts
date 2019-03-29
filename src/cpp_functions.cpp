#include <Rcpp.h>
#include <math.h>
#include <string>
using namespace Rcpp;

double f_k(double zk, double Ik, double zk1, double Ik1, double theta) {
  // implements equation (19.4) from Jennison-Turbull p.347 
  double delta = Ik - Ik1;
  double x = (zk * sqrt(Ik) - zk1 * sqrt(Ik1) - theta * delta) / sqrt(delta);
  double pi = 3.14159265359;
  double out;
  out = sqrt(Ik / delta) * exp(-x * x / 2) / sqrt(2*pi);
  return(out);
}


double cppnorm(double x) {
  return 0.5 * erfc(-x * M_SQRT1_2);
}

// side="lower": solves for c: P(lower[1] < X1 < upper[1], ..., -Inf < X < c) = error_spend
// side="upper": solves for c: P(lower[1] < X1 < upper[1], ...,  c < X < Inf) = error_spend
// [[Rcpp::export]]
double cpp_calc_critical(int r, NumericVector lower, NumericVector upper, double error_spend, NumericVector information, double theta, String side) {
  
  int k = lower.size();
  int m = 12 * r - 3;
  int j, l;
  double integral = 0;
  double shift;
  NumericVector grid_base(m);
  NumericVector weight_base(m);
  NumericMatrix grid(k, m);
  NumericMatrix weight(k, m);
  NumericMatrix h(k, m);
  // Variables for determing the limit
  double x;
  double delta = information(k) - information(k-1);
  double sqrt_delta = sqrt(delta);
  double sqrt_info_k = sqrt(information(k));
  double sqrt_info_k1 = sqrt(information(k-1));
  double crit;
  double f = 1, f_grad = 0, eps = 0.000001;
  double f_left, f_right; // variables for starting value search
  int n_iterate = 0;
  
  if (error_spend == 0 && side == "lower")
    return(R_NegInf);
  if (error_spend == 0 && side == "upper")
    return(R_PosInf);
  
  // base grid points according to Jennison-Turbull (JB) p.349f.
  grid_base[0] = -3 - 4 * log((double) r);
  for(int i = 1; i < 6*r-1; i++) {
    if (i <= r-2) {
      grid_base[2*i] = -3 - 4 * log((double) r / (i+1));
    }
    if ((i >= r-1) && (i <= 5*r-1)) {
      grid_base[2*i] = -3 + 3 * (double) (i + 1 - r) / (2*r);
    }
    if (i >= 5*r) {
      grid_base[2*i] = 3 + 4 * log((double) r / (6*r-i-1));
    }
    grid_base[2*i-1] = (grid_base[2*i] + grid_base[2*i-2])/2;
  }
  
  // weights to grid points according to Jennison-Turbull (JB) p.350.
  // weights to grid points according to Jennison-Turbull (JB) p.350.
  weight_base[0] = (grid_base[2] - grid_base[0])/6;
  weight_base[m-2] = 4*(grid_base[m-1] - grid_base[m-3])/6;
  weight_base[m-1] = (grid_base[m-1] - grid_base[m-3])/6;
  for(int i = 1; i < 6*r-2; i++) {
    weight_base[2*i] = (grid_base[2*i+2] - grid_base[2*i-2])/6;
    weight_base[2*i-1] = 4*(grid_base[2*i] - grid_base[2*i-2])/6;
  }

  // Grid and weight matrices including the shift
  for(int i = 0; i < k; i++) {
    for(int j = 0; j < m; j++) {
      grid(i, j) = theta * sqrt(information(i)) + grid_base[j];
      weight(i, j) = weight_base[j];
    }
  }
  
  // adjustments if integration bounds are within base grid points; confer Jennison-Turbull p. 350
  for(int i = 0; i < k; i++) {
    shift = theta * sqrt(information(i));
    // adjustments to be made if lower boundary is not smaller than smallest grid point
    if (lower(i) > shift + grid_base[0]) {
      for(int j = 0; j < 6*r-3; j++) { // upper bound for j to avoid going out of bounds
        if (lower(i) > shift + grid_base[2*j+2]) {
          grid(i, 2*j) = grid(i, 2*j+1) = 0;
          weight(i, 2*j) = weight(i, 2*j+1) = 0;
        } else {
          grid(i, 2*j) = lower(i);
          grid(i, 2*j+1) = (lower(i) + grid(i, 2*j+2))/2;
          weight(i, 2*j) = (grid(i, 2*j+2) - grid(i, 2*j))/6;
          weight(i, 2*j+1) = 4*(grid(i, 2*j+2) - grid(i, 2*j))/6;
          weight(i, 2*j+2) = (grid(i, 2*j+4) - grid(i, 2*j))/6;
          break;
        }
      }
    }
    // adjustments to be made if upper boundary is not larger than largest grid point
    if (upper(i) < shift + grid_base[m-1]) {
      for(int j = 6*r-2; j > 2; j--) {
        if (upper(i) <= shift + grid_base[2*j-2]) {
          grid(i, 2*j) = grid(i, 2*j-1) = 0;
          weight(i, 2*j) = weight(i, 2*j-1) = 0;
        } else {
          grid(i, 2*j) = upper(i);
          grid(i, 2*j-1) = (upper(i) + grid(i, 2*j-2))/2;
          weight(i, 2*j) = (grid(i, 2*j) - grid(i, 2*j-2))/6;
          weight(i, 2*j-1) = 4*(grid(i, 2*j) - grid(i, 2*j-2))/6;
          weight(i, 2*j-2) = (grid(i, 2*j) - grid(i, 2*j-4))/6;
          break;
        }
      }
    }
  }
  
  
  // Calculate h. Every row of matrix corresponds to one integral (J&T, p.348)
  for(int i = 0; i < m; i++) {
    if (weight(0, i) != 0)
      h(0, i) = weight(0, i) * f_k(grid(0, i), information(0), 0, 0, theta);
  }
  for(int i = 1; i < k; i++) { // number of integrals
    for(int j = 0; j < m; j++) {
      if (weight(i, j) != 0) {
        for(int l = 0; l < m; l++) {
          if (h(i-1, l) != 0) {
            h(i, j) += h(i-1, l) * f_k(grid(i, j), information(i), grid(i-1, l), information(i-1), theta);
          }
        }
        h(i, j) = h(i, j) *  weight(i, j);
      } //else {
      //h(i, j) = 0;
      //}
    }
  }
  
  // Newton-Rapshon algorithm for calculating the limit of interval (-Inf, crit) (side = "lower") or (crit, Inf) (side = "upper") (J&T, p.353f) 
  // starting value search for Newton-Rapshon
  crit = -5 + theta * sqrt_info_k;
  f_left = 0;
  while (crit <= 5 + theta * sqrt_info_k) {
    // calculate f
    f_right = -error_spend;
    for(int i = 0; i < m; i++) {
      if (h(k-1, i) != 0) {
        x = (grid(k-1, i) * sqrt_info_k1 + theta * delta - crit * sqrt_info_k) / sqrt_delta;
        if (side == "lower")
          f_right += h(k-1, i) * (1 - cppnorm(x));
        if (side == "upper")
          f_right += h(k-1, i) * cppnorm(x);
      }
    }
    // stop stat value search if root is considered interval
    if (((f_left < 0) && (f_right > 0)) || ((f_left > 0) && (f_right < 0))) {
      crit = crit - 0.5;
      break;
    }
    else {
      f_left = f_right;
      crit += 1;
    }
  }
  
  // No solution exists if no starting value was found
  if (crit >= 5 + theta * sqrt_info_k) {
    if (side == "lower")
      return(99999);      
    if (side == "upper")
      return(-99999);      
  }
  
  while (fabs(f) > eps) {
    if (n_iterate > 20) {
      Rf_warning("Maximum number of iterations reached");
      break;
    }
    
    // calculate f = function which root must be found
    f = -error_spend;
    f_grad = 0;
    // calculate f and derivation of f
    if (side == "lower") {
      for(int i = 0; i < m; i++) {
        if (h(k-1, i) != 0) {
          x = (grid(k-1, i) * sqrt_info_k1 + theta * delta - crit * sqrt_info_k) / sqrt_delta;
          f += h(k-1, i) * (1 - cppnorm(x));
          f_grad += h(k-1, i) * f_k(crit, information(k), grid(k-1, i), information(k-1), theta);
          continue;
        }
      }
    }
    if (side == "upper") {
      for(int i = 0; i < m; i++) {
        if (h(k-1, i) != 0) {
          x = (grid(k-1, i) * sqrt_info_k1 + theta * delta - crit * sqrt_info_k) / sqrt_delta;
          f += h(k-1, i) * cppnorm(x);
          f_grad -= h(k-1, i) * f_k(crit, information(k), grid(k-1, i), information(k-1), theta);
          continue;
        }
      }
    }
    
    
    // Newton step
    crit = crit - f / f_grad;
    n_iterate++;
  }
  
  
  
  return(crit);
}




// [[Rcpp::export]]
double cpp_pmultinorm(int r, NumericVector lower, NumericVector upper, NumericVector information, double theta) {
  
  int k = lower.size();
  int m = 12 * r - 3;
  int j, l;
  double integral = 0;
  double shift;
  double delta, x_lower, x_upper, sqrt_info_k1, sqrt_info_k2, sqrt_delta;
  NumericVector grid_base(m);
  NumericVector weight_base(m);
  NumericMatrix grid(k-1, m);
  NumericMatrix weight(k-1, m);
  NumericMatrix h(k-1, m);
  
  // grid points according to Jennison-Turbull (JB) p.349f.
  grid_base[0] = -3 - 4 * log((double) r);
  for(int i = 1; i < 6*r-1; i++) {
    if (i <= r-2) {
      grid_base[2*i] = -3 - 4 * log((double) r / (i+1));
    }
    if ((i >= r-1) && (i <= 5*r-1)) {
      grid_base[2*i] = -3 + 3 * (double) (i + 1 - r) / (2*r);
    }
    if (i >= 5*r) {
      grid_base[2*i] = 3 + 4 * log((double) r / (6*r-i-1));
    }
    grid_base[2*i-1] = (grid_base[2*i] + grid_base[2*i-2])/2;
  }
  
  // weights to grid points according to Jennison-Turbull (JB) p.350.
  weight_base[0] = (grid_base[2] - grid_base[0])/6;
  weight_base[m-2] = 4*(grid_base[m-1] - grid_base[m-3])/6;
  weight_base[m-1] = (grid_base[m-1] - grid_base[m-3])/6;
  for(int i = 1; i < 6*r-2; i++) {
    weight_base[2*i] = (grid_base[2*i+2] - grid_base[2*i-2])/6;
    weight_base[2*i-1] = 4*(grid_base[2*i] - grid_base[2*i-2])/6;
  }
  
  // write grid and weight matrices
  for(int i = 0; i < k-1; i++) {
    for(int j = 0; j < m; j++) {
      grid(i, j) = theta * sqrt(information(i)) + grid_base[j];
      weight(i, j) = weight_base[j];
    }
  }
  
  // adjustments if integration bounds are within base grid points
  // JB page 350 mid-page
  for(int i = 0; i < k-1; i++) {
    shift = theta * sqrt(information(i));
    // printf("shift %f %f %f \n", shift, lower(i) , shift + grid_base[0]);
    // printf("shift %f\n", shift);
    // adjustments to be made if lower boundary is not smaller than smallest grid point
    if (lower(i) > shift + grid_base[0]) {
      for(int j = 0; j < 6*r-3; j++) { // upper bound for j to avoid going out of bounds
        if (lower(i) > shift + grid_base[2*j+2]) {
          grid(i, 2*j) = grid(i, 2*j+1) = 0;
          weight(i, 2*j) = weight(i, 2*j+1) = 0;
        } else {
          grid(i, 2*j) = lower(i);
          grid(i, 2*j+1) = (lower(i) + grid(i, 2*j+2))/2;
          weight(i, 2*j) = (grid(i, 2*j+2) - grid(i, 2*j))/6;
          weight(i, 2*j+1) = 4*(grid(i, 2*j+2) - grid(i, 2*j))/6;
          weight(i, 2*j+2) = (grid(i, 2*j+4) - grid(i, 2*j))/6;
          break;
        }
      }
    }
    // adjustments to be made if upper boundary is not larger than largest grid point
    if (upper(i) < shift + grid_base[m-1]) {
      for(int j = 6*r-2; j > 1; j--) {
        if (upper(i) <= shift + grid_base[2*j-2]) {
          grid(i, 2*j) = grid(i, 2*j-1) = 0;
          weight(i, 2*j) = weight(i, 2*j-1) = 0;
        } else {
          grid(i, 2*j) = upper(i);
          grid(i, 2*j-1) = (upper(i) + grid(i, 2*j-2))/2;
          weight(i, 2*j) = (grid(i, 2*j) - grid(i, 2*j-2))/6;
          weight(i, 2*j-1) = 4*(grid(i, 2*j) - grid(i, 2*j-2))/6;
          weight(i, 2*j-2) = (grid(i, 2*j) - grid(i, 2*j-4))/6;
          break;
        }
      }
    }
  }
  
  for(int i = 0; i < m; i++) {
    if (weight(0, i) != 0)
      h(0, i) = weight(0, i) * f_k(grid(0, i), information(0), 0, 0, theta);
  }
  
  for(int i = 1; i < k-1; i++) { // number of integrals
    for(int j = 0; j < m; j++) {
      if (weight(i, j) != 0) {
        for(int l = 0; l < m; l++) {
          if (h(i-1, l) != 0) {
            h(i, j) += h(i-1, l) * f_k(grid(i, j), information(i), grid(i-1, l), information(i-1), theta);
          }
        }
        h(i, j) = h(i, j) *  weight(i, j);
      } else {
        h(i, j) = 0;
      }
    }
  }
  
  // Calculate integral using Formula (19.8) J&T
  delta = information(k-1) - information(k-2);
  sqrt_delta = sqrt(delta);
  sqrt_info_k1 = sqrt(information(k-1));
  sqrt_info_k2 = sqrt(information(k-2));
  
  for(int i = 0; i < m; i++) {
    if (h(k-2, i) != 0) {
      if (lower(k-1) == -INFINITY) {
        x_upper = (grid(k-2, i) * sqrt_info_k2 + theta * delta - upper(k-1) * sqrt_info_k1) / sqrt_delta;
        integral += h(k-2, i) * (1 - cppnorm(x_upper));
        continue;
      }
      if (upper(k-1) == INFINITY) {
        x_lower = (grid(k-2, i) * sqrt_info_k2 + theta * delta - lower(k-1) * sqrt_info_k1) / sqrt_delta;
        integral += h(k-2, i) * cppnorm(x_lower);
        continue;
      }
      x_lower = (grid(k-2, i) * sqrt_info_k2 + theta * delta - lower(k-1) * sqrt_info_k1) / sqrt_delta;
      x_upper = (grid(k-2, i) * sqrt_info_k2 + theta * delta - upper(k-1) * sqrt_info_k1) / sqrt_delta;
      integral += h(k-2, i) * (cppnorm(x_lower) - cppnorm(x_upper));
    }
  }
  return(integral);
}

