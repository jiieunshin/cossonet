#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>

void R_init_markovchain(DllInfo *dll) {
  // Register routines
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  // Enable dynamic loading
  R_useDynamicSymbols(dll, TRUE);
}

// Define the scale function
void scale(double *arr, int n) {
  double mean = 0.0, std_dev = 0.0;

  // Calculate mean
  for (int i = 0; i < n; ++i) {
    mean += arr[i];
  }
  mean /= n;

  // Calculate standard deviation
  for (int i = 0; i < n; ++i) {
    std_dev += pow(arr[i] - mean, 2);
  }
  std_dev = sqrt(std_dev / n);

  // Apply scaling
  for (int i = 0; i < n; ++i) {
    arr[i] = (arr[i] - mean) / std_dev;
  }
}

// Define the sspline_cd function
SEXP Csspline(SEXP zw, SEXP Rw, SEXP cw, SEXP sw, SEXP b, SEXP lambda0) {
  int n = length(zw);
  SEXP result = PROTECT(allocVector(VECSXP, 5)); // Extra space for b_new

  // Convert R vectors to C arrays
  double *zw_c = REAL(zw);
  double *Rw_c = REAL(Rw);
  double *cw_c = REAL(cw);
  double *sw_c = REAL(sw);
  double b_c = REAL(b)[0];
  double lambda0_c = REAL(lambda0)[0];

  // Define variables
  double *cw_new = (double *)malloc(n * sizeof(double));
  double *c_new = (double *)malloc(n * sizeof(double));

  // Main loop
  for (int iter = 0; iter < 10; ++iter) {

    // update cw
    for (int j = 0; j < n; ++j) { // iterate by column

      double V1 = 0.0, V3 = 0.0;
      for (int k = 0; k < n; ++k) { // iterate by row
        double Rc1 = 0.0;
        for (int l = 0; l < n; ++l) {
          if (l != j) {
            Rc1 += Rw_c[l * n + k] * cw_c[l];
          }
        }
        V1 += (zw_c[k] - Rc1 - b_c * sw_c[k]) * Rw_c[j * n + k];
        V3 += pow(Rw_c[j * n + k], 2);
      }

      double V2 = 0.0, Rc2 = 0.0;
      for (int l = 0; l < n; ++l) {
        if (l != j) {
          Rc2 += Rw_c[j * n + l] * cw_c[l];
        }
      }
      V2 += n * lambda0_c * Rc2;
      double V4 = Rw_c[j * n + j];

      cw_new[j] = (V1 - V2) / (V3 + V4);
    }

    // If convergence criteria are met, break the loop
    double max_diff = 1e-5;
    for (int j = 0; j < n; ++j) {
      max_diff = fmax(max_diff, fabs(cw_c[j] - cw_new[j]));
    }
    if (max_diff <= 1e-3) {
      break;
    }

    // Scale cw_new
    scale(cw_new, n);

    // Update cw_c with cw_new values
    for (int j = 0; j < n; ++j) {
      cw_c[j] = cw_new[j];
    }

  } // Update iteration

  // Apply scaling to cw_new
  scale(cw_new, n);

  // Calculate c_new
  for (int i = 0; i < n; ++i) {
    c_new[i] = cw_new[i] * sqrt(sw_c[i]);
  }

  // Calculate b_new
  double sum3 = 0.0, sum4 = 0.0;
  for (int k = 0; k < n; ++k) { // iterate by row
    double Rc = 0.0;
    for (int l = 0; l < n; ++l) { // iterate by col
      Rc += Rw_c[l * n + k] * cw_new[l];
    }
    sum3 += (zw_c[k] - Rc) * sw_c[k];
    sum4 += sw_c[k];
  }

  double b_new = sum3 / sum4;

  // Set values in result SEXP
  SET_VECTOR_ELT(result, 0, allocVector(REALSXP, n));
  SET_VECTOR_ELT(result, 1, ScalarReal(b_new));
  SET_VECTOR_ELT(result, 2, allocVector(REALSXP, n));
  SET_VECTOR_ELT(result, 3, zw);
  SET_VECTOR_ELT(result, 4, sw);

  // Copy values to result SEXP
  for (int i = 0; i < n; ++i) {
    REAL(VECTOR_ELT(result, 0))[i] = cw_new[i];
    REAL(VECTOR_ELT(result, 2))[i] = c_new[i];
  }

  // Set names for the list elements
  SEXP name_ssp = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(name_ssp, 0, mkChar("cw.new"));
  SET_STRING_ELT(name_ssp, 1, mkChar("b.new"));
  SET_STRING_ELT(name_ssp, 2, mkChar("c.new"));
  SET_STRING_ELT(name_ssp, 3, mkChar("zw.new"));
  SET_STRING_ELT(name_ssp, 4, mkChar("sw.new"));
  setAttrib(result, R_NamesSymbol, name_ssp);

  // Free dynamically allocated memory
  free(cw_new);
  free(c_new);

  UNPROTECT(1); // Unprotect result
  return result;
}

SEXP Cnng(SEXP Gw, SEXP uw, SEXP theta, SEXP lambda_theta, SEXP gamma) {
  // Convert R vectors to C arrays
  int d = ncols(Gw);
  int n = nrows(Gw);
  double *uw_c = REAL(uw);
  double *Gw_c = REAL(Gw);
  double *theta_c = REAL(theta);
  double lambda_theta_c = REAL(lambda_theta)[0];
  double gamma_c = REAL(gamma)[0];
  double r = lambda_theta_c * gamma_c / 2;
  // SEXP out = PROTECT(allocVector(VECSXP, 3));

  // Define variables
  double *theta_new = (double *)malloc(d * sizeof(double));
  double *pow_j = (double *)malloc(d * sizeof(double));

  for(int i = 0; i < 10; i++) {

    for(int j = 0; j < d; j++) { // iterate by column
      double V1 = 0.0;
      double V2 = 0.0;
      for(int k = 0; k < n; k++) { // iterate by row
        double GT = 0.0;
        for(int l = 0; l < d; l++) { // iterate by column except j
          if(l != j) {
            GT += Gw_c[l * n + k] * theta_c[l];
          }
        }
        V1 += (uw_c[k] - GT) * Gw_c[j * n + k];
        V2 += pow(Gw_c[j * n + k], 2);
      }
      theta_new[j] += V1;
      pow_j[j] += V2;
    }

    for(int j = 0; j < d; j++) {
      if(theta_new[j] > 0 && r < fabs(theta_new[j])) {
        theta_new[j] = theta_new[j] / (pow_j[j] + lambda_theta_c * (1 - gamma_c));
      } else {
        theta_new[j] = 0;
      }
    }

    // If convergence criteria are met, break the loop
    double max_diff = 1e-5;
    for (int j = 0; j < d; ++j) {
      max_diff = fmax(max_diff, fabs(theta_c[j] - theta_new[j]));
    }
    if (max_diff <= 1e-5) {
      break;
    }

    // Update cw_c with cw_new values
    for (int j = 0; j < d; ++j) {
      theta_c[j] = theta_new[j];
    }

  } // end iteration


  // Convert theta_new to an R numeric vector
  SEXP theta_new_r = PROTECT(allocVector(REALSXP, d));
  for(int i = 0; i < d; ++i) {
    REAL(theta_new_r)[i] = theta_new[i];
  }

  free(theta_new);
  free(pow_j);

  UNPROTECT(1);

  return theta_new_r;
}
