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
  double *L = (double *)malloc(n * sizeof(double));
  double *R = (double *)malloc(n * sizeof(double));

  // Main loop
  for (int iter = 0; iter < 20; ++iter) {

    // update cw
    for (int j = 0; j < n; ++j) { // iterate by column
      double L = 0.0, R = 0.0;
      for (int k = 0; k < n; ++k) { // iterate by row
        double Rc = 0.0;
        for (int l = 0; l < n; ++l) {
          if (l != j) {
            Rc += Rw_c[l * n + k] * cw_c[l];
          }
        }
        L += 2 * (zw_c[k] - Rc - b_c * sw_c[k]) * Rw_c[j * n + k] - lambda0_c * Rc;
        R += 2 * Rw_c[j * n + k] * Rw_c[j * n + k] + lambda0_c * cw_c[j] * Rw_c[j * n + j];
      }
      cw_new[j] = L / R;
    }

    // Scale cw_new
    scale(cw_new, n);
    double max_diff = 0.0;

    for (int j = 0; j < n; ++j) {
      max_diff = fmax(max_diff, fabs(cw_c[j] - cw_new[j]));
    }

    // If convergence criteria are met, break the loop
    if (max_diff < 1e-5) {
      break;
    }

    // Update cw_c with cw_new values
    for (int j = 0; j < n; ++j) {
      cw_c[j] = cw_new[j];
    }

    // Update iteration
  }

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
  SEXP names = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(names, 0, mkChar("cw.new"));
  SET_STRING_ELT(names, 1, mkChar("b.new"));
  SET_STRING_ELT(names, 2, mkChar("c.new"));
  SET_STRING_ELT(names, 3, mkChar("zw.new"));
  SET_STRING_ELT(names, 4, mkChar("sw.new"));
  setAttrib(result, R_NamesSymbol, names);

  // Free dynamically allocated memory
  free(cw_new);
  free(L);
  free(R);
  free(c_new);

  UNPROTECT(1); // Unprotect result
  return result;
}

SEXP Cnng(SEXP Gw, SEXP uw, SEXP theta, SEXP lambda_theta, SEXP gamma) {
  int d = ncols(Gw);
  double r = REAL(lambda_theta)[0] * REAL(gamma)[0] / 2;
  SEXP theta_new = PROTECT(allocVector(REALSXP, d));
  SEXP out;
  int i, j, k;
  double sum_value;

  for(i = 0; i < 10; i++) {
    for(j = 0; j < d; j++) {
      sum_value = 0;
      for(k = 0; k < d; k++) {
        if(k != j) {
          sum_value += (REAL(uw)[k] - REAL(Gw)[j + k * d] * REAL(theta)[k]) * REAL(Gw)[j + k * d];
        }
      }
      REAL(theta_new)[j] = sum_value;
    }

    for(j = 0; j < d; j++) {
      if(REAL(theta_new)[j] > 0 && r < fabs(REAL(theta_new)[j])) {
        REAL(theta_new)[j] = REAL(theta_new)[j];
      } else {
        REAL(theta_new)[j] = 0;
      }
      double sum_sq = 0;
      for(k = 0; k < d; k++) {
        sum_sq += pow(REAL(Gw)[j + k * d], 2);
      }
      REAL(theta_new)[j] = REAL(theta_new)[j] / (sum_sq + REAL(lambda_theta)[0] * (1 - REAL(gamma)[0]));
    }

    double max_diff = 0.0;
    for (int j = 0; j < d; ++j) {
      max_diff = fmax(max_diff, fabs(REAL(theta)[j] - REAL(theta_new)[j]));
    }

    // If convergence criteria are met, break the loop
    if (max_diff < 1e-5) {
      break;
    }

  }

  PROTECT(out = allocVector(VECSXP, 3));
  SET_VECTOR_ELT(out, 0, lambda_theta);
  SET_VECTOR_ELT(out, 1, gamma);
  SET_VECTOR_ELT(out, 2, theta_new);

  // Set names for the list elements
  SEXP names = PROTECT(allocVector(STRSXP, 5));
  SET_STRING_ELT(names, 0, mkChar("lambda_theta"));
  SET_STRING_ELT(names, 1, mkChar("gamma"));
  SET_STRING_ELT(names, 2, mkChar("theta.new"));
  setAttrib(out, R_NamesSymbol, names);

  UNPROTECT(2);

  return out;
}
