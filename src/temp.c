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

// Define the sspline_cd function
SEXP Csspline(SEXP zw, SEXP Rw, SEXP cw, SEXP sw, SEXP lambda0) {
  int n = length(zw);
  SEXP result = PROTECT(allocVector(VECSXP, 5)); // Extra space for b_new

  // Convert R vectors to C arrays
  double *zw_c = REAL(zw);
  double *Rw_c = REAL(Rw);
  double *cw_c = REAL(cw);
  double *sw_c = REAL(sw);
  double lambda0_c = REAL(lambda0)[0];

  // Define variables
  double b_c = 0;
  double *cw_new = (double *)malloc(n * sizeof(double));
  double *c_new = (double *)malloc(n * sizeof(double));
  double *pow_Rc = (double *)malloc(n * sizeof(double));

  double max_diff = 1e-6;
  int iter = 0;
  // calculate square term
  for(int j = 0; j < n; j++) { // iterate by column
    double add = 0.0;
    for(int k = 0; k < n; k++) { // iterate by row
      add += pow(Rw_c[j * n + k], 2);
    }
    pow_Rc[j] += 2 * add;
  }

  // outer loop
  for (iter = 0; iter < 80; ++iter) {

    // update cw
    for (int j = 0; j < n; ++j) { // iterate by column

      double V1 = 0.0;
      for (int k = 0; k < n; ++k) { // iterate by row
        double Rc1 = 0.0;
        for (int l = 0; l < n; ++l) {
          if (l != j) {
            Rc1 += Rw_c[l * n + k] * cw_c[l];
          }
        }
        V1 += (zw_c[k] - Rc1 - b_c * sw_c[k]) * Rw_c[j * n + k];
      }
      V1 = 2 * V1;

      double V2 = 0.0;
      for (int l = 0; l < n; ++l) {
        if (l != j) {
          V2 += Rw_c[l * n + j] * cw_c[l];
        }
      }
      V2 = n * lambda0_c * V2;

      double V4 = n * lambda0_c * Rw_c[j * n + j];

      cw_new[j] = (V1 - V2) / (pow_Rc[j] + V4);


      // If convergence criteria are met, break the loop
      for (int j = 0; j < n; ++j) {
        max_diff = fmax(max_diff, fabs(cw_c[j] - cw_new[j]));
      }
      if (max_diff <= 1e-6) {
        break;
      }

      // If not convergence, update cw
      cw_c[j] = cw_new[j];
    }

  } // end outer iteration

  if (max_diff > 1e-6 && iter == 1){
    memcpy(cw_new, cw_c, n * sizeof(double));
  }

  // Calculate c_new
  for (int i = 0; i < n; ++i) {
    c_new[i] = cw_new[i] * sqrt(sw_c[i]);
  }

  // result
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

  UNPROTECT(2); // Unprotect result
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
  double r = lambda_theta_c * gamma_c * n;
  // SEXP out = PROTECT(allocVector(VECSXP, 3));

  // Define variables
  double *theta_new = (double *)malloc(d * sizeof(double));
  double *pow_theta = (double *)malloc(d * sizeof(double));

  double max_diff = 1e-6;
  int iter = 0;

  for(int j = 0; j < d; j++) { // iterate by column
    double V2 = 0.0;
    for(int k = 0; k < n; k++) { // iterate by row
      V2 += pow(Gw_c[j * n + k], 2);
    }
    pow_theta[j] += V2;
  }

  // outer iteration
  for(iter = 0; iter < 80; iter++) {

    for(int j = 0; j < d; j++) { // iterate by column
      double V1 = 0.0;
      for(int k = 0; k < n; k++) { // iterate by row
        double GT = 0.0;
        for(int l = 0; l < d; l++) { // iterate by column except j
          if(l != j) {
            GT += Gw_c[l * n + k] * theta_c[l];
          }
        }
        V1 += (uw_c[k] - GT) * Gw_c[j * n + k];
      }
      theta_new[j] += 2 * V1;

      for(int j = 0; j < d; j++) {
        if(theta_new[j] > 0 && r < fabs(theta_new[j])) {
          theta_new[j] = (theta_new[j] - r) / (pow_theta[j] + n * lambda_theta_c * (1 - gamma_c));
        } else {
          theta_new[j] = 0;
        }
      }

      // If convergence criteria are met, break the loop
     for (int j = 0; j < d; ++j) {
        max_diff = fmax(max_diff, fabs(theta_c[j] - theta_new[j]));
      }

      if (max_diff <= 1e-6) {
        break;
      }

      // Update theta_c with cw_new values
      for (int j = 0; j < d; ++j) {
        theta_c[j] = theta_new[j];
      }

    }
  } // end outer iteration

  if (max_diff > 1e-6 && iter == 1){
    theta_new = (double *)malloc(d * sizeof(double));
  }

  // result
  SEXP theta_new_r = PROTECT(allocVector(REALSXP, d));
  for(int j = 0; j < d; ++j) {
    REAL(theta_new_r)[j] = theta_new[j];
  }

  free(theta_new);
  free(pow_theta);

  UNPROTECT(1);

  return theta_new_r;
}
