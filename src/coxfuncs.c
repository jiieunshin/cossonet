#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>

// Define the sspline_cd fulention
SEXP cox_c_step(SEXP cinit, SEXP Rtheta, SEXP Rtheta2, SEXP n, SEXP m, SEXP z, SEXP w, SEXP lambda0) {
  int nc = INTEGER(n)[0];
  int mc = INTEGER(m)[0];

  SEXP result = PROTECT(allocVector(VECSXP, 2)); // Extra space for b_new

  // Convert R vectors to C arrays
  double *z_c = REAL(z);
  double *w_c = REAL(w);
  double *R1_c = REAL(Rtheta);
  double *R2_c = REAL(Rtheta2);
  double *cold_c = REAL(cinit);
  double lambda0_c = REAL(lambda0)[0];

  // Define variables
  // double *cw_new = (double *)malloc(nc * sizeof(double));
  // double *c_new = (double *)malloc(nc * sizeof(double));
  double *pow_Rc = (double *)malloc(mc * sizeof(double));
  double cw_new;

  if (pow_Rc == NULL) {
    error("Memory allocation failed");
  }

  for(int k = 0; k < mc; k++) {
    pow_Rc[k] = 0;
  }

  // calculate square term
  for(int j = 0; j < mc; ++j) { // iterate by column
    double add = 0.0;
    for(int k = 0; k < nc; ++k) { // iterate by row
      add += w_c[k] * R1_c[j * nc + k] * R1_c[j * nc + k];
    }
    pow_Rc[j] = add / nc;
  }

  int iter = 0;
  double diff = 0.0;
  double avg_diff;

  // outer loop
  for (iter = 0; iter < 40; ++iter) {
    avg_diff = 0.0;

    // update cw
    for (int j = 0; j < mc; ++j) { // iterate by column

      double V1 = 0.0;
      for (int k = 0; k < nc; ++k) { // iterate by row
        double Rc1 = 0.0;
        for (int l = 0; l < mc; ++l) {
          if (l != j) {
            Rc1 += R1_c[l * nc + k] * cold_c[l];
          }
        }
        V1 += (z_c[k] - Rc1) * R1_c[j * nc + k] * w_c[k];
      }
      V1 = V1 / nc;

      double V4 = 2 * lambda0_c * R2_c[j * mc + j];

      cw_new = V1 / (pow_Rc[j] + V4);

      // If convergence criteria are met, break the loop
      diff = fabs(cold_c[j] - cw_new) / fabs(cold_c[j] + 1);

      if ((diff < 1e-6) | (diff > 10)) {
        cw_new = cold_c[j];
        diff = fabs(cold_c[j] - cw_new);
      }

      avg_diff += diff;

      // Update cw with the new value
      cold_c[j] = cw_new;
    }

    avg_diff /= nc; // Calculate the average difference

    // Check for convergence based on average difference
    if (avg_diff <= 1e-4) {
      break;
    }
  } // End outer iteration

  // Rprintf("min_diff: %f\n", avg_diff);
  // Rprintf("iter: %d\n", iter);


  // Set values in result SEXP
  SET_VECTOR_ELT(result, 0, allocVector(REALSXP, mc));

  // Copy values to result SEXP
  for (int i = 0; i < mc; ++i) {
    REAL(VECTOR_ELT(result, 0))[i] = cold_c[i];
  }

  // Set names for the list elements
  SEXP name_ssp = PROTECT(allocVector(STRSXP, 1));
  SET_STRING_ELT(name_ssp, 0, mkChar("c.new"));
  setAttrib(result, R_NamesSymbol, name_ssp);

  // Free dynamically allocated memory
  free(pow_Rc);

  UNPROTECT(2); // Unprotect result
  return result;
}


SEXP cox_theta_step(SEXP theta_init, SEXP S, SEXP n, SEXP d, SEXP z, SEXP w, SEXP lambda_theta, SEXP gamma) {
  // Convert R vectors to C arrays
  int nc = INTEGER(n)[0]; // Extract the value of n
  int dc = INTEGER(d)[0]; // Extract the value of d
  double *S_c = REAL(S);
  double *z_c = REAL(z);
  double *w_c = REAL(w);
  double *theta_c = REAL(theta_init);
  double lambda_theta_c = REAL(lambda_theta)[0];
  double gamma_c = REAL(gamma)[0];
  double *pow_S = (double *)malloc(dc * sizeof(double));
  double theta_new;

  // 메모리 할당 실패 시 처리 방법
  if (pow_S == NULL) {
    error("Memory allocation failed for pow_theta");
  }

  for (int k = 0; k < dc; ++k){
    pow_S[k] = 0;
  }

  // calculate square term
  for(int j = 0; j < dc; ++j) { // iterate by column
    double add = 0.0;
    for(int k = 0; k < nc; ++k) { // iterate by row
      add += w_c[k] * S_c[j * nc + k] * S_c[j * nc + k];
    }
    pow_S[j] = add;
  }

  // outer iteration
  int iter = 0;
  double *diff = (double *)malloc(dc * sizeof(double));
  double avg_diff;

  for(iter = 0; iter < 40; ++iter) {
    avg_diff = 0;  // Initialize avg_diff for averaging

    for(int j = 0; j < dc; ++j) { // iterate by column

      double V1 = 0.0;
      for (int k = 0; k < nc; ++k) { // iterate by row
        double SS = 0.0;
        for (int l = 0; l < dc; ++l) {
          if (l != j) {
            SS += S_c[l * nc + k] * theta_c[l];
          }
        }
        V1 += (z_c[k] - SS) * S_c[j * nc + k] * w_c[k];
      }
      V1 = V1 / nc;

      if(V1 > 0 && (gamma_c * lambda_theta_c) < fabs(V1)) {
        theta_new = (V1 - gamma_c * lambda_theta_c) / (pow_S[j] + 2 * (1-gamma_c) * lambda_theta_c);
      } else {
        theta_new = 0;
      }

      // If convergence criteria are met, break the loop
      diff[j] = fabs(theta_c[j] - theta_new);

      avg_diff += diff[j];

      if ((theta_new > 0) & (diff[j] <= 1e-4)) {
        break;
      }

      // Update theta with the new value
      theta_c[j] = theta_new;
    }

    avg_diff /= dc;  // Calculate the average difference

    // Check for convergence based on average difference
    if (avg_diff <= 1e-4) {
      break;
    }
  } // end outer iteration

  if (iter == 0){
    for (int k = 0; k < dc; ++k){
      theta_c[k] = 0;
    }
  }

  // result
  SEXP theta_new_r = PROTECT(allocVector(REALSXP, dc));
  for(int j = 0; j < dc; ++j) {
    REAL(theta_new_r)[j] = theta_c[j];
  }

  free(pow_S);

  UNPROTECT(1);

  return theta_new_r;
}
