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


// Define the sspline_cd fulention
SEXP glm_c_step(SEXP zw, SEXP Rw, SEXP Rw2, SEXP cw, SEXP sw, SEXP tr_n, SEXP N, SEXP lambda0) {
  int tr_nc = INTEGER(tr_n)[0];
  int Nc = INTEGER(N)[0];

  SEXP result = PROTECT(allocVector(VECSXP, 2)); // Extra space for b_new

  // Convert R vectors to C arrays
  double *zw_c = REAL(zw);
  double *Rw_c = REAL(Rw);
  double *Rw2_c = REAL(Rw2);
  double *cw_c = REAL(cw);
  double *sw_c = REAL(sw);
  double lambda0_c = REAL(lambda0)[0];

  // Define variables
  // double *cw_new = (double *)malloc(nc * sizeof(double));
  // double *c_new = (double *)malloc(nc * sizeof(double));
  double *pow_Rc = (double *)malloc(Nc * sizeof(double));
  double cw_new;

  if (pow_Rc == NULL) {
    error("Memory allocation failed");
  }

  for(int k = 0; k < Nc; k++) {
    // cw_new[k] = 0;
    // c_new[k] = 0;
    pow_Rc[k] = 0;
  }

  // calculate square term
  for(int j = 0; j < Nc; ++j) { // iterate by column
    double add = 0.0;
    for(int k = 0; k < tr_nc; ++k) { // iterate by row
      add += Rw_c[j * tr_nc + k] * Rw_c[j * tr_nc + k];
    }
    pow_Rc[j] = add;
  }

  int iter = 0;
  // double min_diff = 10;
  double diff = 0.0;
  double avg_diff;
  double sum3;
  double sum4;
  double b_new = 0.0;

  // outer loop
  for (iter = 0; iter < 40; ++iter) {
    avg_diff = 0.0;

    // update cw
    for (int j = 0; j < Nc; ++j) { // iterate by column

      double V1 = 0.0;
      for (int k = 0; k < tr_nc; ++k) { // iterate by row
        double Rc1 = 0.0;
        for (int l = 0; l < Nc; ++l) {
          if (l != j) {
            Rc1 += Rw_c[l * tr_nc + k] * cw_c[l];
          }
        }
        V1 += (zw_c[k] - Rc1 - b_new * sw_c[k]) * Rw_c[j * tr_nc + k];
      }
      // V1 = 2 * V1;

      double V2 = 0.0;
      for (int l = 0; l < Nc; ++l) {
        if (l != j) {
          V2 += Rw2_c[l * Nc + j] * cw_c[l];
        }
      }
      V2 = lambda0_c * V2;

      double V4 = lambda0_c * Rw2_c[j * Nc + j];

      cw_new = (V1 - V2) / (pow_Rc[j] + V4);

      // If convergence criteria are met, break the loop
      diff = fabs(cw_c[j] - cw_new) / fabs(cw_c[j] + 1);

      if ((diff < 1e-6) | (diff > 10)) {
        cw_new = cw_c[j];
        diff = fabs(cw_c[j] - cw_new) / fabs(cw_c[j] + 1);
      }

      avg_diff += diff;

      // Update cw with the new value
      cw_c[j] = cw_new;
    }

    // result
    sum3 = 0.0;
    sum4 = 0.0;
    for (int k = 0; k < tr_nc; ++k) { // iterate by row
      double Rc = 0.0;
      for (int l = 0; l < Nc; ++l) { // iterate by col
        Rc += Rw_c[l * tr_nc + k] * cw_c[l];
      }
      sum3 += (zw_c[k] - Rc) * sw_c[k];
      sum4 += sw_c[k] * sw_c[k];
    }

    b_new = sum3 / sum4;


    avg_diff /= Nc; // Calculate the average difference

    // Check for convergence based on average difference
    if (avg_diff <= 1e-4) {
      break;
    }
  } // End outer iteration

  Rprintf("min_diff: %f\n", avg_diff);
  Rprintf("iter: %d\n", iter);


  // Set values in result SEXP
  SET_VECTOR_ELT(result, 0, allocVector(REALSXP, Nc));
  SET_VECTOR_ELT(result, 1, ScalarReal(b_new));
  // SET_VECTOR_ELT(result, 2, allocVector(REALSXP, nc));
  // SET_VECTOR_ELT(result, 2, zw);
  // SET_VECTOR_ELT(result, 3, sw);

  // Copy values to result SEXP
  for (int i = 0; i < Nc; ++i) {
    REAL(VECTOR_ELT(result, 0))[i] = cw_c[i];
    // REAL(VECTOR_ELT(result, 2))[i] = c_new[i];
  }

  // Set names for the list elements
  SEXP name_ssp = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(name_ssp, 0, mkChar("cw.new"));
  SET_STRING_ELT(name_ssp, 1, mkChar("b.new"));
  // SET_STRING_ELT(name_ssp, 2, mkChar("c.new"));
  // SET_STRING_ELT(name_ssp, 2, mkChar("zw.new"));
  // SET_STRING_ELT(name_ssp, 3, mkChar("sw.new"));
  setAttrib(result, R_NamesSymbol, name_ssp);

  // Free dynamically allocated memory
  // free(cw_new);
  free(pow_Rc);
  // free(c_new);

  UNPROTECT(2); // Unprotect result
  return result;
}


SEXP glm_theta_step(SEXP Gw, SEXP uw, SEXP h, SEXP n, SEXP d, SEXP theta, SEXP r1, SEXP r2) {
  // Convert R vectors to C arrays
  int nc = INTEGER(n)[0]; // Extract the value of n
  int dc = INTEGER(d)[0]; // Extract the value of d
  double *Gw_c = REAL(Gw);
  double *uw_c = REAL(uw);
  double *h_c = REAL(h);
  double *theta_c = REAL(theta);
  double r1_c = REAL(r1)[0];
  double r2_c = REAL(r2)[0];

  // SEXP out = PROTECT(allocVector(VECSXP, 3));

  // Define variables
  // double *theta_new = (double *)malloc(dc * sizeof(double));
  double *pow_theta = (double *)malloc(dc * sizeof(double));
  double theta_new;

  // 메모리 할당 실패 시 처리 방법
  if (pow_theta == NULL) {
    error("Memory allocation failed for pow_theta");
  }

  for (int k = 0; k < dc; ++k){
    // theta_new[k] = 0;
    pow_theta[k] = 0;
  }

  for(int j = 0; j < dc; j++) { // iterate by column
    double add = 0.0;
    for(int k = 0; k < nc; k++) { // iterate by row
      add += Gw_c[j * nc + k] * Gw_c[j * nc + k];
    }
    pow_theta[j] = add;
  }

  // outer iteration
  int iter = 0;
  // double min_diff = 10.0;
  double *diff = (double *)malloc(dc * sizeof(double));
  double avg_diff;

  for(iter = 0; iter < 20; ++iter) {
    avg_diff = 0;  // Initialize avg_diff for averaging

    for(int j = 0; j < dc; ++j) { // iterate by column

      double V1 = 0.0;
      for(int k = 0; k < nc; ++k) { // iterate by row
        double GT = 0.0;
        for(int l = 0; l < dc; ++l) { // iterate by column except j
          if(l != j) {
            GT += Gw_c[l * nc + k] * theta_c[l];
          }
        }
        V1 += (uw_c[k] - GT) * Gw_c[j * nc + k];
      }
      V1 = V1 - h_c[j];
      // V1 = 2 * V1;


      if(V1 > 0 && r1_c < fabs(V1)) {
        theta_new = (V1 - r1_c) / (pow_theta[j] + r2_c);
        // theta_new = V1 / (pow_theta[j] + nc * lambda_theta_c * (1-gamma_c)) / 2;
      } else {
        theta_new = 0;
      }

      // If convergence criteria are met, break the loop
      diff[j] = fabs(theta_c[j] - theta_new) / fabs(theta_c[j] + 1);

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

  // print
  // Rprintf("\n iter: %d \n", iter);
  // for (int k = 0; k < dc; ++k){
  //   Rprintf("theta_new: %g \n", theta_c[k]);
  // }
  // Rprintf("\n max_diff: %d \n", max_diff);

  // Rprintf("min_diff: %f\n", avg_diff);
  // Rprintf("iter: %d\n", iter);

  // result
  SEXP theta_new_r = PROTECT(allocVector(REALSXP, dc));
  for(int j = 0; j < dc; ++j) {
    REAL(theta_new_r)[j] = theta_c[j];
  }

  free(pow_theta);

  UNPROTECT(1);

  return theta_new_r;
}
