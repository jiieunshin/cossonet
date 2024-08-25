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
SEXP glm_c_step(SEXP zw, SEXP Rw, SEXP Rw2, SEXP cw, SEXP sw, SEXP m, SEXP n, SEXP lambda0) {
  int mc = INTEGER(m)[0];
  int nc = INTEGER(n)[0];

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
  double *pow_Rc = (double *)malloc(nc * sizeof(double));
  double cw_new;

  if (pow_Rc == NULL) {
    error("Memory allocation failed");
  }

  for(int k = 0; k < nc; k++) {
    // cw_new[k] = 0;
    // c_new[k] = 0;
    pow_Rc[k] = 0;
  }

  // calculate square term
  for(int j = 0; j < nc; ++j) { // iterate by column
    double add = 0.0;
    for(int k = 0; k < mc; ++k) { // iterate by row
      add += Rw_c[j * mc + k] * Rw_c[j * mc + k];
    }
    pow_Rc[j] = 2 * add;
  }

  int iter = 0;
  double min_diff = 100;
  double diff;

  // outer loop
  for (iter = 0; iter < 80; ++iter) {

    // update cw
    for (int j = 0; j < nc; ++j) { // iterate by column

      double V1 = 0.0;
      for (int k = 0; k < mc; ++k) { // iterate by row
        double Rc1 = 0.0;
        for (int l = 0; l < nc; ++l) {
          if (l != j) {
            Rc1 += Rw_c[l * mc + k] * cw_c[l];
          }
        }
        V1 += (zw_c[k] - Rc1) * Rw_c[j * mc + k];
      }
      V1 = 2 * V1;

      double V2 = 0.0;
      for (int l = 0; l < nc; ++l) {
        if (l != j) {
          V2 += Rw2_c[l * nc + j] * cw_c[l];
        }
      }
      V2 = nc * lambda0_c * V2;

      double V4 = mc * lambda0_c * Rw2_c[j * nc + j];

      cw_new = (V1 - V2) / (pow_Rc[j] + V4);

      // If convergence criteria are met, break the loop
      // for (int k = 0; k < nc; ++k){
        diff = fabs(cw_c[j] - cw_new);

        if (min_diff > diff){
          min_diff = diff;
        }

      // }

      if (min_diff <= 1e-8 || min_diff > 5) {
        break;
      }

      // If not convergence, update cw
      cw_c[j] = cw_new;
    }

    if (min_diff <= 1e-8 || min_diff > 5) {
      break;
    }

  } // end outer iteration

  // if ((min_diff > 1e-8 || min_diff > 5) && iter == 0){
  //   memcpy(cw_new, cw_c, nc * sizeof(double));
  // }

  // Calculate c_new
  // for (int i = 0; i < nc; ++i) {
  //   c_new[i] = cw_new[i] * sqrt(sw_c[i]);
  // }

  // result
  double sum3 = 0.0, sum4 = 0.0;
  for (int k = 0; k < mc; ++k) { // iterate by row
    double Rc = 0.0;
    for (int l = 0; l < nc; ++l) { // iterate by col
      Rc += Rw_c[l * mc + k] * cw_c[l];   // /////////// k와 l 순서 바꾸기
    }
    sum3 += (zw_c[k] - Rc) * sw_c[k];
    sum4 += sw_c[k];
  }

  double b_new = sum3 / sum4;

  // Set values in result SEXP
  SET_VECTOR_ELT(result, 0, allocVector(REALSXP, nc));
  SET_VECTOR_ELT(result, 1, ScalarReal(b_new));
  // SET_VECTOR_ELT(result, 2, allocVector(REALSXP, nc));
  // SET_VECTOR_ELT(result, 2, zw);
  // SET_VECTOR_ELT(result, 3, sw);

  // Copy values to result SEXP
  for (int i = 0; i < nc; ++i) {
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


SEXP glm_theta_step(SEXP Gw, SEXP uw, SEXP n, SEXP d, SEXP theta, SEXP lambda_theta, SEXP gamma) {
  // Convert R vectors to C arrays
  int nc = INTEGER(n)[0]; // Extract the value of n
  int dc = INTEGER(d)[0]; // Extract the value of d
  double *uw_c = REAL(uw);
  double *Gw_c = REAL(Gw);
  double *theta_c = REAL(theta);
  double lambda_theta_c = REAL(lambda_theta)[0];
  double gamma_c = REAL(gamma)[0];
  double r = lambda_theta_c * gamma_c * nc;
  // SEXP out = PROTECT(allocVector(VECSXP, 3));

  // Define variables
  // double *theta_new = (double *)malloc(dc * sizeof(double));
  double *pow_theta = (double *)malloc(dc * sizeof(double));
  double theta_new = 0;

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
  double min_diff = 100;
  double diff;

  for(iter = 0; iter < 80; ++iter) {
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
      V1 = 2 * V1;

      if(V1 > 0 && r < fabs(V1)) {
        theta_new = V1 / (pow_theta[j] + nc * lambda_theta_c * (1-gamma_c)) / 2;
      } else {
        theta_new = 0;
      }

      // If convergence criteria are met, break the loop
      diff = fabs(theta_c[j] - theta_new);

      if (min_diff > diff){
        min_diff = diff;
      }

      // if (min_diff <= 1e-8 || min_diff > 10) {
      //   break;
      // }
      // }

      if (min_diff <= 1e-8) {
        break;
      }

      // If not convergence, Update theta_new.
      theta_c[j] = theta_new;
    }

    if (min_diff <= 1e-8) {
      break;
    }
  } // end outer iteration

  if (min_diff > 1e-8 && iter == 0){
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


  // result
  SEXP theta_new_r = PROTECT(allocVector(REALSXP, dc));
  for(int j = 0; j < dc; ++j) {
    REAL(theta_new_r)[j] = theta_c[j];
  }

  free(pow_theta);

  UNPROTECT(1);

  return theta_new_r;
}
