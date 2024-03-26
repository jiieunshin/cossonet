#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>


double* cal_expRc(double *R1, int m, int n, double *c_old, double *Risk) {
  double *expRc = (double*)malloc(m * sizeof(double));

  for (int j = 0; j < m; j++) {
    int *ids = (int *)malloc(m * sizeof(int));
    double expRc_j = 0.0;

    for(int k = 0; k < m; k++){
      ids[k] = Risk[j * m + k];
      Rprintf("ids");
      Rprintf("%d\t", ids[k]);
      int exists = 1;
      for(int l = 0; l < k; l++){
        if(ids[k] == ids[l] || ids[k] == 0){
          exists = 0;
          break;
        }
      }

      double add = 0.0;
      if(exists == 1){
        int id = ids[k];
        Rprintf("id");
          Rprintf("%d\t", id);
        for(int l = 0; l < n; l++){
          add += R1[l * m + id] * c_old[l];
        }
        expRc_j += exp(add);
      } else if(exists == 0){
        expRc_j += 0;
      }
      Rprintf("\n");
    }
    expRc[j] = expRc_j;
  }

  return expRc;
}

// SEXP time, SEXP status,
SEXP Cget_c(SEXP R1, SEXP R2, SEXP n, SEXP m, SEXP Risk, SEXP c_init, SEXP lambda0) {
  int nc = INTEGER(n)[0];
  int mc = INTEGER(m)[0];

  // Convert R vectors to C arrays
  double *R1_c = REAL(R1);
  // double *R2_c = REAL(R2);
  // double *time_c = REAL(time);
  // double *status_c = REAL(status);
  double *Risk_c = REAL(Risk);
  double *c_old = REAL(c_init);
  // double lambda0_c = REAL(lambda0)[0];

  // Define variables
  // double *c_new = (double *)malloc(nc * sizeof(double));
  double *expRc = (double *)malloc(nc * sizeof(double));
  // double *RexpRc = (double *)malloc(mc * sizeof(double));
  // double *RRexpRc = (double *)malloc(mc * sizeof(double));

  // int iter = 0;
  // double max_diff = 1e-6;
  // outer loop
  // for (iter = 0; iter < 20; ++iter) {

    expRc = cal_expRc(R1_c, mc, nc, c_old, Risk_c);

    // calculate expRc
    // for (int j = 0; j < mc; j++) {  // iterate by column of Rist set
    //
    // }
  // } // end outer iteration

  // if (max_diff > 1e-6 && iter == 1){
  //   memcpy(cw_new, cw_c, nc * sizeof(double));
  // } else{
  //   double cw_new_sd = sd(cw_new, nc); // cw_new의 표준편차 계산
  //   for (int i = 0; i < nc; ++i) {
  //     cw_new[i] /= cw_new_sd; // cw_new를 표준편차로 나누어 줌
  //   }
  // }
  //
  // if (max_diff > 1e-6 && iter == 1){
  //   theta_new = (double *)malloc(dc * sizeof(double));
  // }
  //
  // // result
  SEXP expRc_r = PROTECT(allocVector(REALSXP, nc));
  for(int j = 0; j < nc; ++j) {
    REAL(expRc_r)[j] = expRc[j];
  }

  free(expRc);

  UNPROTECT(1);

  return expRc_r;
}

