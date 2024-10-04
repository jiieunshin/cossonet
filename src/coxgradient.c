#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Helper function to allocate a 3D array as a flat 1D array
double* alloc3DArray_flat(int dim1, int dim2, int dim3) {
  return (double*) calloc(dim1 * dim2 * dim3, sizeof(double));
}

// Define the gradient_Hessian_C
SEXP gradient_Hessian_C(SEXP initC, SEXP n, SEXP m, SEXP nc_rs, SEXP eta, SEXP Rtheta1, SEXP Rtheta2, SEXP time,
                        SEXP status, SEXP mscale, SEXP lambda0,
                        SEXP riskset, SEXP tie_size) {

  // Extract input parameters
  int nc = INTEGER(n)[0]; // number of samples
  int mc = INTEGER(m)[0]; // number of coefficients
  int n_riskset = INTEGER(nc_rs)[0]; // number of risk sets

  double *r1 = REAL(Rtheta1); // Rtheta1 is mc x nc (column-major)
  double *r2 = REAL(Rtheta2); // Rtheta2 is mc x mc (column-major)
  double *init_c = REAL(initC); // initC is mc
  int *tie_size_c = INTEGER(tie_size); // tie_size is n_riskset
  int *riskset_c = INTEGER(riskset); // riskset is nc x n_riskset (column-major)
  double *eta_c = REAL(eta); // eta is nc

  // Grad.Term1 = -t(Rtheta1) %*% status / n
  SEXP Grad_Term1_SEXP = PROTECT(allocVector(REALSXP, mc));
  double *Grad_Term1 = REAL(Grad_Term1_SEXP);

  for (int j = 0; j < mc; j++) {
    Grad_Term1[j] = 0.0;
  }

  int *status_ptr = INTEGER(status); // status is nc

  for (int j = 0; j < mc; j++) {   // col
    for (int i = 0; i < nc; i++) {     // row
      Grad_Term1[j] -= r1[i + j * nc] * status_ptr[i] / (double) nc; // Rtheta1[j, i] * status[i] / n
    }
  }

  // Grad.FullNumer = t(Rtheta1) %*% diag(exp(eta))
  // Allocate Grad.FullNumer as mc x nc
  SEXP Grad_FullNumer_SEXP = PROTECT(allocMatrix(REALSXP, mc, nc));
  double *grad_fullnumer = REAL(Grad_FullNumer_SEXP);

  // Compute exp(eta)
  double *exp_eta = (double *) R_alloc(nc, sizeof(double));
  for (int i = 0; i < nc; i++) {
    exp_eta[i] = exp(eta_c[i]);
  }

  // Compute Grad.FullNumer = t(Rtheta1) %*% diag(exp(eta))
  for (int j = 0; j < mc; j++) { // columns
    for (int i = 0; i < nc; i++) { // rows
      grad_fullnumer[i + j * nc] = r1[i + j * nc] * exp_eta[i];
    }
  }

  // Grad.Term3 = 2 * lambda0 * Rtheta2 %*% initC
  SEXP Grad_Term3_SEXP = PROTECT(allocVector(REALSXP, mc));
  double *Grad_Term3 = REAL(Grad_Term3_SEXP);
  double lambda0_val = REAL(lambda0)[0];

  for (int i = 0; i < mc; i++) {
    Grad_Term3[i] = 0.0;
    for (int j = 0; j < mc; j++) {
      Grad_Term3[i] += 2.0 * lambda0_val * r2[i + j * mc] * init_c[j];
    }
  }

  // Grad.FullDenom = Hess.FullDenom = exp(eta)
  SEXP Grad_FullDenom_SEXP = PROTECT(allocVector(REALSXP, nc));
  SEXP Hess_FullDenom_SEXP = PROTECT(allocVector(REALSXP, nc));
  double *grad_fulldenom = REAL(Grad_FullDenom_SEXP);
  double *hess_fulldenom = REAL(Hess_FullDenom_SEXP);

  for (int j = 0; j < nc; j++) {
    grad_fulldenom[j] = exp_eta[j];
    hess_fulldenom[j] = exp_eta[j];
  }

  // Hess.FullNumer = Hess.FullNumer.unScale * array(rep(exp(eta), each = mc^2), dim = c(mc, mc, nc))
  // Allocate Hess.FullNumer as mc x mc x nc (flattened as mc * mc * nc)
  SEXP Hess_FullNumer_SEXP = PROTECT(allocVector(REALSXP, mc * mc * nc));
  double *hess_fullnumer = REAL(Hess_FullNumer_SEXP);


  // Allocate Hess_FullNumer_unScale as mc x mc x nc (flattened)
  SEXP Hess_FullNumer_unScale_SEXP = PROTECT(allocVector(REALSXP, mc * mc * nc));
  double *hess_fullnumer_unScale = REAL(Hess_FullNumer_unScale_SEXP);

  // Initialize Hess_FullNumer_unScale as identity matrices for each slice k
  for (int k = 0; k < nc; k++) {
    for (int a = 0; a < mc; a++) {
      for (int b = 0; b < mc; b++) {
        if (a == b) {
          hess_fullnumer_unScale[a + b * mc + k * mc * mc] = 1.0;
        } else {
          hess_fullnumer_unScale[a + b * mc + k * mc * mc] = 0.0;
        }
      }
    }
  }

  // Compute Hess.FullNumer = Hess.FullNumer.unScale * exp(eta) for each slice k
  for (int k = 0; k < nc; k++) {
    for (int a = 0; a < mc; a++) {
      for (int b = 0; b < mc; b++) {
        hess_fullnumer[a + b * mc + k * mc * mc] = hess_fullnumer_unScale[a + b * mc + k * mc * mc] * exp_eta[k];
      }
    }
  }

  // Hess.Term1 = Hess.Term2 = array(NA, dim = c(mc, mc, n_riskset))
  // Allocate Hess_Term1 and Hess_Term2 as mc x mc x n_riskset (flattened)
  SEXP Hess_Term1_SEXP = PROTECT(allocVector(REALSXP, mc * mc * n_riskset));
  SEXP Hess_Term2_SEXP = PROTECT(allocVector(REALSXP, mc * mc * n_riskset));
  double *hess_term1 = REAL(Hess_Term1_SEXP);
  double *hess_term2 = REAL(Hess_Term2_SEXP);

  // Initialize Hess.Term1 and Hess.Term2 to zero
  memset(hess_term1, 0, sizeof(double) * mc * mc * n_riskset);
  memset(hess_term2, 0, sizeof(double) * mc * mc * n_riskset);

  // temp.Hessian.numer = apply(Hess.FullNumer[, , riskset[, k]], c(1, 2), sum, na.rm = TRUE)
  // Allocate temp_Hessian_numer as mc x mc matrix
  SEXP temp_Hessian_numer_SEXP = PROTECT(allocMatrix(REALSXP, mc, mc));
  double *temp_hessian_numer = REAL(temp_Hessian_numer_SEXP);

  // Initialize temp_hessian_numer to zero
  for (int i = 0; i < mc * mc; i++) {
    temp_hessian_numer[i] = 0.0;
  }

  // Compute initial temp_Hessian_numer and tempSum_exp_eta for k=0
  double tempSum_exp_eta = 0.0;
  for (int i = 0; i < nc; i++) { // for k=0
    int idx = riskset_c[i + 0 * nc] - 1; // 1-based to 0-based
    if (idx >= 0 && idx < nc) {
      tempSum_exp_eta += exp_eta[idx];
      // Accumulate Hess.FullNumer into temp_Hessian_numer
      for (int a = 0; a < mc; a++) {
        for (int b = 0; b < mc; b++) {
          temp_hessian_numer[a + b * mc] += hess_fullnumer[a + b * mc + idx * mc * mc];
        }
      }
    }
  }

  // Grad.Term2 = matrix(NA, ncol = ncol(riskset), nrow = length(initC))
  // Allocate Grad_Term2 as mc x n_riskset
  SEXP temp_Gradient_numer = PROTECT(allocVector(REALSXP, mc));
  double *temp_gradient_numer = REAL(temp_Gradient_numer);

  // 초기화
  for (int i = 0; i < mc; i++) {
    temp_gradient_numer[i] = 0.0;
  }

  SEXP Grad_Term2_SEXP = PROTECT(allocMatrix(REALSXP, mc, n_riskset));
  double *grad_term2 = REAL(Grad_Term2_SEXP);

  // Grad.Term2[, 0] = tie.size[0] * temp.Gradient.numer / tempSum_exp_eta
  for (int j = 0; j < mc; j++) {
    grad_term2[j + 0 * mc] = (tie_size_c[0] * temp_gradient_numer[j]) / tempSum_exp_eta;
  }

  // Hess.Term1[, , 0] = temp_Hessian_numer / tempSum_exp_eta
  for (int a = 0; a < mc; a++) {
    for (int b = 0; b < mc; b++) {
      hess_term1[a + b * mc + 0 * mc * mc] = temp_hessian_numer[a + b * mc] / tempSum_exp_eta;
    }
  }

  // Hess.Term2[, , 0] = 1 / tie.size[0] * Grad.Term2[, 0] %*% t(Grad.Term2[, 0])
  if (tie_size_c[0] > 0) {
    for (int a = 0; a < mc; a++) {
      for (int b = 0; b < mc; b++) {
        hess_term2[a + b * mc + 0 * mc * mc] = (grad_term2[a + 0 * mc] * grad_term2[b + 0 * mc]) / tie_size_c[0];
      }
    }
  }

  // Loop for k=1 to n_riskset-1
  for (int k = 1; k < n_riskset; k++) {
    // Identify excludeID: riskset[, k-1] not in riskset[, k]
    // Implement via a boolean array to mark inclusion
    int *included = (int *) calloc(nc, sizeof(int));
    if (included == NULL) {
      error("Memory allocation failed for included");
    }

    // Mark all IDs in riskset[, k] as included
    for (int i = 0; i < nc; i++) {
      int current_id = riskset_c[i + k * nc] - 1;
      if (current_id >= 0 && current_id < nc) {
        included[current_id] = 1;
      }
    }

    // Identify excludeID: IDs in riskset[, k-1] not in riskset[, k]
    int exclude_ids[mc];
    int exclude_count = 0;
    for (int i = 0; i < nc; i++) {
      int prev_id = riskset_c[i + (k - 1) * nc] - 1;
      if (prev_id >= 0 && prev_id < nc && !included[prev_id]) {
        exclude_ids[exclude_count++] = prev_id;
      }
    }

    // Update tempSum_exp_eta by subtracting sum(exp(eta[excludeID]))
    for (int i = 0; i < exclude_count; i++) {
      int idx = exclude_ids[i];
      tempSum_exp_eta -= exp_eta[idx];
    }

    // Update temp_gradient_numer by subtracting sum(Grad.FullNumer[, excludeID])
    for (int i = 0; i < exclude_count; i++) {
      int idx = exclude_ids[i];
      for (int j = 0; j < mc; j++) {
        temp_gradient_numer[j] -= grad_fullnumer[j + idx * mc];
      }
    }

    // Update temp_hessian_numer by subtracting sum(Hess.FullNumer[, , excludeID])
    for (int i = 0; i < exclude_count; i++) {
      int idx = exclude_ids[i];
      for (int a = 0; a < mc; a++) {
        for (int b = 0; b < mc; b++) {
          temp_hessian_numer[a + b * mc] -= hess_fullnumer[a + b * mc + idx * mc * mc];
        }
      }
    }

    // Grad.Term2[, k] = tie.size[k] * temp.Gradient.numer / tempSum_exp_eta
    if (tempSum_exp_eta > 0) {
      for (int j = 0; j < mc; j++) {
        grad_term2[j + k * mc] = (tie_size_c[k] * temp_gradient_numer[j]) / tempSum_exp_eta;
      }
    } else {
      for (int j = 0; j < mc; j++) {
        grad_term2[j + k * mc] = 0.0;
      }
    }

    // Hess.Term1[, , k] = temp_Hessian_numer / tempSum_exp_eta
    if (tempSum_exp_eta > 0) {
      for (int a = 0; a < mc; a++) {
        for (int b = 0; b < mc; b++) {
          hess_term1[a + b * mc + k * mc * mc] = temp_hessian_numer[a + b * mc] / tempSum_exp_eta;
        }
      }
    } else {
      for (int a = 0; a < mc; a++) {
        for (int b = 0; b < mc; b++) {
          hess_term1[a + b * mc + k * mc * mc] = 0.0;
        }
      }
    }

    // Hess.Term2[, , k] = 1 / tie.size[k] * Grad.Term2[, k] %*% t(Grad.Term2[, k])
    if (tie_size_c[k] > 0) {
      for (int a = 0; a < mc; a++) {
        for (int b = 0; b < mc; b++) {
          hess_term2[a + b * mc + k * mc * mc] = (grad_term2[a + k * mc] * grad_term2[b + k * mc]) / tie_size_c[k];
        }
      }
    } else {
      for (int a = 0; a < mc; a++) {
        for (int b = 0; b < mc; b++) {
          hess_term2[a + b * mc + k * mc * mc] = 0.0;
        }
      }
    }

    // Free included array
    free(included);
  }

  // Grad.Term2 = apply(Grad.Term2, 1, sum) / n
  SEXP Grad_Term2_Sum_SEXP = PROTECT(allocVector(REALSXP, mc));
  double *Grad_Term2_Sum = REAL(Grad_Term2_Sum_SEXP);

  for (int j = 0; j < mc; j++) {
    Grad_Term2_Sum[j] = 0.0;
    for (int k = 0; k < n_riskset; k++) {
      Grad_Term2_Sum[j] += grad_term2[j + k * mc];
    }
    Grad_Term2_Sum[j] /= (double) nc; // /n
  }

  // Gradient = Grad.Term1 + Grad.Term2 + Grad.Term3
  SEXP Gradient_SEXP = PROTECT(allocVector(REALSXP, mc));
  double *Gradient = REAL(Gradient_SEXP);

  for (int j = 0; j < mc; j++) {
    Gradient[j] = Grad_Term1[j] + Grad_Term2_Sum[j] + Grad_Term3[j];
  }

  // Compute Hessian = apply(Hess.Term1, c(1, 2), sum)/n - apply(Hess.Term2, c(1, 2), sum)/n + 2 * lambda0 * Rtheta2
  // Sum Hess.Term1 across the 3rd dimension
  SEXP sum_Hess_Term1_SEXP = PROTECT(allocMatrix(REALSXP, mc, mc));
  double *sum_Hess_Term1 = REAL(sum_Hess_Term1_SEXP);
  memset(sum_Hess_Term1, 0, sizeof(double) * mc * mc);

  for (int k = 0; k < n_riskset; k++) {
    for (int a = 0; a < mc; a++) {
      for (int b = 0; b < mc; b++) {
        sum_Hess_Term1[a + b * mc] += hess_term1[a + b * mc + k * mc * mc];
      }
    }
  }

  for (int i = 0; i < mc * mc; i++) {
    sum_Hess_Term1[i] /= (double) nc; // /n
  }

  // Sum Hess.Term2 across the 3rd dimension
  SEXP sum_Hess_Term2_SEXP = PROTECT(allocMatrix(REALSXP, mc, mc));
  double *sum_Hess_Term2 = REAL(sum_Hess_Term2_SEXP);
  memset(sum_Hess_Term2, 0, sizeof(double) * mc * mc);

  for (int k = 0; k < n_riskset; k++) {
    for (int a = 0; a < mc; a++) {
      for (int b = 0; b < mc; b++) {
        sum_Hess_Term2[a + b * mc] += hess_term2[a + b * mc + k * mc * mc];
      }
    }
  }

  for (int i = 0; i < mc * mc; i++) {
    sum_Hess_Term2[i] /= (double) nc; // /n
  }

  // Compute 2 * lambda0 * Rtheta2
  SEXP term_lambda_Rtheta2_SEXP = PROTECT(allocMatrix(REALSXP, mc, mc));
  double *term_lambda_Rtheta2 = REAL(term_lambda_Rtheta2_SEXP);

  for (int i = 0; i < mc * mc; i++) {
    term_lambda_Rtheta2[i] = 2.0 * lambda0_val * r2[i];
  }

  // Compute Hessian = sum_Hess_Term1 - sum_Hess_Term2 + 2 * lambda0 * Rtheta2
  SEXP Hessian_SEXP = PROTECT(allocMatrix(REALSXP, mc, mc));
  double *Hessian = REAL(Hessian_SEXP);

  for (int i = 0; i < mc * mc; i++) {
    Hessian[i] = sum_Hess_Term1[i] - sum_Hess_Term2[i] + term_lambda_Rtheta2[i];
  }

  // Create a list to return Gradient and Hessian
  SEXP result = PROTECT(allocVector(VECSXP, 2));
  SET_VECTOR_ELT(result, 0, Gradient_SEXP);
  SET_VECTOR_ELT(result, 1, Hessian_SEXP);

  // Set names
  SEXP names = PROTECT(allocVector(STRSXP, 2));
  SET_STRING_ELT(names, 0, mkChar("Gradient"));
  SET_STRING_ELT(names, 1, mkChar("Hessian"));
  setAttrib(result, R_NamesSymbol, names);

  // Clean up
  UNPROTECT(18);

  return result;
}
