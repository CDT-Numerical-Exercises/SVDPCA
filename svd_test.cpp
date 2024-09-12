#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#include "svd.h"
#include "helpers.h"

int main() {
  double X_vals[8] = {
    2, 4,
    1, 3,
    1, 1,
    2, 2
  };
  gsl_matrix_view X_view = gsl_matrix_view_array(X_vals, 4, 2);
  gsl_matrix *X = &X_view.matrix;

  gsl_matrix *U, *V, *S;
  do_svd(X, U, V, S);

  gsl_matrix *SV = gsl_matrix_alloc(S->size1, V->size1);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, S, V, 0, SV);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, U, SV, 0, S); // overwrite S to avoid allocating a new matrix

  std::cout << "Result of USV^T:" << std::endl;
  print_matrix(S, 9);

  gsl_matrix_free(SV);
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_matrix_free(S);

  // compare with GSL's implementation
  gsl_vector *work = gsl_vector_alloc(X->size2);
  U = gsl_matrix_calloc(X->size1, X->size2);
  gsl_matrix_memcpy(U, X);
  S = gsl_matrix_calloc(X->size2, X->size2);
  gsl_vector_view s_diag_view = gsl_matrix_diagonal(S);
  gsl_vector *s = &s_diag_view.vector;
  V = gsl_matrix_calloc(X->size2, X->size2);

  // can't just overwrite S now as S will be the wrong shape
  gsl_matrix *USV = gsl_matrix_calloc(X->size1, X->size2);

  gsl_linalg_SV_decomp(U, V, s, work);
  gsl_vector_free(work);

  SV = gsl_matrix_calloc(S->size1, V->size1);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, S, V, 0, SV);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, U, SV, 0, USV);

  std::cout << "Result of USV^T (GSL implementation):" << std::endl;
  print_matrix(USV, 2);

  gsl_matrix_free(USV);
  gsl_matrix_free(SV);
  gsl_matrix_free(U);
  gsl_matrix_free(V);
  gsl_matrix_free(S);
  
  return 0;
}
