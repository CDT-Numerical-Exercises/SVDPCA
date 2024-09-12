#include "svd.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>

#define _USE_MATH_DEFINES
#include <cmath>

#ifdef VERBOSE_SVD
#include <iostream>
#include "helpers.h"
#endif

// determine, from a given input matrix X, the matrix dimensions
// needed to store U, V and S
void get_svd_dims(gsl_matrix *X, size_t &u1, size_t &u2, size_t &v1, size_t &v2, size_t &s1, size_t &s2) {
  // for an DxN matrix:
  //  U is DxD
  //  S is DxN
  //  V is NxN
  u1 = X->size1;
  u2 = X->size2;
  v1 = X->size2;
  v2 = X->size2;
  s1 = X->size2;
  s2 = X->size2;
}

// Perform singular value decomposition for matrix X. The output
// matrices U, V, S are allocated on the heap using GSL's matrix
// malloc, and must be freed manually by the caller using
// gsl_matrix_free.
//
// For simplicity of use, this function allocates working memory on
// the heap as required, and all of this will be cleaned up when
// finished. For an implementation with more control over memory
// usage, use the SVD implementation in GSL's linalg library.
void do_svd(gsl_matrix *X, gsl_matrix *&U, gsl_matrix *&V, gsl_matrix *&S) {
  // calculate A = X @ X^T
  // X is MxN => A is MxM
  gsl_matrix *A = gsl_matrix_alloc(X->size1, X->size1);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, X, X, 0, A);

  #ifdef VERBOSE_SVD
  std::cout << "A: " << std::endl;
  print_matrix(A, 2);
  #endif

  // compute eigenvalues and eigenvectors of A
  // A is mangled in the process, so we can dealloc
  gsl_vector *LU = gsl_vector_alloc(A->size1);
  U = gsl_matrix_alloc(A->size1, A->size1);
  {
    gsl_eigen_symmv_workspace *eig_work_A = gsl_eigen_symmv_alloc(A->size1);
    gsl_eigen_symmv(A, LU, U, eig_work_A);
    gsl_eigen_symmv_free(eig_work_A);
  }
  gsl_matrix_free(A);

  // calculate B = X^T @ X
  // X is MxN => B is NxN
  gsl_matrix *B = gsl_matrix_alloc(X->size2, X->size2);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, X, X, 0, B);

  #ifdef VERBOSE_SVD
  std::cout << "B: " << std::endl;
  print_matrix(B, 2);
  #endif

  // compute eigenvalues and eigenvectors of B
  // B is mangled in the process, so we can dealloc
  gsl_vector *LV = gsl_vector_alloc(B->size1);
  V = gsl_matrix_alloc(B->size2, B->size2);
  {
    gsl_eigen_symmv_workspace *eig_work_B = gsl_eigen_symmv_alloc(B->size1);
    gsl_eigen_symmv(B, LV, V, eig_work_B);
    gsl_eigen_symmv_free(eig_work_B);
  }
  gsl_matrix_free(B);

  // simultaneously sort (LU, U) and (LV, V)
  // sorting should be from largest to smallest (descending)
  gsl_eigen_symmv_sort(LU, U, GSL_EIGEN_SORT_VAL_DESC);
  gsl_eigen_symmv_sort(LV, V, GSL_EIGEN_SORT_VAL_DESC);

  #ifdef VERBOSE_SVD
  std::cout << "LU: ";
  print_vector(LU);
  std::cout << "U: " << std::endl;
  print_matrix(U, 10);
  std::cout << std::endl;
  std::cout << "LV: ";
  print_vector(LV);
  std::cout << "V: " << std::endl;
  print_matrix(V, 10);
  #endif

  // build S from the smallest set of eigenvalues
  S = gsl_matrix_calloc(X->size1, X->size2);
  {
    gsl_vector_view S_diag = gsl_matrix_diagonal(S);
    for (int i = 0; i < S_diag.vector.size; ++i) {
      gsl_vector_set(&S_diag.vector, i, sqrt(gsl_vector_get(LU, i)));
    }
  }
  
  #ifdef VERBOSE_SVD
  std::cout << "S: " << std::endl;
  print_matrix(S, 8);
  #endif

  // all done! time to clean up
  gsl_vector_free(LU);
  gsl_vector_free(LV);
}
