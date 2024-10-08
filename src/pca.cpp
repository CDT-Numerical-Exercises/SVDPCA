#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include "pca.h"

// Perform principle component analysis on a dataset X. This proceeds
// via SVD. X must be an NxD matrix, where N >= D. You may store the
// data either as row vectors or column vectors, however is necessary
// to ensure N >= D. You must specify the data format to retrieve the
// correct data centre and eigenvectors.
//
// The result will be a matrix where each column represents a basis
// vector for the dataset. The centre of the data will be stored in a
// vector at centre. The eigenvalues will be stored in a vector at
// eigenvals.
//
// The returned matrix and returned vectors are allocated on the heap,
// and it is the responsibility of the caller to free these using
// gsl_matrix_free and gsl_vector_free.
gsl_matrix *do_pca(const gsl_matrix *X, const DataFormat format, gsl_vector *&centre, gsl_vector *&eigenvals) {
  // set the right variables dependent on the format
  const size_t N = X->size1;
  const size_t D = X->size2;
  size_t datum_length;
  size_t datapoints;
  switch (format) {
  case RowVector:
    datum_length = D;
    datapoints = N;
    break;
  case ColumnVector:
    datum_length = N;
    datapoints = D;
    break;
  default:
    return NULL; // we don't need a default case after this now
  }

  // copy the data
  gsl_matrix *U = gsl_matrix_alloc(N, D);
  gsl_matrix_memcpy(U, X);

  // determine the centre by summing the indices of each datum and
  // averaging
  centre = gsl_vector_alloc(datum_length);
  for (int i = 0; i < datum_length; ++i) {
    // gsl_vector_sum not in 2.5, so sum manually
    double s = 0;
    for (int datum = 0; datum < datapoints; ++datum) {
      if (format == RowVector) {
        s += gsl_matrix_get(X, datum, i);
      } else {
        s += gsl_matrix_get(X, i, datum);
      }
    }
    s /= datapoints;
    gsl_vector_set(centre, i, s);

    // offset the data to ensure it is central
    gsl_vector_view index_view;
    if (format == RowVector) {
      index_view = gsl_matrix_column(U, i);
    } else {
      index_view = gsl_matrix_row(U, i);
    }
    gsl_vector_add_constant(&index_view.vector, -s);
  }

  // apply SVD to the data
  // GSL's SVD requires N >= D, and performs the thin SVD.
  // If using column vectors, we get that U contains the eigenvectors
  // If using row vectors, we get that V contains the eigenvectors
  gsl_matrix *V = gsl_matrix_alloc(D, D);
  eigenvals = gsl_vector_calloc(D);
  gsl_vector *work = gsl_vector_alloc(D);
  gsl_linalg_SV_decomp(U, V, eigenvals, work);

  // clean up the unneeded parts
  gsl_vector_free(work);
  if (format == RowVector) {
    gsl_matrix_free(U);
    return V;
  }
  gsl_matrix_free(V);
  return U;
}

// Using a complete set of PCA eigenvectors, project the data X to
// reduce its dimensionality. The eigenvectors should be stored as
// column vectors. Returns a heap-allocated vector, which must be
// freed by the caller.
gsl_vector *pca_project(const gsl_matrix *eigenvecs, const gsl_vector *Xc, const int dims) {
  const gsl_matrix_const_view reduced = gsl_matrix_const_submatrix(eigenvecs, 0, 0, eigenvecs->size1, dims);
  const gsl_matrix *U = &reduced.matrix;

  // allocate the return vector
  // this tells us the proportions of each eigenvector
  gsl_vector *A = gsl_vector_alloc(U->size2);

  // go go gadget BLAS
  gsl_blas_dgemv(CblasTrans, 1, U, Xc, 0, A);

  // calculate the projection
  gsl_vector *X_p = gsl_vector_alloc(Xc->size);
  gsl_blas_dgemv(CblasNoTrans, 1, U, A, 0, X_p);

  gsl_vector_free(A);

  return X_p;
}
