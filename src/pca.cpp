#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>

// Perform principle component analysis on a dataset X. This proceeds
// via SVD. X will be an NxD matrix, where N >= D, and where each
// datapoint is formatted as a row of X. The result will be a DxD
// matrix, where each column represents a basis vector for the
// dataset. The centre of the data will be stored in a vector at
// centre. The eigenvalues will be stored in a vector at
// eigenvals. The returned matrix and returned vectors are allocated
// on the heap, and it is the responsibility of the caller to free
// these using gsl_matrix_free and gsl_vector_free.
gsl_matrix *do_pca(const gsl_matrix *X, gsl_vector *&centre, gsl_vector *&eigenvals) {
  // centre the data
  const size_t N = X->size1;
  const size_t D = X->size2;
  gsl_matrix *U = gsl_matrix_alloc(N, D);
  gsl_matrix_memcpy(U, X);
  centre = gsl_vector_alloc(D);
  for (int col = 0; col < D; ++col) {
    // gsl_vector_sum not in 2.5, so sum manually
    double s = 0;
    for (int row = 0; row < N; ++row) {
      s += gsl_matrix_get(X, row, col);
    }
    s /= N;
    gsl_vector_set(centre, col, s);
    gsl_vector_view col_v = gsl_matrix_column(U, col);
    gsl_vector_add_constant(&col_v.vector, -s);
  }

  // apply SVD to the data
  // GSL's SVD requires N >= D, and performs the thin SVD.
  // The result is that V will be the DxD eigenvector matrix.
  gsl_matrix *eigenvecs = gsl_matrix_alloc(D, D);
  eigenvals = gsl_vector_calloc(D);
  gsl_vector *work = gsl_vector_alloc(D);
  gsl_linalg_SV_decomp(U, eigenvecs, eigenvals, work);
  gsl_matrix_free(U);
  gsl_vector_free(work);

  // V will have the eigenvectors as rows because of how we are doing
  // this calculation, so we can transpose it here.
  gsl_matrix_transpose(eigenvecs);

  return eigenvecs;
}

// provided for compatibility with earlier code
gsl_matrix *do_pca(const gsl_matrix *X, gsl_vector *&centre) {
  gsl_vector *eigenvals;
  gsl_matrix *eigenvecs = do_pca(X, centre, eigenvals);
  // discard the unwanted eigenvals
  gsl_vector_free(eigenvals);
  return eigenvecs;
}
