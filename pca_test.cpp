#include <iostream>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "gnuplot-iostream/gnuplot-iostream.h"

#include "pca.h"
#include "helpers.h"

gsl_matrix *gen_test_data(const int points) {
  constexpr double means[2] = {0, 0};
  constexpr double stddevs[2] = {1, 2};
  // each datum (x,y) is stored as a row to make it easier to apply
  // SVD
  gsl_matrix *M = gsl_matrix_alloc(points, 2);
  for (int i = 0; i < M->size1; ++i) {
    for (int j = 0; j < 2; ++j) {
      gsl_matrix_set(M, i, j, norm(means[j], stddevs[j]));
    }
  }

  // apply a rotation
  // ordinarily, we would do X' = R@X
  // this assumes X uses column vectors
  // we must instead do X' = X@R^T
  const double theta = M_PI/4+M_PI/8;
  double R_data[4] = {
    cos(theta), sin(theta),
    -sin(theta), cos(theta)
  };
  gsl_matrix_view R_view = gsl_matrix_view_array(R_data, 2, 2);

  gsl_matrix *X = gsl_matrix_alloc(points, 2);
  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, M, &R_view.matrix, 0, X);
  gsl_matrix_free(M);

  // apply the shifts
  constexpr double shifts[2] = { 1, 0.5 };
  gsl_vector_view x = gsl_matrix_column(X, 0);
  gsl_vector_view y = gsl_matrix_column(X, 1);
  gsl_vector_add_constant(&x.vector, shifts[0]);
  gsl_vector_add_constant(&y.vector, shifts[1]);

  return X;
}

int main() {
  gsl_matrix *X = gen_test_data(2000);

  gsl_vector *centre, *eigenvals;
  gsl_matrix *eigenvecs = do_pca(X, RowVector, centre, eigenvals);
  print_matrix(eigenvecs);

  // scale the eigenvectors by S
  // I'm using GSL 2.5, so gsl_matrix_scale_columns is not available.
  // we can construct a matrix as diag(eigenvals) and multiply by it
  gsl_matrix *scaled_vecs = gsl_matrix_alloc(eigenvecs->size1, eigenvecs->size2);
  gsl_matrix *S = gsl_matrix_calloc(eigenvals->size, eigenvals->size);
  gsl_vector_view s_view = gsl_matrix_diagonal(S);
  gsl_vector_memcpy(&s_view.vector, eigenvals);
  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, eigenvecs, S, 0, scaled_vecs);

  // free up the bits we no longer care about
  gsl_matrix_free(eigenvecs);
  gsl_matrix_free(S);
  gsl_vector_free(eigenvals);

  // scale everything to make the arrows visible
  gsl_vector_view vec1 = gsl_matrix_column(scaled_vecs, 0);
  gsl_matrix_scale(scaled_vecs, 3/gsl_vector_length(&vec1.vector));

  // draw the vectors
  Gnuplot gp;
  gp << "set xrange [-8:8]\nset yrange [-8:8]\n";
  gp << "set size square\n";
  gp << "set style arrow 2 front lc 'red' lw 4\n";
  gp << "set style arrow 1 front lc 'blue' lw 4\n";
  
  gsl_vector *endpoint = gsl_vector_alloc(scaled_vecs->size1);
  for (int col = 0; col < scaled_vecs->size2; ++col) {
    gsl_vector_view eigenvec = gsl_matrix_column(scaled_vecs, col);
    gsl_vector_memcpy(endpoint, centre);
    gsl_vector_add(endpoint, &eigenvec.vector);
    gp << "set arrow as " << col+1 << " from " << gsl_vector_get(centre, 0) << "," << gsl_vector_get(centre, 1);
    gp << " to " << gsl_vector_get(endpoint, 0) << "," << gsl_vector_get(endpoint, 1) << "\n";
  }
  gsl_vector_free(endpoint);
  gsl_matrix_free(scaled_vecs);
  gsl_vector_free(centre);

  gp << "plot '-' with points pt 7\n";
  for (int i = 0; i < X->size1; ++i) {
    for (int j = 0; j < 2; ++j) {
      gp << gsl_matrix_get(X, i, j) << " ";
    }
    gp << "\n";
  }
  // gp << "set arrow from 0,0 to 1,1\n";
  // gp << "plot NaN t ''\n";

  gsl_matrix_free(X);

  return 0;
}
