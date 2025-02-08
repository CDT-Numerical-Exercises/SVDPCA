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
#include "csv.h"

constexpr int KEEP_COMPONENTS = 2;

int main() {
  // data is stored as row vectors, i.e.
  //   [x0, x1, x2...]
  //   [x0, x1, x2...]
  // which is the format we need
  gsl_matrix *X = load_csv_to_dmatrix("svdpca-problem1-data.csv");

  gsl_vector *centre, *eigenvals;
  gsl_matrix *eigenvecs = do_pca(X, RowVector, centre, eigenvals);
  if (eigenvecs == NULL) {
    std::cout << "PCA failed" << std::endl;
    return -1;
  }

  // turn the eigenvals into percentages of the variance
  double s2_sum = 0;
  for (int i = 0; i < eigenvals->size; ++i) {
    double s = gsl_vector_get(eigenvals, i);
    s2_sum += s*s;
  }

  // output the eigenvectors
  std::cout << "== Best Eigenvectors ==" << std::endl;
  for (int col = 0; col < eigenvecs->size2; ++col) {
    if (col == KEEP_COMPONENTS) {
      std::cout << "== Discarded Eigenvectors ==" << std::endl;
    }
    gsl_vector_view view = gsl_matrix_column(eigenvecs, col);
    std::cout << "Component " << col + 1 << ": ";
    print_vector(&view.vector);
  }

  // plot the scree plot
  {
  Gnuplot gp;
  gp << "plot '-' with linespoints\n";
  for (int i = 0; i < eigenvals->size; ++i) {
    double s = gsl_vector_get(eigenvals, i);
    gp << i << " " << s * s / s2_sum << "\n";
  }
  gp << "e\n";
  }

  // project each datum into 2D and plot it
  Gnuplot gp;
  gp << "set size ratio -1\n";
  gp << "plot '-' with points\n";
  // we want to keep the first two eigenvectors
  for (int i = 0; i < X->size1; ++i) {
    // centre the datum
    const gsl_vector_view row = gsl_matrix_row(X, i);
    gsl_vector *Xc = gsl_vector_alloc(X->size2);
    gsl_vector_memcpy(Xc, &row.vector);
    gsl_vector_sub(Xc, centre);

    // find the projection
    gsl_vector *proj = pca_project(eigenvecs, Xc, KEEP_COMPONENTS);

    // output
    for (int j = 0; j < proj->size - 1; ++j) {
      gp << gsl_vector_get(proj, j) << " ";
    }
    gp << gsl_vector_get(proj, proj->size-1) << "\n";

    // clean up
    gsl_vector_free(proj);
    gsl_vector_free(Xc);
  }
  gp << "e\n";

  // clean up
  gsl_matrix_free(X);
  gsl_matrix_free(eigenvecs);
  gsl_vector_free(centre);
  gsl_vector_free(eigenvals);
}
