#ifndef PCA_H
#define PCA_H 1

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

enum DataFormat {
  RowVector,
  ColumnVector
};

gsl_matrix *do_pca(const gsl_matrix *X, DataFormat format, gsl_vector *&centre, gsl_vector *&eigenvals);
gsl_vector *pca_project(const gsl_matrix *eigenvecs, const gsl_vector *X, const int dims);

#endif
