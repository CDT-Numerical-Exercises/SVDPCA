#ifndef PCA_H
#define PCA_H 1

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

gsl_matrix *do_pca(const gsl_matrix *X, gsl_vector *&centre);

#endif
