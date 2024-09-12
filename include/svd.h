#ifndef SVD_H
#define SVD_H 1

#include <gsl/gsl_matrix.h>

void get_svd_dims(gsl_matrix *X, size_t &u1, size_t &u2, size_t &v1, size_t &v2, size_t &s1, size_t &s2);
void do_svd(gsl_matrix *X, gsl_matrix *&U, gsl_matrix *&V, gsl_matrix *&S);

#endif 
