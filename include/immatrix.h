#ifndef IMMATRIX_H
#define IMMATRIX_H 1

#include "gsl/gsl_matrix.h"
#include <filesystem>

bool load_channel_to_matrix(gsl_matrix *&out, const std::filesystem::path imfile, const int channel);
bool save_matrix_to_image(const gsl_matrix *immatrix, const std::filesystem::path imfile);

#endif 
