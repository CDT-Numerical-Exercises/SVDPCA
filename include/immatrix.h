#ifndef IMMATRIX_H
#define IMMATRIX_H 1

#include "gsl/gsl_matrix.h"
#include <filesystem>

bool load_channel_to_matrix(gsl_matrix *&out, std::filesystem::path imfile, int channel);
bool save_matrix_to_image(gsl_matrix *immatrix, std::filesystem::path imfile);

#endif 
