#ifndef IMMATRIX_H
#define IMMATRIX_H 1

#include "gsl/gsl_matrix.h"
#include <filesystem>

bool load_channel_to_matrix(gsl_matrix *&out, std::filesystem::path imfile, int channel);

#endif 
