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
#include "immatrix.h"

constexpr int SUBIM_ROWS = 2;
constexpr int SUBIM_COLS = 5;
constexpr int SUBIMS = SUBIM_ROWS * SUBIM_COLS;
constexpr int SUBIM_WIDTH = 128;
constexpr int SUBIM_HEIGHT = 155;
constexpr int SUBIM_PIXELS = SUBIM_WIDTH * SUBIM_HEIGHT;

int main() {
  gsl_matrix *im_matrix_raw;
  if (!load_channel_to_matrix(im_matrix_raw, "Faces.png", 0)) {
    std::cout << "Error loading image file!" << std::endl;
    return -1;
  }

  // the images at the top are 156 pixels tall, whereas at the bottom
  // they are 155 tall
  // let's crop off the first row of pixels
  gsl_matrix_view im_matrix_view = gsl_matrix_submatrix(im_matrix_raw, 1, 0, SUBIM_HEIGHT*SUBIM_ROWS, SUBIM_WIDTH*SUBIM_COLS);
  gsl_matrix *im_matrix = &im_matrix_view.matrix;

  // reformat the data into a bunch of row vectors, one for each image
  gsl_matrix *X = gsl_matrix_alloc(SUBIMS, SUBIM_PIXELS);
  for (int y = 0; y < SUBIM_ROWS; ++y) {
    for (int x = 0; x < SUBIM_COLS; ++x) {
      const int subim = y*SUBIM_COLS + x;
      const int topleft_x = x*SUBIM_WIDTH;
      const int topleft_y = y*SUBIM_HEIGHT;
      std::cout << "Getting image " << subim << " @ (i,j) = " << topleft_y << " " << topleft_x << std::endl;
      const gsl_matrix_const_view im_view = gsl_matrix_const_submatrix(im_matrix, topleft_y, topleft_x, SUBIM_HEIGHT, SUBIM_WIDTH);
      for (int i = 0; i < SUBIM_HEIGHT; ++i) {
        for (int j = 0; j < SUBIM_WIDTH; ++j) {
          const int pixel = SUBIM_WIDTH*i + j;
          gsl_matrix_set(X, subim, pixel, gsl_matrix_get(&im_view.matrix, i, j));
        }
      }
    }
  }

  // dealloc the original matrix now, since it's redundant
  gsl_matrix_free(im_matrix_raw);
}
