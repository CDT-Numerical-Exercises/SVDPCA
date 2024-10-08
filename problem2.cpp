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

void write_image(std::filesystem::path out_file, gsl_vector *flat, int width, int height) {
  // convert to a matrix via a view
  gsl_matrix_view view = gsl_matrix_view_vector(flat, height, width);
  save_matrix_to_image(&view.matrix, out_file);
}

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

  // reformat the data into a bunch of column vectors, one for each image
  // we have to use column vectors to extract the correct vectors
  // otherwise, the SVD function will not be compatible with the data
  gsl_matrix *X = gsl_matrix_alloc(SUBIM_PIXELS, SUBIMS);
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
          gsl_matrix_set(X, pixel, subim, gsl_matrix_get(&im_view.matrix, i, j));
        }
      }
    }
  }

  // dealloc the original matrix now, since it's redundant
  // im_matrix_view and im_matrix are now dead
  gsl_matrix_free(im_matrix_raw);

  // do the principle component analysis
  // we will use the components later, but the centre also tells us
  // the mean vector!
  gsl_vector *centre, *eigenvals;
  gsl_matrix *eigenvecs = do_pca(X, ColumnVector, centre, eigenvals);
  if (eigenvecs == NULL) {
    std::cout << "PCA failed" << std::endl;
    return -1;
  }

  std::cout << "Got " << eigenvecs->size2 << " eigenvecs of length " << eigenvecs->size1 << std::endl;
  std::cout << "Got " << eigenvals->size << " eigenvalues" << std::endl;
  std::cout << "Centre length: " << centre->size << std::endl;
  std::cout << "Centre max: " << gsl_vector_max(centre) << std::endl;
  std::cout << "Centre min: " << gsl_vector_min(centre) << std::endl;
  // save the mean vector as an image
  write_image("x0.png", centre, SUBIM_WIDTH, SUBIM_HEIGHT);

  // create the scree plot
  double s2_sum = 0;
  for (int i = 0; i < eigenvals->size; ++i) {
    const double s = gsl_vector_get(eigenvals, i);
    s2_sum += s*s;
  }
  Gnuplot gp;
  gp << "plot '-' using 1:2 with linespoints title 'eigenval', '-' using 1:3 with linespoints title 'cumulative'\n";
  for (int j = 0; j < 2; ++j) {
    double running = 0.0;
    for (int i = 0; i < eigenvals->size; ++i) {
      const double s = gsl_vector_get(eigenvals, i);
      const double pc = s * s / s2_sum;
      running += pc;
      gp << i << " " << pc << " " << running << "\n";
    }
    gp << "E\n";
  }

  // create projections of the faces and export them
  const int dims[2] = { 5, 10 };
  constexpr int len = sizeof(dims)/sizeof(dims[0]);
  constexpr int bufsize = 100;
  char buf[bufsize];
  for (int i = 0; i < len; ++i) {
    const int D = dims[i];
    std::cout << D << std::endl;
    for (int face = 0; face < X->size2; ++face) {
      gsl_vector_view view = gsl_matrix_column(X, face);
      // centre the data
      gsl_vector *Xc = gsl_vector_alloc(view.vector.size);
      gsl_vector_memcpy(Xc, &view.vector);
      gsl_vector_sub(Xc, centre);

      // find the projection
      gsl_vector *proj = pca_project(eigenvecs, Xc, D);
      gsl_vector_add(proj, centre);

      // output
      snprintf(buf, bufsize, "face%d_%d.png", face, D);
      write_image(buf, proj, SUBIM_WIDTH, SUBIM_HEIGHT);
      gsl_vector_free(proj);
      gsl_vector_free(Xc);
      std::cout << "Reduced face " << face << " to " << D << " dimensions" << std::endl;
    }
  }

  gsl_vector_free(centre);
  gsl_matrix_free(eigenvecs);
  gsl_vector_free(eigenvals);
  gsl_matrix_free(X);
}
