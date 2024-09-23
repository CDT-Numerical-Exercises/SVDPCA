#include <filesystem>
#include "immatrix.h"
#include "gsl/gsl_matrix.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"

bool load_channel_to_matrix(gsl_matrix *&out, std::filesystem::path imfile, int channel) {
  // get the details of the image
  // we can use this to verify the choice of channel
  int x,y,n,ok;
  const std::string fn_string = imfile.string(); 
  const char *filename = fn_string.c_str(); 
  ok = stbi_info(filename, &x, &y, &n);
  if (!ok) {
    return false;
  }
  if (channel > n - 1) {
    return false;
  }

  // load the full image as a heap-allocated float array
  float *im_array = stbi_loadf(filename, &x, &y, &n, n);

  // copy the contents into a matrix
  // I don't think there's any clever way of doing this in GSL as the
  // different image channels are stored interleaved, so dropping a
  // channel from the matrix is not so straightforward
  out = gsl_matrix_alloc(y, x);
  for (int i = 0; i < y; ++i) {
    for (int j = 0; j < x; ++j) {
      // for n = 4:
      //  i,j=0 -> 4*0 + c
      //  i=0, j=1 -> 4*1 + c
      //  i=1, j=0 -> 4*width + c
      gsl_matrix_set(out, i, j, im_array[n*(x*i + j) + channel]);
    }
  }

  // free stb's version of the image
  stbi_image_free(im_array);

  return true;
}
