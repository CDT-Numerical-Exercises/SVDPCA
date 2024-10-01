#include <filesystem>
#include <algorithm>
#include "immatrix.h"
#include "gsl/gsl_matrix.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb/stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

double clamp(double x, double min, double max) {
  if (x < min) return min;
  if (x > max) return max;
  return x;
}

bool load_channel_to_matrix(gsl_matrix *&out, const std::filesystem::path imfile, const int channel) {
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
  unsigned char *im_array = stbi_load(filename, &x, &y, &n, n);

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
      gsl_matrix_set(out, i, j, im_array[n*(x*i + j) + channel]/255.);
    }
  }

  // free stb's version of the image
  stbi_image_free(im_array);

  return true;
}

bool save_matrix_to_image(const gsl_matrix *immatrix, const std::filesystem::path imfile) {
  const int w = immatrix->size2;
  const int h = immatrix->size1;
  const int comp = 1;

  // convert the matrix to an array
  // we need to do this manually as we have to convert from float to char
  unsigned char *data = new unsigned char[w*h];
  for (int i = 0; i < h; ++i) {
    for (int j = 0; j < w; ++j) {
      const int pixel = w*i + j;
      data[pixel] = (unsigned char)(255*clamp(gsl_matrix_get(immatrix, i, j), 0, 1));
    }
  }

  std::string ext = imfile.extension();
  // convert to uppercase
  std::transform(ext.begin(), ext.end(), ext.begin(), ::toupper);

  const std::string fn = imfile.string();
  int outcode = 0;
  if (ext == ".PNG") {
    outcode = stbi_write_png(fn.c_str(), w, h, comp, data, 0);
  } else if (ext == ".JPG") {
    outcode = stbi_write_jpg(fn.c_str(), w, h, comp, data, 100);
  } else if (ext == ".BMP") {
    outcode = stbi_write_bmp(fn.c_str(), w, h, comp, data);
  } else {
    return false;
  }

  delete[] data;

  return true;
}
