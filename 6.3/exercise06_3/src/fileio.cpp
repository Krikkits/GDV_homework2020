/**
 * Computer Graphics Lecture WS 2009/2010 Ulm University
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#include "fileio.h"
#include "rgbe.h"
#include <iostream>
#include <fstream>

void save_image_ppm(const char * file, float * image, const int ResX, const int ResY) {
  FILE *out = fopen(file, "wb");
  if ( out ) {
    // convert to unsigned chars
    const unsigned int size = ResX * ResY * 3;
    unsigned char *buffer_c = new unsigned char[size];
    for ( int j = 0; j < ResY; j++ ) {
      for ( int i = 0; i < ResX * 3; i++ ) {
        float tmp = image[i + (ResY - j - 1) * ResX * 3] * 255.f;
        buffer_c[i + j * ResX * 3] = tmp < 255.f ? (unsigned char)tmp : 255;
      }
    }

    fprintf(out, "P6\n%d\n%d\n255\n", ResX, ResY);
    if ( fwrite(buffer_c, ResX * ResY * 3, 1, out) != 1 ) {
      std::cout << "Could not write image! " << std::endl;
    }
    fclose(out);
    delete [] buffer_c;
  }
}

bool load_image_ppm(const char * file, float *&image, int &ResX, int &ResY) {
  FILE *in = fopen(file, "rb");
  if ( in ) {
    if ( !fscanf(in, "P6\n%d\n%d\n255\n", &ResX, &ResY) ) return false;
    image = new float[ResX * ResY * 3];
    unsigned char *buffer_c = new unsigned char[ResX * ResY * 3];
    if ( fread(buffer_c, ResX * ResY * 3, 1, in) != 1 ) return false;
    for ( int j = 0; j < ResY; j++ ) {
      for ( int i = 0; i < ResX * 3; i++ ) {
        float tmp = (float)buffer_c[i + (ResY - j - 1) * ResX * 3] / 255.0f;
        image[i + j * ResX * 3] = tmp < 1.0f ? tmp : 1.0f;
      }
    }

    std::cout << ResX << " " << ResY << std::endl;
    fclose(in);
    delete buffer_c;
    return true;
  }
  return false;
}

void save_image_hdr(const char * file, float * image, const int ResX, const int ResY) {
  FILE *out_hdr = fopen(file, "wb");
  if ( out_hdr ) {
    float *buffer_f = new float[ResX * ResY * 3];
    for ( int j = 0; j < ResY; j++ )
      for ( int i = 0; i < ResX * 3; i++ )
        buffer_f[i + j * ResX * 3] = image[i + (ResY - j - 1) * ResX * 3];
    RGBE_WriteHeader(out_hdr, ResX, ResY, 0);
    RGBE_WritePixels_RLE(out_hdr, buffer_f, ResX, ResY);
    fclose(out_hdr);

    delete [] buffer_f;
  }
}

bool load_image_hdr(const char * file, float *&image, int &ResX, int &ResY) {
  FILE *in_hdr = fopen(file, "rb");
  if ( in_hdr ) {
    RGBE_ReadHeader(in_hdr, &ResX, &ResY, 0);
    image = new float[ResX * ResY * 3];
    float *buffer_f = new float[ResX * ResY * 3];
    RGBE_ReadPixels_RLE(in_hdr, buffer_f, ResX, ResY);
    for ( int j = 0; j < ResY; j++ )
      for ( int i = 0; i < ResX * 3; i++ )
        image[i + j * ResX * 3] = buffer_f[i + (ResY - j - 1) * ResX * 3];
    fclose(in_hdr);

    delete [] buffer_f;
    return true;
  }
  return false;
}


int load_float_data(const char * file, float * &data) {
  std::fstream f;
  f.open(file, std::ios::in | std::ios::binary);
  if ( !f ) return 0;

  f.seekg(0, std::ios::end);
  int length = f.tellg();
  f.seekg(0, std::ios::beg);

  data = new float[length/sizeof(float)];
  f.read((char*)data, length);

  f.close();
  return length / sizeof(float);
}

void save_data(const char * file, void * data, unsigned int num) {
  FILE *out = fopen(file, "wb");
  if ( out ) {
    if ( fwrite(data, 1, num, out) != num ) {
      std::cout << "Could not write data! " << std::endl;
    }
  }
}
