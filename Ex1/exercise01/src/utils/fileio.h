#ifndef FILEIO_H
#define FILEIO_H

void save_image_ppm(const char * file, float * image, const int ResX, const int ResY);
bool load_image_ppm(const char * file, float *&image, int &ResX, int &ResY);

int load_float_data(const char * file, float * &data);

#endif

