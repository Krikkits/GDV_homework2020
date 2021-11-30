#include <iostream>
#include "fileio.h"
#include "vec.h"
#include <stdio.h>
#include <stdlib.h>
#include <complex>

using namespace std;

/**
 * Squares a float.
 * @param v Input value.
 * @returns v*v
 */
inline float sq(float v)
{
	return v * v;
}

/**
 * Discrete Fourier Transformation.
 * @remarks Both in and out contain one element per pixel. Pixels are stored
 * in row major order.
 * @param out The fourier transformed input image is stored here.
 * @param in Input image.
 * @param resX Input/output image width in pixels.
 * @param resY Input/output image height in pixels.
 */
void DFT(complex<float>* out, complex<float>* in, unsigned int resX,
		unsigned int resY)
{
	// TODO 6.3 a) Implement the DFT according to the function's documentation.
	// To make OpenMP work, uncomment the following macros and place your most outer for-loop directly behind it.
	float dividedBy=1/sqrt(resX*resY);
	complex<float> complexnumber1, complexnumber2;
	complex<float> input;
#ifdef OPENMP
#pragma omp parallel for
#endif
	// Most outer for-loop here!


	for(unsigned int k1=0;k1<resX;k1++){
		//cout<<"k1: "<<k1<<endl;
		for(unsigned int k2=0; k2<resY;k2++){
			for(unsigned int n1=0;n1<resX;n1++){
				for(unsigned int n2=0;n2<resY;n2++){
					complexnumber1=-1i*2*M_PI*(k1/resX)*n1;
					complexnumber2=-1i*2*M_PI*(k2/resY)*n2;
					input=in[n2*resX+n1];
					out[k2*resX+k1]+= input * exp(complexnumber1)*exp(complexnumber2)*dividedBy;
					//cout<<"output: "<<out[k2*resX+k1]<<endl;
				}
			}

		}
	}


}

/**
 * Inverse Discrete Fourier Transformation.
 * @remarks Both in and out contain one element per pixel. Pixels are stored
 * in row major order.
 * @param out The restored image will be stored here.
 * @param in Input image containing complex fourier coefficients.
 * @param resX Input/output image width in pixels.
 * @param resY Input/output image height in pixels.
 */
void IDFT(complex<float>* out, complex<float>* in, unsigned int resX,
		unsigned int resY)
{
	// TODO 6.3 b) Implement the IDFT according to the function's documentation.
	// To make OpenMP work, uncomment the following macros and place your most outer for-loop directly behind it.
	float dividedBy=1/sqrt(resX*resY);
	complex<float> complexnumber1, complexnumber2;
	complex<float> input;
#ifdef OPENMP
#pragma omp parallel for
#endif
	// Most outer for-loop here!
	for(unsigned int n1=0;n1<resX;n1++){
		//cout<<"n1: "<<n1<<endl;
		for(unsigned int n2=0; n2<resY;n2++){
			for(unsigned int k1=0;k1<resX;k1++){
				for(unsigned int k2=0;k2<resY;k2++){
					complexnumber1=1i*2*M_PI*(k1/resX)*n1;
					complexnumber2=1i*2*M_PI*(k2/resY)*n2;
					input=in[k2*resX+k1];
					out[n2*resX+n1]+= input*exp(complexnumber1)*exp(complexnumber2)*dividedBy;
					//cout<<"output: "<<out[n2*resX+n1]<<endl;
				}
			}
		}
	}
}

/**
 * Converts a array of complex numbers to polar coordinate representation.
 * @param amplitudePhase The amplitude/phase representation of the data
 * contained in realImag will be stored here.
 * @param realImag Array of complex numbers which are converted to
 * polar coordinate (amplitude/phase) representation.
 * @param resX Input/output array width in pixels.
 * @param resY Input/output array height in pixels.
 */
void Complex2AmplitudePhase(Vec2* amplitudePhase, complex<float>* realImag,
		unsigned int resX, unsigned int resY)
{
	// TODO 6.3 c) Implement the conversion from Cartesian complex numbers (realImag[n]) to polar coordinate representation (amplitudePhase[n]).
	// amplitudePhase[n].x = amplitude; amplitudePhase[n].y = phase.

	for(unsigned int i=0;i<resX*resY;i++){

		float x=realImag[i].real();
		float y=realImag[i].imag();
		amplitudePhase[i].x=sqrt(sq(x)+sq(y)); //amplitude
		amplitudePhase[i].y=atan2f(y,x); //phase
		//cout<<"amplitude: "<<amplitudePhase[i].x<<endl;
	}
}

/**
 * Converts a array of polar coordinates (amplitude/phase) to cartesian
 * complex number representation.
 * @param realImag The cartesian representation of the data will be stored
 * here.
 * @param amplitudePhase Array containing the polar coordinate representation
 * of the data to be converted.
 * @param resX Input/output array width in pixels.
 * @param resY Input/output array height in pixels.
 */
void AmplitudePhase2Complex(complex<float>* realImag, Vec2* amplitudePhase,
		unsigned int resX, unsigned int resY)
{
	// TODO 6.3 d) Implement the conversion from polar coordinate representation (amplitudePhase[n]) to Cartesian complex numbers (realImag[n]).
	// amplitudePhase[n].x = amplitude; amplitudePhase[n].y = phase.
	for(unsigned int i=0;i<resX*resY;i++){

			float amplitude=amplitudePhase[i].x;
			float phase=amplitudePhase[i].y;
			realImag[i].real(cos(phase)*amplitude);
			realImag[i].imag(sin(phase)*amplitude);
		}

}

/**
 * Converts an array of real numbers to complex number representation.
 * @out The complex number representation (value from in + 0i) will be
 * stored here.
 * @in Array of real input numbers.
 * @param resX Input/output array width in pixels.
 * @param resY Input/output array height in pixels.
 */
void Real2Complex(complex<float>* out, float* in, unsigned int resX,
		unsigned int resY)
{
	for (unsigned int p = 0; p < resX * resY; p++)
	{
		out[p].real(in[p]);
		out[p].imag(0.0f);
	}
}

/**
 * Converts an array of complex numbers to real numbers by dropping the
 * imaginary component.
 * @out The resulting real numbers will be stored here.
 * @in Array of complex input numbers.
 * @param resX Input/output array width in pixels.
 * @param resY Input/output array height in pixels.
 */
void Complex2Real(float* out, complex<float>* in, unsigned int resX,
		unsigned int resY)
{
	for (unsigned int p = 0; p < resX * resY; p++)
		out[p] = in[p].real();
}

/**
 * Converts an array of 2d vectors to real numbers by dropping the
 * second component.
 * @out The resulting real numbers will be stored here.
 * @in Array of the 2d input vectors.
 * @param resX Input/output array width in pixels.
 * @param resY Input/output array height in pixels.
 */
void Vec2Scalar(float* out, Vec2* in, unsigned int resX, unsigned int resY)
{
	for (unsigned int p = 0; p < resX * resY; p++)
		out[p] = in[p].x;
}

/**
 * Converts an RGB image to grayscale.
 * @remarks Takes into account, how noticable the different color channels are
 * to humans.
 * @param out The grayscale pixel brightnesses will be stored here.
 * @param in Input RGB pixels.
 * @param resX Width of the input/output image.
 * @param resY Height of the input/output image.
 */
void rgb2Grayscale(float* out, Vec3* in, unsigned int resX, unsigned int resY)
{
	for (unsigned int p = 0; p < resX * resY; p++)
		out[p] = in[p].x * 0.3f + in[p].y * 0.59f + in[p].z * 0.11f;
}

/**
 * Converts a grayscale image to rgb representation.
 * @param out The rgb image pixels will be stored here.
 * @param in Input grayscale pixels.
 * @param resX Width of the input/output image.
 * @param resY Height of the input/output image.
 */
void grayscale2Rgb(Vec3* out, float* in, unsigned int resX, unsigned int resY)
{
	for (unsigned int p = 0; p < resX * resY; p++)
	{
		out[p].x = in[p];
		out[p].y = in[p];
		out[p].z = in[p];
	}
}

/**
 * Normalizes a grayscale image by dividing all pixels by the value of the
 * brightest pixel.
 * @params image Input/output image.
 * @params resX Width of image.
 * @params resY Height of image.
 */
void normalize(float* image, unsigned int resX, unsigned int resY)
{
	float max = 0.0f;
	for (unsigned int p = 0; p < resX * resY; p++)
		max = image[p] > max ? image[p] : max;
	if (max == 0.0f)
		return;
	for (unsigned int p = 0; p < resX * resY; p++)
		image[p] /= max;
}

/**
 * Transforms the pixel values of an image to a logarithmic scale
 * to make dark areas more visible.
 * @params image Input/output image.
 * @params resX Width of image.
 * @params resY Height of image.
 */
void log(float* image, unsigned int resX, unsigned int resY)
{
	for (unsigned int p = 0; p < resX * resY; p++)
	{
		image[p] = log(image[p]) * 0.1f + 1.0f;
		if (image[p] < 0.0f)
			image[p] = 0.0f;
	}
}

/**
 * Rotates an image by half its size in vertical and horizontal direction.
 * @remarks Can be used to shift the lowest frequency fourier coefficients
 * to the center of the image and back.
 * @param out The shifted image will be stored here.
 * @param in Input image.
 * @param resX Width of the input/output image.
 * @param resY Height of the input/output image.
 */
void shiftHalf(float* out, float* in, unsigned int resX, unsigned int resY)
{
	for (unsigned int yi = 0; yi < resY; yi++)
	{
		for (unsigned int xi = 0; xi < resX; xi++)
		{
			unsigned int xo = (xi + (resX / 2)) % resX;
			unsigned int yo = (yi + (resY / 2)) % resY;

			out[yo * resX + xo] = in[yi * resX + xi];
		}
	}
}

/**
 * Main function. Loads the image whose filename is given as parameter
 * and processes it:
 * <ul>
 *  <li>Conversion to grayscale *</li>
 *  <li>Conversion to complex numbers</li>
 *  <li>Conversion to fourier space</li>
 *  <li>Conversion to polar coordinate (amplitude/phase) representation *</li>
 *  <li>Back-conversion to cartesian fourier space</li>
 *  <li>Back-conversion to image space</li>
 *  <li>Back-conversion to grayscale (non-complex)*</li>
 * </ul>
 * An image is saved for each step marked with a *.
 * @param argc Number of program parameters.
 * @param argv Array of program parameters.
 */
int main(int argc, char **argv)
{
	if (argc < 2)
		return 0;

	float*floatRgbImage;

	int resX, resY;
	load_image_ppm(argv[1], floatRgbImage, resX, resY);

	float* tmpImage1 = new float[resX * resY];
	float* tmpImage1_ = new float[resX * resY];
	complex<float>* tmpImage2 = new complex<float> [resX * resY];
	Vec3* tmpImage3 = new Vec3[resX * resY];

	float* grayImage = new float[resX * resY];
	complex<float>* gray2complexImage = new complex<float> [resX * resY];
	complex<float>* gray2complex2fourierImage = new complex<float> [resX * resY];
	Vec2* gray2complex2fourier2ampPhaseImage = new Vec2[resX * resY];

	complex<float>* gray2complex2fourier2ampPhase2fourierImage = new complex<
			float> [resX * resY];
	complex<float>* gray2complex2fourier2ampPhase2fourier2complexImage =
			new complex<float> [resX * resY];
	float* gray2complex2fourier2ampPhase2fourier2complex2grayImage =
			new float[resX * resY];

	// Convert to grayscale
	rgb2Grayscale(grayImage, (Vec3*) floatRgbImage, resX, resY);
	// Output
	grayscale2Rgb(tmpImage3, grayImage, resX, resY);
	save_image_ppm("grayImage.ppm", (float*) tmpImage3, resX, resY);

	// Convert to fourier space
	Real2Complex(gray2complexImage, grayImage, resX, resY);
	DFT(gray2complex2fourierImage, gray2complexImage, resX, resY);

	// Convert to amplitude/phase description
	Complex2AmplitudePhase(gray2complex2fourier2ampPhaseImage,
			gray2complex2fourierImage, resX, resY);
	// Output
	Vec2Scalar(tmpImage1, gray2complex2fourier2ampPhaseImage, resX, resY);
	shiftHalf(tmpImage1_, tmpImage1, resX, resY);
	normalize(tmpImage1_, resX, resY);
	log(tmpImage1_, resX, resY);
	grayscale2Rgb(tmpImage3, tmpImage1_, resX, resY);
	save_image_ppm("gray2complex2fourier2ampPhaseImage.ppm", (float*) tmpImage3,
			resX, resY);

	// Convert back to real/imaginary vector fourier base description
	AmplitudePhase2Complex(gray2complex2fourier2ampPhase2fourierImage,
			gray2complex2fourier2ampPhaseImage, resX, resY);

	// Convert back to grayscale
	IDFT(gray2complex2fourier2ampPhase2fourier2complexImage,
			gray2complex2fourier2ampPhase2fourierImage, resX, resY);
	Complex2Real(gray2complex2fourier2ampPhase2fourier2complex2grayImage,
			gray2complex2fourier2ampPhase2fourier2complexImage, resX, resY);
	// Output
	grayscale2Rgb(tmpImage3,
			gray2complex2fourier2ampPhase2fourier2complex2grayImage, resX,
			resY);
	save_image_ppm(
			"gray2complex2fourier2ampPhase2fourier2complex2grayImage.ppm",
			(float*) tmpImage3, resX, resY);

	delete[] floatRgbImage;

	delete[] tmpImage1;
	delete[] tmpImage1_;
	delete[] tmpImage2;
	delete[] tmpImage3;

	delete[] grayImage;
	delete[] gray2complexImage;
	delete[] gray2complex2fourierImage;
	delete[] gray2complex2fourier2ampPhaseImage;

	delete[] gray2complex2fourier2ampPhase2fourierImage;
	delete[] gray2complex2fourier2ampPhase2fourier2complexImage;
	delete[] gray2complex2fourier2ampPhase2fourier2complex2grayImage;

	return 0;
}

