#include <iostream>
#include <fstream>
#include <vector>
#include <omp.h>
#include "cam.h"
#include "rtStructs.h"
#include "bvh.h"
#include "utils/vec.h"
#include "utils/fileio.h"
#include "utils/MersenneTwister.h"

using namespace std;

PerspCam *cam;

#ifdef INTERACTIVE
#include <SDL/SDL.h>
#include <SDL_opengl.h>

bool finished = false;
bool switched = true;
int SamplingMode = 1;

SDL_Surface *screen;
void initScreen(int ResX, int ResY)
{
	SDL_Init(SDL_INIT_VIDEO);
	SDL_SetVideoMode(ResX * 3, ResY * 3, 32, SDL_OPENGL | SDL_GL_DOUBLEBUFFER);
}

void eventHandling()
{
	SDL_Event event;
	while (SDL_PollEvent(&event))
	{
		switch (event.type)
		{
		case SDL_KEYDOWN:
			switch (event.key.keysym.sym)
			{
			case SDLK_ESCAPE:
				finished = true;
				break;
			case SDLK_q:
				finished = true;
				break;
			case SDLK_1:
				SamplingMode = 1;
				switched = true;
				break;
			case SDLK_2:
				SamplingMode = 2;
				switched = true;
				break;
			case SDLK_3:
				SamplingMode = 3;
				switched = true;
				break;
			case SDLK_4:
				SamplingMode = 4;
				switched = true;
				break;
			case SDLK_5:
				SamplingMode = 5;
				switched = true;
				break;
			case SDLK_6:
				SamplingMode = 6;
				switched = true;
				break;
			case SDLK_7:
				SamplingMode = 7;
				switched = true;
				break;
			default:
				break;
			}
			break;
		case SDL_QUIT:
			finished = true;
			break;
		default:
			break;
		}
	}
}

#endif

MTRand* mtrand;

/// Width of the rendered image.
int ResX = 200;
/// Height of the rendered image.
int ResY = 200;

/**
 * Buffer containing the rendered images.
 * @remarks Pixels are stored rowwise, so pixel positions can be calculated by
 * y*ResX+x.
 */
Vec3 *buffer;

/// Stores an image of the sampling pattern of one pixel.
Vec3 *patternBuffer;
/// BVH tree used for accelerating the rendering.
BVH *bvh;

const float diskRadius = 0.25f;

//Triangle *tris; data is now in the bvh structure
//int num_tris;

Vec3 rayTrace(const Ray &ray)
{
	HitRec rec = bvh->intersect(ray);

	if (rec.id != -1)
	{
		Vec3 normal = bvh->tris[rec.id].getNormal();
		return fabsf(ray.dir * normal);
	}

	return Vec3(0.0f, 0.0f, 0.0f);
}

/**
 * Renders one pixel with multiple samples.
 * @param x X coordinate of the pixel.
 * @param y Y coordinate of the pixel.
 * @param offsets Contains one element per sample that should be rendered.
 * Sampling offsets are relative to the bottom left corner of the pixel and assume
 * a pixel size of 1x1 unit. (x=right, y=up)
 * @returns Pixel color averaged from multiple samples.
 */
Vec3 renderPixel(int x, int y, const vector<Vec2> &offsets)
{
	Vec3 color(0.0f);
	for (unsigned int r = 0; r < offsets.size(); r++)
	{
		Ray ray = cam->getRay((float) x + offsets[r].x,
				(float) y + offsets[r].y);
		color += rayTrace(ray);
	}
	if (offsets.size() > 0)
		color *= (1.0f / (float) offsets.size());

	return color;
}

/**
 * Renders a small debugging image which shows the distribution of the samples.
 * @param offsets Contains one element per sample that should be rendered.
 * Sampling offsets are relative to the bottom left corner of the pixel and assume
 * a pixel size of 1x1 unit. (x=right, y=up)
 */
void renderPatternBuffer(const vector<Vec2> &offsets)
{
	for (unsigned int p = 0; p < 100 * 100; p++)
		patternBuffer[p] = Vec3(1.0f);

	for (unsigned int s = 0; s < offsets.size(); s++)
	{
		unsigned int x = (int) ((float) offsets[s].x * 99.9f);
		unsigned int y = (int) ((float) offsets[s].y * 99.9f);
		if (x > 99 || y > 99) // Avoid segfaults...
			continue;
		patternBuffer[y * 100 + x] = Vec3(1.0f, 0.0f, 0.0f);
	}
}

/**
 * Renders an image of ResX * ResY size into buffer.
 */
void render(int sampling)
{
	for (int y = 0; y < ResY; y++)
	{
		for (int x = 0; x < ResX; x++)
		{
			vector<Vec2> samples;
			Vec2 randSamples[4];
			Vec2 nextrandSample;
			float diskR=diskRadius;

			switch (sampling)
			{

			case 1: // No supersampling
				samples.push_back(Vec2(0.5f, 0.5f));
				break;
			case 2: // Regular grid

				// TODO 6.4 a) Push back 4 regular grid samples to samples
				samples.push_back(Vec2(0.25, 0.25));
				samples.push_back(Vec2(0.25, 0.75));
				samples.push_back(Vec2(0.75, 0.25));
				samples.push_back(Vec2(0.75, 0.75));
				break;
			case 3: // Rotated grid

				// TODO 6.4 b) Push back 4 rotated grid samples to samples.

				//using values according to wikipedia
				//round up 26.6 degrees
				samples.push_back(Vec2((0.25 * cos(atan2f(2,1)) - 0.25 * sin(atan2f(2,1))), (0.25 * sin(atan2f(2,1)) + 0.25 * cos(atan2f(2,1))))*sqrt(5)/2);
				samples.push_back(Vec2((0.25 * cos(atan2f(2,1)) - 0.75 * sin(atan2f(2,1))), (0.25 * sin(atan2f(2,1)) + 0.75 * cos(atan2f(2,1))))*sqrt(5)/2);
				samples.push_back(Vec2((0.75 * cos(atan2f(2,1)) - 0.25 * sin(atan2f(2,1))), (0.75 * sin(atan2f(2,1)) + 0.25 * cos(atan2f(2,1))))*sqrt(5)/2);
				samples.push_back(Vec2((0.75 * cos(atan2f(2,1)) - 0.75 * sin(atan2f(2,1))), (0.75 * sin(atan2f(2,1)) + 0.75 * cos(atan2f(2,1))))*sqrt(5)/2);

				break;
			case 4: // Random

				// TODO 6.4 c) Push back 4 random samples to samples.
				samples.push_back(Vec2((mtrand->rand()), (mtrand->rand())));
				samples.push_back(Vec2((mtrand->rand()), (mtrand->rand())));
				samples.push_back(Vec2((mtrand->rand()), (mtrand->rand())));
				samples.push_back(Vec2((mtrand->rand()), (mtrand->rand())));

				break;
			case 5: // Poisson disk

				// TODO 6.4 d) Push back 4 samples which have at least 2*diskRadius distance to each other.
				randSamples[0]=(mtrand->rand()), (mtrand->rand());
				samples.push_back(randSamples[0]);

				for(int i=0;i<4;i++){
					nextrandSample=(mtrand->rand()), (mtrand->rand());
					if((nextrandSample-randSamples[i]).normalize()<=2*diskR){
						nextrandSample=(mtrand->rand()), (mtrand->rand());
					}else{
						samples.push_back(nextrandSample);
					}
				}


				break;
			case 6: // Jittered

				// TODO 6.4 e) Push back 4 samples, which are jittered up to 0.2 pixels around the regular grid sampling positions.
				samples.push_back(Vec2(0.25*mtrand->rand(), 0.25*mtrand->rand()));
				samples.push_back(Vec2(0.25*mtrand->rand(), 0.75*mtrand->rand()));
				samples.push_back(Vec2(0.75*mtrand->rand(), 0.25*mtrand->rand()));
				samples.push_back(Vec2(0.75*mtrand->rand(), 0.75*mtrand->rand()));
				break;
			case 7: // Stratified

				// TODO 6.4 f) Push back 4 stratified samples.
				samples.push_back(Vec2(0.5f, 0.5f));
				break;
			}

			buffer[y * ResX + x] = renderPixel(x, y, samples);

			if (x == 0 && y == 0)
				renderPatternBuffer(samples);
		}
	}
}

int main(int argc, char **argv)
{
	float *data;
	int num_data;

	mtrand = new MTRand(1234);

	if (argc == 2)
		num_data = load_float_data(argv[1], data);
	else
		num_data = load_float_data("Test.ra2", data);

	Triangle *tris = (Triangle*) data;
	int num_tris = num_data / 9;

	std::cout << "#Triangles " << num_tris << std::endl;

	bvh = new BVH(tris, num_tris);

	std::cout << bvh->bbox.bounds[0][0] << " " << bvh->bbox.bounds[0][1] << " "
			<< bvh->bbox.bounds[0][2] << std::endl;
	std::cout << bvh->bbox.bounds[1][0] << " " << bvh->bbox.bounds[1][1] << " "
			<< bvh->bbox.bounds[1][2] << std::endl;

	buffer = new Vec3[ResX * ResY];
	patternBuffer = new Vec3[100 * 100];

	cam = new PerspCam(Vec3(0.0f, 0.0f, -1000.0f), Vec3(0.0f, 0.0f, 1.0f),
			Vec3(0.0f, 1.0f, 0.0f), 45, ResX, ResY);

#ifdef INTERACTIVE
	initScreen(ResX, ResY);

	while (!finished)
	{
		do
			eventHandling();
		while (!switched && !finished);

		render(SamplingMode);

		glPixelZoom(3.0f, 3.0f);
		glDrawPixels(ResX, ResY, GL_RGB, GL_FLOAT, (float*) buffer);
		glPixelZoom(1.0f, 1.0f);
		glDrawPixels(100, 100, GL_RGB, GL_FLOAT, (float*) patternBuffer);
		SDL_GL_SwapBuffers();
		switched = false;
	}
#else
	render();
#endif
	save_image_ppm("image.ppm", (float*) buffer, ResX, ResY);

	delete[] patternBuffer;
	delete[] buffer;
	delete mtrand;

	return 0;
}

