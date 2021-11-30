#include <iostream>
#include <fstream>
#include <omp.h>
#include "cam.h"
#include "rtStructs.h"
#include "bvh.h"
#include "utils/vec.h"
#include "utils/fileio.h"

PerspCam *cam;

#ifdef INTERACTIVE
#include <SDL/SDL.h>
#include <SDL_opengl.h>

bool finished = false;

SDL_Surface *screen;
void initScreen(int ResX, int ResY)
{
	SDL_Init(SDL_INIT_VIDEO);
	SDL_SetVideoMode(ResX, ResY, 32, SDL_OPENGL | SDL_GL_DOUBLEBUFFER);
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
			case SDLK_w:
				break;
			case SDLK_a:
				break;
			case SDLK_s:
				break;
			case SDLK_d:
				break;
			default:
				break;
			}
			break;
		case SDL_KEYUP:
			switch (event.key.keysym.sym)
			{
			case SDLK_w:
				break;
			case SDLK_a:
				break;
			case SDLK_s:
				break;
			case SDLK_d:
				break;
			default:
				break;
			}
			break;
		case SDL_MOUSEMOTION:
			if (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(1))
			{
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

/// Width of the rendered image.
int ResX = 300;
/// Height of the rendered image.
int ResY = 300;

/**
 * Buffer containing the rendered images.
 * @remarks Pixels are stored rowwise, so pixel positions can be calculated by
 * y*ResX+x.
 */
Vec3 *buffer;
/// BVH tree used for accelerating the rendering.
BVH *bvh;

//Triangle *tris; data is now in the bvh structure
//int num_tris;

/**
 * Renders an image of ResX * ResY size into buffer.
 */
void render()
{
	for (int y = 0; y < ResY; y++)
	{
		for (int x = 0; x < ResX; x++)
		{
			Ray ray = cam->getRay(x, y);

			HitRec rec = bvh->intersect(ray);

			if (rec.id != -1)
			{
				Vec3 normal = bvh->tris[rec.id].getNormal();
				buffer[x + y * ResX] = fabsf(ray.dir * normal);
			}
			else
				buffer[x + y * ResX] = Vec3(0.0f, 0.0f, 0.0f);
		}
	}
}

int main(int argc, char **argv)
{
	float *data;
	int num_data;

	if (argc == 2)
		num_data = load_float_data(argv[1], data);
	else
		num_data = load_float_data("Test.ra2", data);

	Triangle *tris = (Triangle*) data;
	int num_tris = num_data / 9;

	std::cout << "#Triangles " << num_tris << std::endl;

	bvh = new BVH(tris, num_tris);

	buffer = new Vec3[ResX * ResY];

	cam = new PerspCam(Vec3(0.0f, 0.0f, 120.0f), Vec3(0.0f, 0.0f, -1.0f),
			Vec3(0.0f, 1.0f, 0.0f), ResX, ResX, ResY);

#ifdef INTERACTIVE
	initScreen(ResX, ResY);
	char title[256];
	float fps = 0.0f;
	int frame = 0;
	unsigned int t0 = SDL_GetTicks();

	while (!finished)
	{
		eventHandling();

		render();

		glDrawPixels(ResX, ResY, GL_RGB, GL_FLOAT, (float*) buffer);
		SDL_GL_SwapBuffers();
		sprintf(title, "%g fps", fps);
		SDL_WM_SetCaption(title, NULL);
		unsigned int t1 = SDL_GetTicks();
		frame++;
		if (t1 - t0 > 500.0f)
		{
			fps = 1000.0f / (t1 - t0) * frame;
			t0 = t1;
			frame = 0;
		}
	}
#else
	render();
#endif
	save_image_ppm("image.ppm", (float*) buffer, ResX, ResY);

	return 0;
}

