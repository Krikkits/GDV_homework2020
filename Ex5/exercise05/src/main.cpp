/**
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#include <iostream>
#include <fstream>
#include <omp.h>
#include "cam.h"
#include "scene.h"
#include "render.h"
#include "utils/vec.h"
#include "utils/fileio.h"
#include "utils/MersenneTwister.h"

/// Renderer object
Render *render;
/// Camera object
Cam *cam;

/// Currently selected shader
int shader = 0;

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
				cam->cam_w(true);
				break;
			case SDLK_a:
				cam->cam_a(true);
				break;
			case SDLK_s:
				cam->cam_s(true);
				break;
			case SDLK_d:
				cam->cam_d(true);
				break;
			case SDLK_1:
				shader = 1;
				render->accum_index = 0;
				break;
			case SDLK_2:
				shader = 2;
				render->accum_index = 0;
				break;
			case SDLK_3:
				shader = 3;
				render->accum_index = 0;
				break;
			case SDLK_4:
				shader = 4;
				render->accum_index = 0;
				break;
			case SDLK_5:
				shader = 5;
				render->accum_index = 0;
				break;
			case SDLK_6:
				shader = 6;
				render->accum_index = 0;
				break;
			default:
				break;
			}
			break;
		case SDL_KEYUP:
			switch (event.key.keysym.sym)
			{
			case SDLK_w:
				cam->cam_w(false);
				break;
			case SDLK_a:
				cam->cam_a(false);
				break;
			case SDLK_s:
				cam->cam_s(false);
				break;
			case SDLK_d:
				cam->cam_d(false);
				break;
			default:
				break;
			}
			break;
		case SDL_MOUSEMOTION:
			if (SDL_GetMouseState(NULL, NULL) & SDL_BUTTON(1))
			{
				cam->cam_rx(event.motion.yrel);
				cam->cam_ry(event.motion.xrel);
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
int ResX;
/// Height of the rendered image.
int ResY;

/**
 * Main program routine.
 */
int main(int argc, char **argv)
{
	const char* sceneFile = "CornellBox";
	const char* envFile = 0;

	if (argc >= 2)
		sceneFile = argv[1];
	if (argc >= 3)
		envFile = argv[2];

	Scene *scene = new Scene(sceneFile, envFile);

	cam = scene->cam;
	ResX = cam->ResX;
	ResY = cam->ResY;

	render = new Render(scene);

	// just for illustration how the MersenneTwister pseudo random number generator is used
	MTRand *drand = new MTRand(1337); // the initialisation is arbitrary, but never initialize randomly (i.e. time or /dev/ramdom), this would make
									  // debugging much more complicated (i.e. during rendering a certain sequence of random numbers generates
									  // a ray which results in a segmentation fault during ray traversal -> this can not be debugged with
									  // a randomly initialized PRNG.
	std::cout << "Pseudo random number generation example, see main.cpp"
			<< std::endl << std::endl;
	std::cout << "10 real pseudo random numbers in [0,1]: " << std::endl;
	for (int i = 0; i < 10; i++)
	{
		std::cout << drand->rand() << " ";
	}
	std::cout << std::endl;
	std::cout << "10 integer pseudo random numbers in [0,2^32-1]: "
			<< std::endl;
	for (int i = 0; i < 10; i++)
	{
		std::cout << drand->randInt() << " ";
	}
	std::cout << std::endl;
	std::cout << "10 integer pseudo random numbers in [0,777]: " << std::endl;
	for (int i = 0; i < 10; i++)
	{
		std::cout << drand->randInt(777) << " ";
	}
	std::cout << std::endl;

#ifdef INTERACTIVE
	initScreen(ResX, ResY);
	char title[256];
	float fps = 0.0f;
	int frame = 0;
	unsigned int t0 = SDL_GetTicks();

	while (!finished)
	{
		eventHandling();
		cam->cam_move();

		render->render(shader);

		glDrawPixels(ResX, ResY, GL_RGB, GL_FLOAT, (float*) render->image);
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
	render->render(shader);
#endif
	save_image_ppm("image.ppm", (float*) render->image, ResX, ResY);

	return 0;
}

