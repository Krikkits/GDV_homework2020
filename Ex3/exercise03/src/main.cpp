/**
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#include <iostream>
#include <fstream>
#include <omp.h>
#include "cam.h"
#include "rtStructs.h"
#include "bvh.h"
#include "utils/vec.h"
#include "utils/fileio.h"
#include "utils/MersenneTwister.h"

Cam *cam;

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

/// Random number generator.
MTRand *drand;

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

/// For pointlights to light the scene
Pointlight pointlights[4];
/// Area light used to light the scene
AreaLight areaLight;
/// Index of the first triangle that should be rendered as mirroring instead of diffuse.
int firstMirrorTriangle;

// Uncomment exactly one of the following lines!
//#define FLAT_SHADER
//#define DEBUG_SHADER
//#define SMOOTH_SHADER
#define POINTLIGHT_SHADER
//#define AREALIGHT_SHADER

/**
 * Determines the color seen along a certain ray.
 * @param ray This ray is used for RayTracing.
 * @returns Color seen along the given ray.
 */
Vec3 rayTrace(const Ray &ray)
{
	HitRec rec = bvh->intersect(ray);

	if (rec.id != -1)
	{
#ifdef FLAT_SHADER
		Vec3 normal = bvh->tris[rec.id].getNormal();
		return fabsf(ray.dir * normal);
#endif

#ifdef DEBUG_SHADER
		// TODO 3.3 d) Implement a debug shader here.
		Vec3 normal = bvh->tris[rec.id].getInterpolatedNormal(ray);
		Vec3 RGBs;
		RGBs[0]=fabsf(normal[0]);
		RGBs[1]=fabsf(normal[1]);
		RGBs[2]=fabsf(normal[2]);
		return RGBs;// Replace with debug shader shade color.
#endif

#ifdef SMOOTH_SHADER
		// TODO 3.3 e) Implement diffuse shading using the interpolated vertex normal here.
		Vec3 direction = ray.dir;
		Vec3 dirnorm = ray.dir * direction.normalize();
		Vec3 hitpoint = ray.origin + rec.dist * dirnorm;

		//diffuse function from first exercise
		Vec3 diffuse=(bvh->tris[rec.id].getInterpolatedNormal(ray)*hitpoint)/(bvh->tris[rec.id].getInterpolatedNormal(ray).length()*hitpoint.length());
		return diffuse;// Replace with smooth shader shade color.
#endif

#ifdef POINTLIGHT_SHADER
		// TODO 3.4 a) Calculate the shade by summing up the contributions of all four pointlights.
		// Take the direction towards the light source, the surface normal and the distance into account.

		//calculates the hitpoint on the bunny
		Vec3 direction = ray.dir;
		Vec3 dirnorm = ray.dir * direction.normalize();
		Vec3 hitpoint = ray.origin + rec.dist * direction;
		//Vec3 lightDir=hitpoint-pointlights->pos;
		Vec3 interpolNormal = bvh->tris[rec.id].getInterpolatedNormal(ray);
		Vec3 interpolNormalNorm=interpolNormal.normalize();

		//calculates the vector between hitpoint and  each pointlight
		Vec3 ptol[4];

		//distance between hitpoint and pointlight
		float disttol[4];

		Vec3 colour; //our return

		//normailzed direction from hitpoint to pointlight
		Vec3 normptol[4];

		for (int i = 0; i < 4; i++){
			ptol[i] = pointlights[i].pos - hitpoint;
			disttol[i] = ptol[i].normalize();
			normptol[i] = ptol[i]*disttol[i];
			Vec3 interpolColour=Vec3(pointlights[i].color.x * interpolNormal.x, pointlights[i].color.y * interpolNormal.y, pointlights[i].color.z * interpolNormal.z);
			colour += interpolColour * (1/(pow(disttol[i],2)));
		}




		// TODO 3.4 b) Check if each lightsource is visible from the shaded surface at all.
		// Only add its contribution if it is.
		//Rays from hitpoint to pointlights
		Ray rays[4];
		for (int i=0; i<4; i++){
			rays[i].dir = ptol[i];
			rays[i].origin = hitpoint;
			rays[i].tmin = RAY_EPS;
			rays[i].tmax = RAY_EPS;
			//my attempt to make shadows
			if (bvh->tris[rec.id].intersect(rays[i], rec, rec.id)){
				colour = colour-(pointlights[i].color * 1/(pow(disttol[i],2)));
			}
		}


		// TODO 3.5 Shade all triangles with rec.id >= firstMirrorTriangle with mirror shading instead of pointlight shading.
		 if (rec.id >= firstMirrorTriangle){
			 Vec3 reflectionDir;
			 reflectionDir = ray.dir - 2 * (ray.dir * bvh->tris[rec.id].getInterpolatedNormal(ray)) * bvh->tris[rec.id].getInterpolatedNormal(ray);
			 Ray reflection(hitpoint,reflectionDir, RAY_EPS + ray.tmin,rec.dist);
			 colour = rayTrace(reflection);
			 return colour;
		 }

		return colour;// Replace with point light / mirror shader shade color.
#endif

#ifdef AREALIGHT_SHADER
		// TODO 3.6 b) Implement area light shading by picking a random point light on the area light's surface.
		//calculates the hitpoint on the bunny
		Vec3 direction = ray.dir;
		Vec3 dirnorm = ray.dir * direction.normalize();
		Vec3 randpoint = ray.origin + rec.dist * drand->rand();
		Vec3 neighbor1=areaLight.pos[0]+areaLight.extent1;
		Vec3 neighbor2=areaLight.pos[1]+areaLight.extent2;
		Vec3 colour; //our return

		/* theoretically:
		Vec3 newLocation;
		newLocation[0] = pos[0] + area * drand->randInt(1000);
		newLocation[1] = pos[1] + area * drand->randInt(1000);
		newLocation[2] = pos[2] + area * drand->randInt(1000);
		*/

		// TODO 3.6 c) Implement shading by averaging the shade from 1000 random point lights from the area light's surface.
		colour = areaLight.radiance * areaLight.getArea();

		return colour;// Replace with area light shader shade color.
#endif
	}

	return Vec3(0.0f, 0.0f, 0.0f);
}

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

			buffer[x + y * ResX] = rayTrace(ray);
		}
	}
}

/**
 * Main program routine.
 */
int main(int argc, char **argv)
{

	//load bunny.n
	float *nData;
	int n_data;

	const char* nFilename = argc >= 2 ? argv[1] : "bunny.n";
	n_data = load_float_data(nFilename, nData);



	float *meshData;
	int num_data;

	const char* meshFilename = argc >= 2 ? argv[1] : "bunny.ra2";
	num_data = load_float_data(meshFilename, meshData);

	// The two additional triangles are for the ground...
	int num_tris = num_data / 9 + 2;
	Triangle *tris = new Triangle[num_tris];

	for (int t = 0; t < num_tris; t++)
	{
		for (unsigned int v = 0; v < 3; v++)
		{
			for (int d = 0; d < 3; d++)
			{
				tris[t].v[v][d] = meshData[t * 9 + v * 3 + d];
				tris[t].n[v][d] = nData[t * 9 + v * 3 + d];
			}
		}
	}
	delete[] meshData;
	delete[] nData;



	// Bounding box around the loaded part of the scene.
	AABB sceneBox = Triangle::getAABB(tris, num_tris - 2);

	// Initialization of the point lights.
	pointlights[0].pos = Vec3(sceneBox.bounds[0].x, sceneBox.bounds[1].y,
			sceneBox.bounds[0].z);
	pointlights[0].color = Vec3(1500.0f, 1500.0f, 1500.0f);
	pointlights[1].pos = Vec3(sceneBox.bounds[1].x, sceneBox.bounds[1].y,
			sceneBox.bounds[0].z);
	pointlights[1].color = Vec3(1500.0f, 0.0f, 0.0f);
	pointlights[2].pos = Vec3(sceneBox.bounds[0].x, sceneBox.bounds[1].y,
			sceneBox.bounds[1].z);
	pointlights[2].color = Vec3(0.0f, 1500.0f, 0.0f);
	pointlights[3].pos = Vec3(sceneBox.bounds[1].x, sceneBox.bounds[1].y,
			sceneBox.bounds[1].z);
	pointlights[3].color = Vec3(0.0f, 0.0f, 1500.0f);

	// Initialization of the area light.
	areaLight.pos = Vec3(sceneBox.bounds[0].x, sceneBox.bounds[1].y,
			sceneBox.bounds[0].z);
	areaLight.extent1 = Vec3(sceneBox.bounds[1].x - sceneBox.bounds[0].x, 0.0f,
			0.0f);
	areaLight.extent2 = Vec3(0.0f, 0.0f,
			sceneBox.bounds[1].z - sceneBox.bounds[0].z);
	areaLight.radiance = Vec3(1.0f, 1.0f, 1.0f);

	// Adding hardcoded ground plane which should be rendered as mirror in 3.5
	float boxWidth = sceneBox.bounds[1].x - sceneBox.bounds[0].x;
	float boxHeight = sceneBox.bounds[1].z - sceneBox.bounds[0].z;
	sceneBox.bounds[0].x -= boxWidth;
	sceneBox.bounds[1].x += boxWidth;
	sceneBox.bounds[0].z -= boxHeight;
	sceneBox.bounds[1].z += boxHeight;
	tris[num_tris - 2].v[0] = Vec3(sceneBox.bounds[0].x, sceneBox.bounds[0].y,
			sceneBox.bounds[1].z);
	tris[num_tris - 2].v[1] = Vec3(sceneBox.bounds[1].x, sceneBox.bounds[0].y,
			sceneBox.bounds[1].z);
	tris[num_tris - 2].v[2] = Vec3(sceneBox.bounds[0].x, sceneBox.bounds[0].y,
			sceneBox.bounds[0].z);
	tris[num_tris - 1].v[0] = Vec3(sceneBox.bounds[1].x, sceneBox.bounds[0].y,
			sceneBox.bounds[1].z);
	tris[num_tris - 1].v[1] = Vec3(sceneBox.bounds[1].x, sceneBox.bounds[0].y,
			sceneBox.bounds[0].z);
	tris[num_tris - 1].v[2] = Vec3(sceneBox.bounds[0].x, sceneBox.bounds[0].y,
			sceneBox.bounds[0].z);

	// TODO 3.3 a) Initialize the normals of the ground triangles with (0, 1, 0).
	// If your normals in "Vec3 Triangle::n", you may just uncomment the following lines.

	 tris[num_tris - 2].n[0] = Vec3(0.0f, 1.0f, 0.0f);
	 tris[num_tris - 2].n[1] = Vec3(0.0f, 1.0f, 0.0f);
	 tris[num_tris - 2].n[2] = Vec3(0.0f, 1.0f, 0.0f);
	 tris[num_tris - 1].n[0] = Vec3(0.0f, 1.0f, 0.0f);
	 tris[num_tris - 1].n[1] = Vec3(0.0f, 1.0f, 0.0f);
	 tris[num_tris - 1].n[2] = Vec3(0.0f, 1.0f, 0.0f);

	firstMirrorTriangle = num_tris - 2;

	std::cout << "#Triangles " << num_tris << std::endl;

	bvh = new BVH(tris, num_tris);

	std::cout << bvh->bbox.bounds[0][0] << " " << bvh->bbox.bounds[0][1] << " "
			<< bvh->bbox.bounds[0][2] << std::endl;
	std::cout << bvh->bbox.bounds[1][0] << " " << bvh->bbox.bounds[1][1] << " "
			<< bvh->bbox.bounds[1][2] << std::endl << std::endl;

	buffer = new Vec3[ResX * ResY];

	cam = new Cam(bvh->bbox, ResX, ResY);

	drand = new MTRand(1337); // the initialisation is arbitrary, but never initialize randomly (i.e. time or /dev/ramdom), this would make
							  // debugging much more complicated (i.e. during rendering a scertain sequence of random numbers generates
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

