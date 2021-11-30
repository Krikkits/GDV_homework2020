/**
 * Computer Graphics Lecture WS 2009/2010 Ulm University
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#include "render.h"

#include <omp.h>

Render::Render(Scene *scene) :
		scene(scene)
{
	accel = new BVH(scene->triangles, scene->num_tris);

	cam = scene->cam;
	ResX = cam->ResX;
	ResY = cam->ResY;

	image = new Vec3[ResX * ResY];
	for (int i = 0; i < ResX * ResY; i++)
		image[i] = Vec3(0.0f, 0.0f, 0.0f);

#ifndef OPENMP
	mtrand = new MTRand*[0];
	mtrand[0] = new MTRand(1337);
#else
	mtrand = new MTRand*[omp_get_max_threads()];
	for (int t = 0; t < omp_get_max_threads(); t++)
		mtrand[t] = new MTRand(1337 + t);
#endif
	accum_index = 0;
}

Render::~Render()
{
	delete accel;
	delete[] image;
#ifndef OPENMP
	delete mtrand[0];
#else
	for (int t = 0; t < omp_get_max_threads(); t++)
		delete mtrand[t];
#endif
	delete[] mtrand;
}

void Render::render(int shader)
{
	float inv_accum, shrink;
	shrink_accum(inv_accum, shrink);

#ifdef OPENMP
#pragma omp parallel for
#endif
	for (int y = 0; y < ResY; y++)
	{
#ifndef OPENMP
		int thread = 0;
#else
		int thread = omp_get_thread_num();
#endif
		for (int x = 0; x < ResX; x++)
		{
			image[x + y * ResX] *= shrink;
			Ray ray = cam->getRay((float) x + mtrand[thread]->rand(),
					(float) y + mtrand[thread]->rand());
			HitRec rec = accel->intersect(ray);

			Vec3 color;
			if (rec.id != -1)
			{
				switch (shader)
				{
				case 1:
					color = shade_debug(ray, rec);
					break;
				case 2:
					color = shade_simple(ray, rec);
					break;
				case 3:
					color = shade_path(ray, rec, 0, thread);
					break;
				default:
					color = shade_simple(ray, rec);
					break;
				}
			}
			else
			{
				color = scene->getEnvironment(ray.dir);
			}
			image[x + y * ResX] += color * inv_accum;
		}
	}
}

Vec3 Render::shade_debug(Ray &ray, HitRec &rec)
{
	Vec3 normal = scene->getShadingNormal(ray, rec.id);
	normal.abs();
	return normal;
}

Vec3 Render::shade_simple(Ray &ray, HitRec &rec)
{
	Vec3 normal = scene->getShadingNormal(ray, rec.id);
	float cos = fabsf(normal * ray.dir);
	Vec3 color = cos * scene->material[scene->mat_index[rec.id]].color_d;
	return color;
}

Vec3 Render::shade_path(Ray &ray, HitRec &rec, int depth, int thread)
{
	// TODO 4.3 Implement a simple path tracer.

	Material material=scene->material[scene->mat_index[rec.id]];
	Vec3 hitpoint=ray.origin+ray.dir*rec.dist;
	//same call as in shade_simple
	Vec3 interpolNormal=scene->getShadingNormal(ray, rec.id);

	Ray newRay;
	newRay.origin=hitpoint;
	newRay.tmin = RAY_EPS;
	newRay.tmax = RAY_MAX;
	//shoot new ray like above (line 68)
	HitRec rayShoot = accel->intersect(newRay);

	//if depth too big
	if(depth>10){
		return Vec3(0.0f);
	}


	//if triangle hit
	if(rec.id!=-1){
		if (material.color_e[0] != 0 && material.color_e[1]!=0 && material.color_e[2]!=0) {
				return material.color_e;
			}

		//front or back facing
		if (interpolNormal * ray.dir > 0.0f){
				interpolNormal[0] = interpolNormal[0]*-1;
				interpolNormal[1] = interpolNormal[1]*-1;
				interpolNormal[2] = interpolNormal[2]*-1;
			}
		material.diffuse(newRay.dir, interpolNormal, mtrand[thread]->rand(),mtrand[thread]->rand());


		if (rayShoot.id == -1){
			return Vec3::product(material.color_d, scene->getEnvironment(newRay.dir));
		}
	}


	return Vec3::product(material.color_d, shade_path(newRay, rayShoot, depth++, thread));
}

inline void Render::shrink_accum(float &inv_accum, float &shrink)
{
	if (!cam->moved)
	{
		accum_index++;
		inv_accum = 1.0f / (float) accum_index;
		shrink = (accum_index - 1) * inv_accum;
	}
	else
	{
		accum_index = 1;
		inv_accum = 1.0f;
		shrink = 0.0f;
	}
}

