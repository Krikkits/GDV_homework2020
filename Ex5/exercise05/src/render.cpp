/**
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
	mtrand = new MTRand*[1];
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

			// TODO 5.4 e) Remove the random pixel sampling to have only one sample per pixel.
			Ray ray = cam->getRay((float) x,(float) y);
			//original:
			/*
			Ray ray = cam->getRay((float) x + mtrand[thread]->rand(),
					(float) y + mtrand[thread]->rand());
			*/

			HitRec rec = accel->intersect(ray);

			Vec3 color(0.0f, 0.0f, 0.0f);
			if (rec.id != -1)
			{
				switch (shader)
				{
				case 1:
					color = shade_debug_normal(ray, rec);
					break;
				case 2:
					color = shade_debug_uv(ray, rec);
					break;
				case 3:
					color = shade_debug_miplevel(ray, rec);
					break;
				case 4:
					color = shade_noshading(ray, rec);
					break;
				case 5:
					color = shade_simple(ray, rec);
					break;
				case 6:
					color = shade_path(ray, rec, 0, thread);
					break;
				default:
					color = shade_noshading(ray, rec);
					break;
				}
			}
			else if (scene->environment != NULL)
			{
				color = scene->getEnvironment(ray.dir);
			}
			image[x + y * ResX] += color * inv_accum;
		}
	}
}

float Render::getMipLevel(float distance)
{
	// TODO 5.4 d) Calculate the mip level from the given distance and use the result for the getTextureColor(...) calls.


	return log2(distance);
}

Vec3 Render::shade_debug_normal(Ray &ray, HitRec &rec)
{
	Vec3 normal = scene->getShadingNormal(ray, rec.id);
	normal.abs();
	return normal;
}

Vec3 Render::shade_debug_uv(Ray &ray, HitRec &rec)
{
	Vec2 uv = scene->getTextureCoordinates(ray, rec.id);
	return Vec3(uv.x, uv.y, 0.0f);
}

Vec3 Render::shade_debug_miplevel(Ray &ray, HitRec &rec)
{
	int level = (int)getMipLevel(rec.dist);

	Vec3 color(0.0f);
	color[level % 3] = 1.0f;

	return color;
}

Vec3 Render::shade_noshading(Ray &ray, HitRec &rec)
{
	Vec3 color(scene->material[scene->mat_index[rec.id]].color_d);
	Vec2 coords=scene->getTextureCoordinates(ray, rec.id);
	Vec3 texColor=scene->material[scene->mat_index[rec.id]].Material::GetTextureColor(coords);

	// TODO 5.3 b) Multiply color with the texture color by calling Material::getTextureColor(coords).
	Vec3 newcolor=Vec3::product(color, texColor);

	// TODO 5.4 d) Add the second parameter to Material::getTextureColor(...).
	int c=4;
	//int level = (int)getMipLevel(rec.dist);

	//set level for test
	int level=3;
	Vec3 mipColor=scene->material[scene->mat_index[rec.id]].Material::GetTextureColor(coords, log2(level*c));

	return Vec3::product(color, mipColor);
	//return newcolor;
}

Vec3 Render::shade_simple(Ray &ray, HitRec &rec)
{
	Vec3 normal = scene->getShadingNormal(ray, rec.id);
	float cos = fabsf(normal * ray.dir);
	Vec3 color = cos * scene->material[scene->mat_index[rec.id]].color_d;
	Vec2 coords=scene->getTextureCoordinates(ray, rec.id);
	Vec3 texColor=scene->material[scene->mat_index[rec.id]].Material::GetTextureColor(coords);

	// TODO 5.3 b) Multiply color with the texture color by calling Material::getTextureColor(coords).
	Vec3 newcolor=Vec3::product(color, texColor);

	// TODO 5.4 d) Add the second parameter to Material::getTextureColor(...).
	int c=4;
	//int level = (int)getMipLevel(rec.dist);
	//set level for test
	int level=3;
	Vec3 mipColor=scene->material[scene->mat_index[rec.id]].Material::GetTextureColor(coords, log2(level*c));

	return Vec3::product(color, mipColor);
	//return newcolor;
}

Vec3 Render::shade_path(Ray &ray, HitRec &rec, int depth, int thread)
{
	if (depth > 5)
		return Vec3(0.0f);

	Material &mat = scene->material[scene->mat_index[rec.id]];

	if (mat.color_e[0] != 0.0f && mat.color_e[1] != 0.0f
			&& mat.color_e[2] != 0.0f)
		return mat.color_e;

	Ray newRay;
	newRay.origin = ray.origin + ray.dir * rec.dist;
	Vec3 hitNormal = scene->getShadingNormal(ray, rec.id);
	if (hitNormal * ray.dir > 0.0f)
		hitNormal *= -1.0f;
	mat.diffuse(newRay.dir, hitNormal, mtrand[thread]->rand(),
			mtrand[thread]->rand());
	newRay.tmin = RAY_EPS;
	newRay.tmax = RAY_MAX;

	Vec3 color = mat.color_d;

	// TODO 5.3 b) Multiply color with the texture color by calling Material::getTextureColor(coords).
	Vec2 coords=scene->getTextureCoordinates(ray, rec.id);
	Vec3 texColor=scene->material[scene->mat_index[rec.id]].Material::GetTextureColor(coords);
	Vec3 newcolor=Vec3::product(color, texColor);

	// TODO 5.4 d) Add the second parameter to Material::getTextureColor(...).
	int c=4;
	//int level = (int)getMipLevel(rec.dist);
	//set level for test
	int level=3;
	Vec3 mipColor=scene->material[scene->mat_index[rec.id]].Material::GetTextureColor(coords, log2(level*c));
	mipColor= Vec3::product(color, mipColor);

	HitRec newRec = accel->intersect(newRay);
	if (newRec.id == -1)
		return Vec3::product(color, scene->getEnvironment(newRay.dir));

	//return Vec3::product(mipColor, shade_path(newRay, newRec, depth + 1, thread));
	return Vec3::product(newcolor, shade_path(newRay, newRec, depth + 1, thread));
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

