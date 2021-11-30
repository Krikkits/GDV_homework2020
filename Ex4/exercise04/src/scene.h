/**
 * Computer Graphics Lecture WS 2009/2010 Ulm University
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#ifndef SCENE_H
#define SCENE_H

#include "material.h"
#include "cam.h"
#include "rtStructs.h"
#include <iostream>

/**
 * Scene to be rendered.
 */
struct Scene
{
	/// All triangles contained in the scene.
	Triangle * triangles;
	/**
	 * All normals contained in the scene.
	 * normals[t*3 + n] contains the nth normal of the tth triangle.
	 */
	Vec3 * normals;
	/// Number of triangles in the scene.
	int num_tris;

	/// All materials in the scene.
	Material * material;
	/// Number of materials in material.
	int num_material;
	/**
	 * mat_index[t] contains the index of the material in materials which
	 * should be used for rendering the tth triangle.
	 */
	int * mat_index;

	/// Environment map of the scene.
	Vec3 * environment;
	/// Width of the environment map in pixels.
	int env_x;
	/// Height of the environment map in pixels.
	int env_y;

	/// Camera that will be used for rendering.
	Cam * cam;

	/**
	 * Initializes the scene from files.
	 * @param sceneFile Relative path to the scene files without extension (and ".").
	 * @param envFile Relative path to the environment map file (.hdr format).
	 * May be NULL if no environment map should be loaded.
	 */
	Scene(const char* sceneFile, const char* envFile);
	~Scene();

	/**
	 * Loads the scene.
	 * @param file Relative path to the scene files without extension (and ".").
	 */
	bool LoadScn(const char * file);
	/**
	 * Loads an environment map.
	 * @param file Relative path to the environment map file (.hdr format).
	 * May be NULL if no environment map should be loaded.
	 */
	bool LoadEnv(const char * file);

	/**
	 * Retrieves a smooth shading normal from the scene.
	 * @param ray Ray that hits the shaded point.
	 * @param tri_id Index of the hit triangle.
	 * @returns Smooth shading normal.
	 */
	Vec3 getShadingNormal(const Ray &ray, const int tri_id) const;
	/**
	 * Retrieves the environment color for a given direction.
	 * @param d Ray direction.
	 * @returns HDR environment map color for the given direction.
	 */
	Vec3 getEnvironment(const Vec3 &d) const;
};

inline Vec3 Scene::getShadingNormal(const Ray &ray, const int tri_id) const
{
	float alpha, beta, gamma;
	triangles[tri_id].getBarycentric(ray, alpha, beta);
	gamma = 1.0f - alpha - beta;
	Vec3 n = gamma * normals[tri_id * 3] + alpha * normals[tri_id * 3 + 1]
			+ beta * normals[tri_id * 3 + 2];
	n.normalize();
	return n;
}

inline float sign(float a)
{
	if (a >= 0.0f)
		return 1.0f;
	else
		return -1.0f;
}

inline Vec3 Scene::getEnvironment(const Vec3 &drot90) const
{
	if (environment == 0)
		return Vec3(0.0f);

	Vec3 d(drot90.x, -drot90.z, drot90.y);

	float inv_sum = 1.0f / (fabs(d[0]) + fabs(d[1]) + fabs(d[2]));
	float px = d[0] * inv_sum;
	float py = d[1] * inv_sum;
	float pz = d[2] * inv_sum;

	float u, v;
	if (pz >= 0.0f)
	{
		u = px * 0.5f + 0.5f;
		v = py * 0.5f + 0.5f;
	}
	else
	{
		u = sign(px) * (1.0f - fabs(py)) * 0.5f + 0.5f;
		v = sign(py) * (1.0f - fabs(px)) * 0.5f + 0.5f;
	}
	int x = (int) (u * (float) env_x);
	int y = (int) (v * (float) env_y);

	x = x > 0 ? x : 0;
	x = x < env_x ? x : env_x - 1;
	y = y > 0 ? y : 0;
	y = y < env_y ? y : env_y - 1;

	return environment[x + y * env_x];
}

#endif

