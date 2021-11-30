/**
 * Computer Graphics Lecture WS 2009/2010 Ulm University
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include "utils/vec.h"
#include "utils/MersenneTwister.h"
#include "rtStructs.h"
#include <iostream>

/**
 * Simple material class
 */
struct Material
{
	/// Diffuse color
	Vec3 color_d;
	/// Emitting color
	Vec3 color_e;

	/**
	 * Default constructor not initializing anything.
	 */
	Material()
	{
	}

	/**
	 * Initializes the object completely.
	 */
	Material(const Vec3 &d, const Vec3 &e)
	{
		color_d = d;
		color_e = e;
	}

	/**
	 * Generates a bounce ray direction which is statistically distributed
	 * according to a diffuse surface.
	 * @remarks t0 and t1 have to be independent!
	 * @param out Out parameter, the bounced ray direction.
	 * @param n Normal of the surface.
	 * @param t0 random variable between 0 and 1.
	 * @param t1 random variable between 0 and 1.
	 */
	inline static void diffuse(Vec3 &out, const Vec3 &n, const float t0,
			const float t1);

	/**
	 * Generates a bounce ray direction for a mirroring surface.
	 * @param out Out parameter, the bounced ray direction.
	 * @param in Direction of the incoming ray.
	 * @param n Normal of the surface.
	 */
	inline static void mirror(Vec3 &out, const Vec3 &in, const Vec3 &n);
};

inline void Material::diffuse(Vec3 &out, const Vec3 &n, const float t0,
		const float t1)
{
	Vec3 u, v, w;
	Vec3::onb(u, v, w, n);

	const float phi = 2.0f * M_PI * t0;
	const float sqrtu = sqrtf(t1);
	const float x = sqrtu * cosf(phi);
	const float y = sqrtu * sinf(phi);
	const float z = sqrtf(1.0f - t1);

	out = x * u + y * v + z * w;
}

inline void Material::mirror(Vec3 &out, const Vec3 &in, const Vec3 &n)
{
	const float dot2 = 2.0f * (in * n); // assumes in points onto surface -> (in * normal < 0.0f)

	out = in - (dot2 * n);
}

#endif

