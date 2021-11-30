/**
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#ifndef RTSTRUCTS_H
#define RTSTRUCTS_H

#include "utils/vec.h"

#define RAY_MAX FLT_MAX
#define RAY_EPS 0.0001f

inline float minf(const float a, const float b)
{
	return a < b ? a : b;
}
inline float maxf(const float a, const float b)
{
	return a > b ? a : b;
}

/**
 * Ray containing an intersection interval.
 */
struct Ray
{
	/// Origin of the ray.
	Vec3 origin;
	// Direction of the ray
	Vec3 dir;
	/// Minimum distance for intersection (inclusive).
	float tmin;
	/// Maximum distance for intersection (inclusive).
	float tmax;

	/**
	 * Standard constructor which does not initialize the ray.
	 */
	inline Ray()
	{
	}

	/**
	 * Initializes the ray.
	 * @param origin See origin.
	 * @param dir See dir.
	 * @param tmin See tmin.
	 * @param tmax See tmax.
	 */
	inline Ray(const Vec3 &origin, const Vec3 &dir, const float tmin,
			const float tmax) :
			origin(origin), dir(dir), tmin(tmin), tmax(tmax)
	{
	}
};

/**
 * Stores information about a hit of a ray against a surface.
 */
struct HitRec
{
	/// Distance from the ray origin to the intersection point.
	float dist;
	/// Id of the hit surface (e.g. triangle).
	int id;

	/// Initializes the record as not hitting anything.
	inline HitRec()
	{
		dist = RAY_MAX;
		id = -1;
	}

	/** Initializes the record with given values.
	 * @param d See dist.
	 * @param i See id.
	 */
	inline HitRec(float d, int i)
	{
		dist = d;
		id = i;
	}
};

/**
 * Axis aligned bounding box.
 */
struct AABB
{
	/// Corners of the bounding box with the minimum and the maximum elements.
	Vec3 bounds[2];

	/**
	 * Initializes the bounding box with -max volume. This way, it can be
	 * extended by vertices or other AABBs easily.
	 */
	inline AABB()
	{ // empty box
		bounds[0] = Vec3(FLT_MAX);
		bounds[1] = Vec3(-FLT_MAX);
	}

	/**
	 * Initializes the bounding box with given values.
	 * @param bmin Initializes bounds[0]
	 * @param bmax Initializes bounds[1]
	 */
	inline AABB(Vec3 bmin, Vec3 bmax)
	{
		bounds[0] = bmin;
		bounds[1] = bmax;
	}

	/**
	 * Extends the AABB to contain another AABB.
	 * @param bb Other AABB to be contained.
	 */
	inline void extend(const AABB &bb)
	{
		bounds[0].minf(bb.bounds[0]);
		bounds[1].maxf(bb.bounds[1]);
	}

	/**
	 * Intersects a ray with the bounding box.
	 * @remarks outsources costly if operations for managing different
	 * ray directions so that they only have to be done once per ray
	 * and not once per node during BVH traversal.
	 * @param r Ray to intersect.
	 * @param intervalMin In/out parameter: Contains the minimum of the
	 * interval to check before the call and the minimum of the interval
	 * in which the ray is inside the check-interval and inside the AABB
	 * after the call.
	 * @param intervalMax In/out parameter: Contains the minimum of the
	 * interval to check before the call and the minimum of the interval
	 * in which the ray is inside the check-interval and inside the AABB
	 * after the call.
	 * @param invRayDir Inverse ray direction used to save costly divisions.
	 * @param raySign Contain different signs depending on the ray
	 * direction to be able to do an if-operation free box intersection.
	 * @returns True if the ray intersects with the AABB in the given
	 * interval.
	 */
	inline bool intersect(const Ray &r, float &intervalMin, float &intervalMax,
			const Vec3 &invRayDir, const unsigned int raySign[3][2]) const
	{
		const float xMin = (bounds[raySign[0][0]][0] - r.origin[0])
				* invRayDir[0];
		const float xMax = (bounds[raySign[0][1]][0] - r.origin[0])
				* invRayDir[0];

		const float yMin = (bounds[raySign[1][0]][1] - r.origin[1])
				* invRayDir[1];
		const float yMax = (bounds[raySign[1][1]][1] - r.origin[1])
				* invRayDir[1];

		const float zMin = (bounds[raySign[2][0]][2] - r.origin[2])
				* invRayDir[2];
		const float zMax = (bounds[raySign[2][1]][2] - r.origin[2])
				* invRayDir[2];

		intervalMin = maxf(intervalMin, xMin);
		intervalMin = maxf(intervalMin, yMin);
		intervalMin = maxf(intervalMin, zMin);

		intervalMax = minf(intervalMax, xMax);
		intervalMax = minf(intervalMax, yMax);
		intervalMax = minf(intervalMax, zMax);

		return !(intervalMin > intervalMax);
	}

	/**
	 * Returns the dimension of the bounding box's greatest extent.
	 * @param Dimesion of the bounding box's greatest extent.
	 */
	inline int getMaxAxis() const
	{
		Vec3 ex = bounds[1] - bounds[0];
		return ex.maxIndex();
	}

	/**
	 * Computes the center of the bounding box.
	 * @returns Center of the bounding box.
	 */
	inline Vec3 getCenter() const
	{
		return (bounds[1] + bounds[0]) * 0.5f;
	}
};

/**
 * Triangle.
 */
struct Triangle
{
	/// Vertices of the triangle (unordered).
	Vec3 v[3];
	// TODO 3.3 a) Extent Triangle to contain also normals.
	Vec3 n[3];

	/**
	 * Returns the extent of a triangle along a certain axis.
	 * @param axis Coordinate axis index for determining the extent.
	 * @param Out parameter: Minimum value covered by the triangle along the
	 * given axis.
	 * @param Out parameter: Maximum value covered by the triangle along the
	 * given axis.
	 */
	inline void getExtents(const int axis, float &xmin, float &xmax) const
	{
		xmin = minf(v[0][axis], v[1][axis]);
		xmax = maxf(v[0][axis], v[1][axis]);
		xmin = minf(xmin, v[2][axis]);
		xmax = maxf(xmax, v[2][axis]);
	}

	/**
	 * Returns the AABB of the triangle.
	 * @returns AABB containing all vertices of the triangle.
	 */
	inline AABB getAABB() const
	{
		AABB bbox;
		for (int axis = 0; axis < 3; axis++)
			getExtents(axis, bbox.bounds[0][axis], bbox.bounds[1][axis]);
		return bbox;
	}

	/**
	 * Returns the AABB of a set of triangles.
	 * @param tris Array containing all triangles which must be inside the
	 * AABB.
	 * @param nTris Number of triangles in tris.
	 * @returns AAB containing all vertices of all triangles.
	 */
	static AABB getAABB(const Triangle * const tris, const unsigned int nTris)
	{
		AABB bbox; // start with empty box

		for (unsigned int i = 0; i < nTris; i++)
		{
			for (int c = 0; c < 3; c++)
			{
				bbox.bounds[0][c] = minf(bbox.bounds[0][c], tris[i].v[0][c]);
				bbox.bounds[1][c] = maxf(bbox.bounds[1][c], tris[i].v[0][c]);
				bbox.bounds[0][c] = minf(bbox.bounds[0][c], tris[i].v[1][c]);
				bbox.bounds[1][c] = maxf(bbox.bounds[1][c], tris[i].v[1][c]);
				bbox.bounds[0][c] = minf(bbox.bounds[0][c], tris[i].v[2][c]);
				bbox.bounds[1][c] = maxf(bbox.bounds[1][c], tris[i].v[2][c]);
			}
		}
		return bbox;
	}

	/**
	 * Intersects a ray with the triangle.
	 * @param ray Ray to intersect with the triangle.
	 * @param rec Hit record containing the id of the intersected object and
	 * the distance to the hit. Should only be updated if the intersection
	 * point of this triangle is closer than the intersection that was stored
	 * in the hit record before.
	 * @param tri_id Triangle id of this triangle.
	 * @returns True if the ray intersects with the triangle in the interval
	 * set in the ray.
	 */
	inline bool intersect(const Ray &ray, HitRec &rec, const int tri_id) const
	{
		// Moeller-Trumbore based Triangle Intersect

		const Vec3 edge1 = v[1] - v[0];
		const Vec3 edge2 = v[2] - v[0];

		const Vec3 pvec = Vec3::cross(ray.dir, edge2);

		const float det = edge1 * pvec;
		const float invDet = 1.0f / det;

		const Vec3 tvec = ray.origin - v[0];

		const float alpha = (tvec * pvec) * invDet;

		if (!(0.0f <= alpha) || alpha > 1.0f)
			return false;

		const Vec3 qvec = Vec3::cross(tvec, edge1);

		const float beta = (ray.dir * qvec) * invDet;

		if (!(0.0f <= beta) || alpha + beta > 1.0f)
			return false;

		const float t = (edge2 * qvec) * invDet;

		if (!(ray.tmin < t) || rec.dist < t)
			return false;

		rec.dist = t;
		rec.id = tri_id;

		return true;
	}

	/**
	 * Returns the normal of the triangle.
	 * @remarks It is undefined if the normal points to the front or the back
	 * side of the triangle because the vertex order is undefined.
	 */
	inline Vec3 getNormal() const
	{
		const Vec3 edge1 = v[1] - v[0];
		const Vec3 edge2 = v[2] - v[0];

		Vec3 normal = Vec3::cross(edge1, edge2);
		normal.normalize();
		return normal;
	}

	/**
	 * Calculates the barycentric coordinates of a given point.
	 * @param ray Ray that hits the triangle.
	 * @returns Barycentric coordinates b.
	 * v[0] * b[0] + v[1] * b[1] + v[2] * b[2] = point
	 */
	inline Vec3 calculateBarycentric(const Ray &ray) const
	{

		// TODO: 3.3 b) Calculate the barycentric coordinates of the intersection of the given ray and the triangle.
		float lambda1, lambda2, lambda3;
		Vec3 b = v[1] - v[0];
		Vec3 c = v[2] - v[0];

		Vec3 p = Vec3::cross(ray.dir, c);

		float d = b * p;


		Vec3 t = ray.origin - v[0];
		Vec3 q = Vec3::cross(t, b);


		lambda2 = (t * p) / d;
		lambda3 = (ray.dir * q) /d;
		lambda1=1-lambda2-lambda3;

		return Vec3(lambda1, lambda2, lambda3);
	}

	/**
	 * Calculates a normal by interpolating between the normals of the
	 * vertices according to a given hit point.
	 * @param hitPoint Point where the triangle was hit by a ray.
	 * @returns Interpolated vertex normal.
	 */
	inline Vec3 getInterpolatedNormal(const Ray &ray) const
	{
		// TODO 3.3 c) Calculate an interpolated hitpoint normal using the barycentric coordinates.
		Vec3 bary=calculateBarycentric(ray);

		Vec3 interpol=(bary[0]*n[0]+bary[1]*n[1]+bary[2]*n[2]);
		float l = interpol.normalize();


		return interpol * l;
	}
};

/**
 * Point light
 */
struct Pointlight
{
	/// Position of the point light in the scene.
	Vec3 pos;
	/// Color of the pointlight.
	Vec3 color;
};

/**
 * Parallelogram shaped area light
 * @remarks Calculate a sample on the area light with
 * pos + extent1 * rand() + extent2 * rand().
 */
struct AreaLight
{
	/// One corner of the area light.
	Vec3 pos;
	/// Vector from pos to a neighbor corner.
	Vec3 extent1;
	/// Vector from pos to the other neighbor corner.
	Vec3 extent2;
	/// Radiance of the area light.
	Vec3 radiance;

	/**
	 * Returns the area of the area light with respect to the extent vectors.
	 */
	float getArea() const
	{
		// TODO 3.6 a) Calculate the area of the light source.
		float height=(extent1*(1-(extent1*extent2))).normalize();
		float extent2L=extent2.length();

		float area=height*extent2L;




		return area;
	}
};

#endif

