#ifndef CAM_H
#define CAM_H

#include "utils/vec.h"
#include "rtStructs.h"

/**
 * Orthographic camera watching a AABB into z direction.
 */
struct Cam
{
	/// Volume watched by the camera.
	AABB box;
	/// Resolution of the image plane in horizontal direction.
	int ResX;
	/// Resolution of the image plane in vertical direction.
	int ResY;
	/// Extent of the observed volume.
	Vec3 extent;
	/// Viewing direction of the camera (always (0,0,1)).
	Vec3 dir;

	/// Contains the step size per pixel in horizontal direction.
	float xstep;
	/// Contains the step size per pixel in vertical direction.
	float ystep;

	/**
	 * Initializes the camera.
	 * @param box See box.
	 * @param ResX See ResX.
	 * @param ResY See ResY.
	 */
	Cam(const AABB &box, const int ResX, const int ResY) :
			box(box), ResX(ResX), ResY(ResY)
	{
		extent = box.bounds[1] - box.bounds[0];

		xstep = extent[0] / (float) ResX;
		ystep = extent[1] / (float) ResY;

		dir = Vec3(0.0f, 0.0f, 1.0f);
	}

	/**
	 * Returns a camera ray given pixel coordinates on the image plane.
	 * @param x X coordinate on the image plane in [0, ResX-1].
	 * @param y Y coordinate on the image plane in [0, ResY-1].
	 * @returns Ray according to the given pixel coordinates.
	 */
	Ray getRay(const float x, const float y) const;
};

inline Ray Cam::getRay(const float x, const float y) const
{
	Vec3 eye = box.bounds[0];
	eye[0] += x * xstep;
	eye[1] += y * ystep;
	return Ray(eye, dir, 0.0f, RAY_MAX);
}

/**
 * Perspective camera
 */
struct PerspCam
{
	/// Position of the camera
	Vec3 pos;
	/// Viewing direction
	Vec3 dir;
	/// Up direction
	Vec3 up;
	/// Right direction
	Vec3 right;
	/// Resolution of the rendertarget
	int ResX, ResY;
	/// Move distance on the image plane per pixel
	float xstep, ystep;
	/// Edge of the image plane
	Vec3 planeEdge;

	/**
	 * Initializes the perspective cam for ray generation.
	 * @param pos See pos.
	 * @param dir See dir.
	 * @param up See up.
	 * @param fov Vertical field of view in degrees.
	 * @param resX See ResX.
	 * @param resY See ResY.
	 */
	PerspCam(const Vec3 &pos, const Vec3 &dir, const Vec3 &up, const float &fov,
			const int &resX, const int &resY) :
			pos(pos), dir(dir), up(up), ResX(resX), ResY(resY)
	{
		right = Vec3::cross(up, dir);
		float maxY = tanf((float) fov * 3.1415f / 360.0f);
		float maxX = maxY * (float) ResX / float(ResY);
		xstep = maxX * 2.0f / (float) ResX;
		ystep = maxY * 2.0f / (float) ResY;
		planeEdge = pos + dir - right * xstep * ((float) ResX / 2.0f)
				- up * ystep * ((float) ResY / 2.0f);
	}

	/**
	 * Returns a camera ray given pixel coordinates on the image plane.
	 * @param x X coordinate on the image plane in [0, ResX-1].
	 * @param y Y coordinate on the image plane in [0, ResY-1].
	 * @returns Ray according to the given pixel coordinates.
	 */
	Ray getRay(const float x, const float y) const;
};

inline Ray PerspCam::getRay(const float x, const float y) const
{
	return Ray(pos, planeEdge + right * xstep * x + up * ystep * y - pos, 0.0f,
			RAY_MAX);
}

#endif

