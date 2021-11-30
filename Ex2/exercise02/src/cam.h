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
	Vec3 origin;
	Vec3 dir;
	Vec3 up;
	Vec3 right;
	float focalLength;
	int resX;
	int resY;

	/**
	 * Initializes the camera.
	 * @param box See box.
	 * @param ResX See ResX.
	 * @param ResY See ResY.
	 */
	PerspCam(const Vec3 &origin, const Vec3 &dir, const Vec3 &up,
			float focalLength, int resX, int resY) :
			origin(origin), dir(dir), up(up), focalLength(focalLength), resX(
					resX), resY(resY)
	{
		right = Vec3::cross(dir, up);
	}

	/**
	 * Returns a gamera ray given pixel coordinates on the image plane.
	 * @param x X coordinate on the image plane in [0, ResX-1].
	 * @param y Y coordinate on the image plane in [0, ResY-1].
	 * @returns Ray according to the given pixel coordinates.
	 */
	Ray getRay(const float x, const float y) const;
};

inline Ray PerspCam::getRay(const float x, const float y) const
{
	float xOffset = x - (float) resX / 2.0f;
	float yOffset = y - (float) resY / 2.0f;

	Vec3 rayDir = dir * focalLength + right * xOffset + up * yOffset;
	rayDir.normalize();
	return Ray(origin, rayDir, 0.0f, RAY_MAX);
}

#endif

