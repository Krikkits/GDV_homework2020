/**
 * Computer Graphics Lecture WS 2009/2010 Ulm University
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#ifndef RENDER_H
#define RENDER_H

#include "bvh.h"
#include "rtStructs.h"
#include "utils/vec.h"
#include "utils/MersenneTwister.h"

#include "scene.h"
#include "cam.h"
#include "material.h"

/**
 * Renderer
 */
struct Render
{
	/// Scene to be rendered.
	Scene *scene;
	/// Camera from which the scene is rendered.
	Cam *cam;
	/**
	 * Acceleration structure which supports generating HitRecords by
	 * shooting rays into the scene.
	 */
	BVH *accel;

	/// Array of random number generators - one per thread (since they are not reentrant).
	MTRand **mtrand;
	/// Number of threads / elements in mtrand.
	int numThreads;

	/// The image that is rendered.
	Vec3 *image;
	/// Width of the image in pixels.
	int ResX;
	/// Height of the image in pixels.
	int ResY;
	/// Number of pictures that have been accumulated in image.
	int accum_index;

	/**
	 * Initializes the scene completely according to a given scene.
	 */
	Render(Scene *scene);
	/**
	 * Tidies up everything.
	 */
	~Render();

	/**
	 * Updates image: Improves it if the camera was not moved or refreshes it
	 * otherwise.
	 * @param shader: Shader to use for the rendering process:
	 * 1=debug, 2=simple, 3=path
	 */
	void render(int shader);

	/** Calculates the shading of the surface in debugging colors.
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 */
	inline Vec3 shade_debug(Ray &ray, HitRec &rec);

	/** Calculates the shading with simple smooth shading.
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 */
	inline Vec3 shade_simple(Ray &ray, HitRec &rec);\
	/** Recursively calculates the shading with a path tracer which
	 * incorporates indirect light.
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 * @param recursion depth of the path tracing.
	 * @param thread Thread id used for choosing the right Twister in mtrand.
	 */
	inline Vec3 shade_path(Ray &ray, HitRec &rec, int depth, int thread);

	/**
	 * Updates the accum_index.
	 * @param inv_accom Output parameter: Contribution of one single rendering to image.
	 * @param shrink Output parameter: Shrink factor that has to be multiplied to
	 * image to keep it normalized when the next rendering is added.
	 */
	inline void shrink_accum(float &inv_accum, float &shrink);
};

#endif

