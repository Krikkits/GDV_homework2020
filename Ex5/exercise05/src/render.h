/**
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

	/**
	 * Calculates the mip level for a certain distance.
	 * @remarks does not adapt to different texture sizes in world space.
	 * @param distance Distance from the ray origin to the texture.
	 * @returns Estimate which is at least proportional to the real mip level that
	 * should be choosen.
	 */
	inline float getMipLevel(float distance);

	/** Returns the surface normal as shade.
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 * @returns Normal of the surface.
	 */
	inline Vec3 shade_debug_normal(Ray &ray, HitRec &rec);

	/** Returns the uv texture coordinates in the red/green channels.
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 * @returns uv coordinates.
	 */
	inline Vec3 shade_debug_uv(Ray &ray, HitRec &rec);

	/** Color varies with the used mip level (one primary color per mip level).
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 * @returns Blending between primary colors.
	 */
	inline Vec3 shade_debug_miplevel(Ray &ray, HitRec &rec);

	/**
	 * Returns the diffuse color of the surface.
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 * @returns The diffuse color of the surface.
	 */
	inline Vec3 shade_noshading(Ray &ray, HitRec &rec);

	/** Calculates the shading with simple smooth shading.
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 * @returns Diffuse surface as lit by the camera.
	 */
	inline Vec3 shade_simple(Ray &ray, HitRec &rec);

	/** Recursively calculates the shading with a path tracer which
	 * incorporates indirect light.
	 * @param Ray along which the the shading has to be calculated.
	 * @param Defines where the ray hits the scene. Must be a valid hit!
	 * @param recursion depth of the path tracing.
	 * @param thread Thread id used for choosing the right Twister in mtrand.
	 * @returns Surface color as lit directly and indirectly by light sources
	 * in the scene and by the environment.
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

