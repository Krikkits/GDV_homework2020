/**
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#ifndef MATERIAL_H
#define MATERIAL_H

#include "utils/vec.h"
#include "utils/MersenneTwister.h"
#include "rtStructs.h"
#include <iostream>
#include <assert.h>

/**
 * Texture which supports bilinear and trilinar sampling with mipmaps.
 */
struct Texture
{
	/// Width of the texture in pixels.
	int ResX;
	/// Height of the texture in pixels.
	int ResY;
	/// Number of mip levels.
	int MipLevels;
	/**
	 * RGB data of all mip levels of the texture.
	 * data[l][y*mipWidth+x] contains the color of the sample at texel (x, y) of mip level l.
	 * Mip level 0 has the full texture resolution while the highest level is 1x1 pixel in size.
	 */
	Vec3 **data;

	/**
	 * Calculates the width of a certain mip level.
	 * @remarks Only defined for 1..MipLevels-1 !
	 * @param level Mip level for which the width should be calculated.
	 * @returns Width in pixels of the mip level.
	 */
	inline int MipResX(int level)
	{
		// Assures precondition in debug mode.
		assert(level >= 0 && level < MipLevels);

		// TODO 5.4 a) Calculate a mip level's width given the level.
		if(level==0){
			return ResX;
		}
		for(int i=0;i<level;i++){
			ResX=ResX/2;
		}
		return ResX;
	}

	/**
	 * Calculates the height of a certain mip level.
	 * @remarks Only defined for 1..MipLevels-1 !
	 * @param level Mip level for which the height should be calculated.
	 * @returns Height in pixels of the mip level.
	 */
	inline int MipResY(int level)
	{
		// Assures precondition in debug mode.
		assert(level >= 0 && level < MipLevels);

		// TODO 5.4 a) Calculate a mip level's height given the level.
		if(level==0){
			return ResY;
		}

		for(int i=0;i<level;i++){
			ResY=ResY/2;
		}

		return ResY;
	}

	/**
	 * Initializes the texture and all its mip levels from given data.
	 * @param ResX Width of the texture in texels.
	 * @param ResY Height of the texture in texels.
	 * @param nData Data of the texture. One Vec3 per texel, rowwise.
	 */
	Texture(const int ResX, const int ResY, Vec3 *nData) :
			ResX(ResX), ResY(ResY)
	{
		// TODO 5.4 b) Calculate the amount of mip levels necessary to have 1x1 texel on the smallest level.
		MipLevels = log2(ResX);
		//MipLevels=1;
		data = new Vec3*[MipLevels];
		data[0] = nData;

		// TODO 5.4 b) Calculate all mip levels
		int lastRes;
		int res;


		for (int j=1; j<MipLevels; j++){

			lastRes=MipResX(j-1);
			res=MipResX(j);
			data[j] = new Vec3[res * res]; //allocate memory
			for (int y=0; y<res; y++){
				for (int x=0; x<res; x++){
					//skipping already sampled pixels
					if (((x%2) == 0) && ((y%2) == 0)){

					data[j][y * res + x] = (data[j-1][y * lastRes + x] +
											data[j-1][y * lastRes + x+1] +
											data[j-1][(y+1) * lastRes + x] +
											data[j-1][(y+1) * lastRes + x+1])/4;
					//cout<<"sample: "<<data[j][y * res + x].x<<endl;
					}
				}
			}

		}
	}

	/**
	 * Copy constructor that performs deep copies.
	 */
	Texture(const Texture &original) :
			ResX(original.ResX), ResY(original.ResY), MipLevels(
					original.MipLevels)
	{
		data = new Vec3*[MipLevels];
		for (int l = 0; l < MipLevels; l++)
		{
			if (original.data[l])
			{
				data[l] = new Vec3[MipResX(l) * MipResY(l)];
				for (int p = 0; p < MipResX(l) * MipResY(l); p++)
					data[l][p] = original.data[l][p];
			}
			else
				data[l] = 0;
		}
	}

	/**
	 * Destroys all the data belonging to the texture.
	 */
	~Texture()
	{
		for (int l = 0; l < MipLevels; l++)
			if (data[l])
				delete[] data[l];
		delete[] data;
	}

	/**
	 * Returns a color sample from the texture or its mipmaps.
	 * @remarks The implementation wraps texture coordinates.
	 * @param coords Texture coordinates. (0,0) is the upper left corner, (1,1) is bottom right.
	 * @param level Mip level from which the sample should be calculated.
	 * @returns Color sample at the given position.
	 */
	Vec3 GetMipmappedColor(Vec2 coords, int level)
	{
		// TODO 5.4 c) Return the sample at coords from the given level (from data[level] instead of data[0])

		return data[level][(int) coords[1] * MipResX(level) + (int) coords[0]];
	}

	/**
	 * Returns a color sample from the texture.
	 * @remarks The implementation wraps texture coordinates.
	 * @param coords Texture coordinates. (0,0) is the upper left corner, (1,1) is bottom right.
	 * @param level Mip level which will be used for trilinear interpolation.
	 * @returns Color sample at the given position.
	 */
	Vec3 GetColor(Vec2 coords)
	{
		coords[0] = fmod(coords[0], 1.0f);
		coords[1] = -fmod(coords[1], 1.0f);
		if (coords[0] < 0.0f)
			coords[0] += 1.0f;
		if (coords[1] < 0.0f)
			coords[1] += 1.0f;

		float c0 = coords[0] * (float) ResX;
		float c1 = coords[1] * (float) ResY;

		return data[0][(int) c1 * ResX + (int) c0];
	}
};

/**
 * Simple material class
 */
struct Material
{
	/// Diffuse color
	Vec3 color_d;
	/// Emitting color
	Vec3 color_e;

	/// Texture of the material.
	Texture * tex;

	/**
	 * Default constructor not initializing anything.
	 */
	Material() :
			tex(0)
	{
	}

	/**
	 * Copy constructor which handles deep texture copies.
	 */
	Material(const Material &original) :
			color_d(original.color_d), color_e(original.color_e)
	{
		if (original.tex)
			tex = new Texture(*original.tex);
		else
			tex = 0;
	}

	/**
	 * Initializes the object completely.
	 */
	Material(const Vec3 &d, const Vec3 &e, Texture *tex = 0) : // tex defaults to 0 if not specified
			color_d(d), color_e(e), tex(tex)
	{
	}

	/**
	 * Destroys everything that was referenced from the material.
	 */
	~Material()
	{
		if (tex)
			delete tex;
	}

	/**
	 * Assignment operator which handles deep texture copies.
	 */
	Material& operator=(const Material &original)
	{
		color_d = original.color_d;
		color_e = original.color_e;
		if (tex)
			delete tex;
		if (original.tex)
			tex = new Texture(*original.tex);
		else
			tex = 0;

		return *this;
	}

	/**
	 * Returns the texture color or white if there is no texture.
	 * @remarks This way, multiplication with the return value of this method is always ok.
	 * @param coords Texture coordinates with x, y in (0, 1).
	 * @returns Color at the requested position or white (1, 1, 1) if there is non texture.
	 */
	Vec3 GetTextureColor(const Vec2 &coords) // mipLevel defaults to 0 if not specified.
	{
		if (!tex)
			return Vec3(1.0f);
		return tex->GetColor(coords);
	}

	/**
	 * Returns the mipmapped texture color or white if there is no texture.
	 * @remarks This way, multiplication with the return value of this method is always ok.
	 * @param coords Texture coordinates with x, y in (0, 1).
	 * @param mipLevel Mipmap level to fetch a texture sample from.
	 * @returns Color at the requested position or white (1, 1, 1) if there is non texture.
	 */
	Vec3 GetTextureColor(const Vec2 &coords, int mipLevel)
	{
		if (!tex)
			return Vec3(1.0f);
		return tex->GetMipmappedColor(coords, mipLevel);
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

