/**
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#ifndef VEC_H
#define VEC_H

#include <cmath>
#include <cfloat>

/**
 * 3D vector.
 */
struct Vec3
{
	/// X component
	float x;
	/// Y component
	float y;
	/// Z component
	float z;

	/**
	 * Does not initialize the vector.
	 */
	Vec3();
	/**
	 * Initializes all components of the vector with the same value.
	 * @param a Initializes all components of the vector.
	 */
	Vec3(const float a);
	/**
	 * Initializes the components of the vector from an array.
	 * @param a Array used for initialization.
	 */
	Vec3(const float * const a);
	/**
	 * Initializes the components from three scalars.
	 * @param a Initializes x.
	 * @param b Initializes y.
	 * @param c Initializes z.
	 */
	Vec3(const float a, const float b, const float c);

	/**
	 * Returns the negative vector.
	 * @returns Negated vector.
	 */
	Vec3 operator -() const;
	/**
	 * Calculates the sum of two vectors.
	 * @param v Vector to add.
	 * @returns Sum of the current vector and v.
	 */
	Vec3 operator +(const Vec3 &v) const;
	/**
	 * Adds another vector to the current vector.
	 * @param v Other vector to add.
	 * @returns Reference to the current vector after adding v.
	 */
	Vec3& operator +=(const Vec3 &v);
	/**
	 * Calculates the difference between the current vector and v.
	 * @param v Vector to subtract.
	 * @returns Current vector minus v.
	 */
	Vec3 operator -(const Vec3 &v) const;
	/**
	 * Subtracts another vector from the current vector.
	 * @param v Other vector to add.
	 * @returns Reference to the current vector after subtracting v.
	 */
	Vec3& operator -=(const Vec3 &v);
	/**
	 * Calculates the multiplication of the current vector with a scalar.
	 * @param a Scalar to multiply with the current vector.
	 * @returns Current vector scaled with a.
	 */
	Vec3 operator *(const float a) const;
	/**
	 * Scales the current vector (v * a).
	 * @param a Scaling factor.
	 * @returns Reference to the current vector after scaling.
	 */
	Vec3& operator *=(const float a);
	/**
	 * Calculates the multiplication of a scalar with the current vector.
	 * @param a Scalar to multiply with the current vector.
	 * @param Current vector will be placed here by the compiler.
	 * @returns Current vector scaled with a.
	 */
	friend Vec3 operator *(const float a, const Vec3 &v);
	/**
	 * Calculates a componentwise division of the current vector by a scalar.
	 * @param a Scalar for componentwise division.
	 * @returns Current vector inversely scaled with a.
	 */
	Vec3 operator /(const float a) const;
	/**
	 * Inversely scales the current vector.
	 * @param a Scalar for componentwise division.
	 * @returns Reference to the current vector after inverse scaling.
	 */
	Vec3& operator /=(const float a);

	/**
	 * Calculates the dot product.
	 * @param v Other vector used to calculate the dot product.
	 * @returns Dot product of the current vector and v.
	 */
	float operator *(const Vec3 &v) const; // dot product

	/**
	 * Allows read/write access on the vector components by index.
	 * @param i Index between 0 (x) and 2 (z).
	 * @returns Reference to the vector component corresponding to i.
	 */
	float& operator[](const int i);
	/**
	 * Allows read access on the vector components by index.
	 * @param i Index between 0 (x) and 2 (z).
	 * @returns Vector component corresponding to i.
	 */
	float operator[](const int i) const;

	/**
	 * Calculates the squared length of the vector.
	 * @returns Squared length of the vector.
	 */
	float length_sqr() const;
	/**
	 * Calculates the length of the vector.
	 * @returns Length of the vector.
	 */
	float length() const;
	/**
	 * Normalizes the vector.
	 */
	float normalize();

	/**
	 * Sets each component of the vector to its absolute value.
	 */
	void abs();

	/**
	 * Maximizes each component of the vector individually.
	 * @v Other vector used for maximization of the current vector's
	 * components.
	 */
	void minf(const Vec3 &v);
	/**
	 * Minimizes each component of the vector individually.
	 * @v Other vector used for minimization of the current vector's
	 * components.
	 */
	void maxf(const Vec3 &v);

	/**
	 * Returns the index of the component with the maximum value.
	 * @returns Index of the component with the maximum value.
	 */
	int maxIndex();

	/**
	 * Calculates the cross product of two vectors.
	 * @param v Left vector for the cross product.
	 * @param w Right vector for the cross product.
	 * @returns v x w.
	 */
	static Vec3 cross(const Vec3 &v, const Vec3 &w);
	/**
	 * Calculates the componentwise product of two vectors.
	 * @param v First vector for computing the product.
	 * @param w Second vector for computing the product.
	 * @returns Vector vec with vec[i] = v[i] * w[i].
	 */
	static Vec3 product(const Vec3 &v, const Vec3 &w);

	/**
	 * Computes an orthonormal base.
	 * @param u Out parameter: Normalized and orthogonal to v and w.
	 * @param v Out parameter: Normalized and ortogonal to u and w.
	 * @param w Out parameter: Normalized, orthogonal to u and v and parallel to
	 * n.
	 * @param n Initialization direction for w.
	 */
	static void onb(Vec3 &u, Vec3 &v, Vec3 &w, const Vec3 &n);
};

inline Vec3::Vec3()
{
}

inline Vec3::Vec3(const float a) :
		x(a), y(a), z(a)
{
}

inline Vec3::Vec3(const float * const a)
{
	x = a[0];
	y = a[1];
	z = a[2];
}

inline Vec3::Vec3(const float a, const float b, const float c) :
		x(a), y(b), z(c)
{
}

inline Vec3 Vec3::operator -() const
{
	return Vec3(-x, -y, -z);
}

inline Vec3 Vec3::operator +(const Vec3 &v) const
{
	return Vec3(x + v.x, y + v.y, z + v.z);
}

inline Vec3& Vec3::operator +=(const Vec3 &v)
{
	x += v.x;
	y += v.y;
	z += v.z;

	return *this;
}

inline Vec3 Vec3::operator -(const Vec3 &v) const
{
	return Vec3(x - v.x, y - v.y, z - v.z);
}

inline Vec3& Vec3::operator -=(const Vec3 &v)
{
	x -= v.x;
	y -= v.y;
	z -= v.z;

	return *this;
}

inline Vec3 Vec3::operator *(const float a) const
{
	return Vec3(x * a, y * a, z * a);
}

inline Vec3& Vec3::operator *=(const float a)
{
	x *= a;
	y *= a;
	z *= a;

	return *this;
}

inline Vec3 operator *(const float a, const Vec3 &v)
{
	return v * a;
}

inline float Vec3::operator *(const Vec3 &v) const
{
	return (x * v.x + y * v.y + z * v.z);
}

inline Vec3 Vec3::operator /(const float a) const
{
	const float f = 1.0f / a;
	return Vec3(x * f, y * f, z * f);
}

inline Vec3& Vec3::operator /=(const float a)
{
	*this *= 1.0f / a;

	return *this;
}

inline float& Vec3::operator[](const int i)
{
	return (&x)[i];
}

inline float Vec3::operator[](const int i) const
{
	return (&x)[i];
}

inline float Vec3::length_sqr() const
{
	return (*this) * (*this);
}

inline float Vec3::length() const
{
	return sqrtf((*this) * (*this));
}

inline float Vec3::normalize()
{
	float l = this->length();
	*this /= l;

	return l;
}

inline void Vec3::abs()
{
	x = x < 0.0f ? -x : x;
	y = y < 0.0f ? -y : y;
	z = z < 0.0f ? -z : z;
}

inline void Vec3::minf(const Vec3 &v)
{
	x = x < v.x ? x : v.x;
	y = y < v.y ? y : v.y;
	z = z < v.z ? z : v.z;
}

inline void Vec3::maxf(const Vec3 &v)
{
	x = x > v.x ? x : v.x;
	y = y > v.y ? y : v.y;
	z = z > v.z ? z : v.z;
}

inline int Vec3::maxIndex()
{
	if (x > y)
	{
		if (x > z)
			return 0;
		else
			return 2;
	}
	else if (z > y)
		return 2;
	else
		return 1;
}

inline Vec3 Vec3::cross(const Vec3 &v, const Vec3 &w)
{
	Vec3 u;
	u.x = v.y * w.z - v.z * w.y;
	u.y = v.z * w.x - v.x * w.z;
	u.z = v.x * w.y - v.y * w.x;

	return u;
}

inline Vec3 Vec3::product(const Vec3 &v, const Vec3 &w)
{
	Vec3 u;
	u.x = v.x * w.x;
	u.y = v.y * w.y;
	u.z = v.z * w.z;

	return u;
}

inline void Vec3::onb(Vec3 &u, Vec3 &v, Vec3 &w, const Vec3 &n)
{
	const Vec3 m0(1.0f, 0.0f, 0.0f);
	const Vec3 m1(0.0f, 1.0f, 0.0f);

	w = n;
	w.normalize();

	u = Vec3::cross(w, m0);
	if (u.length_sqr() < 0.1f)
		u = Vec3::cross(w, m1);
	u.normalize();
	v = Vec3::cross(w, u);
	// w, u normalized and orthogonal -> v normalized
}

#endif
