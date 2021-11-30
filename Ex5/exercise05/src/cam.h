/**
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#ifndef CAM_H
#define CAM_H

#include "utils/vec.h"
#include "rtStructs.h"

struct Cam
{
	Vec3 eye, u, v, w;
	Vec3 corner, across, up;

	int ResX, ResY;
	float dist;
	float rX, rY;
	float ratio;

	float u0, u1, v0, v1;

	float speed;
	int action;
	float rotx, roty;
	bool moved;

	Cam(const AABB &box, const int ResX, const int ResY) :
			ResX(ResX), ResY(ResY)
	{
		Vec3 diff = box.bounds[1] - box.bounds[0];
		Vec3 middle = (box.bounds[0] + box.bounds[1]) * 0.5f;
		eye = middle;
		eye[2] += diff.length();

		w = Vec3(0.0f, 0.0f, -1.0f);
		v = Vec3(0.0f, 1.0f, 0.0f);
		u = Vec3::cross(w, v);

		rX = 1.0f / (float) ResX;
		rY = 1.0f / (float) ResY;
		ratio = ResY * rX;

		u0 = -1.0f;
		u1 = 1.0f;
		v0 = -ratio;
		v1 = ratio;
		dist = 2.5f;

		across = (u1 - u0) * u;
		up = (v1 - v0) * v;
		corner = u0 * u + v0 * v + dist * w;

		speed = 0.1f * diff.length();
		action = 0;
		rotx = 0.0f;
		roty = 0.0f;
		moved = false;
	}

	Cam(const Vec3 &eye, const Vec3 &w_, const Vec3 &v_, const float speed,
			const int ResX, const int ResY) :
			eye(eye), v(v_), w(w_), ResX(ResX), ResY(ResY), speed(speed)
	{
		w.normalize();
		v.normalize();
		u = Vec3::cross(w, v);
		u.normalize();
		rX = 1.0f / (float) ResX;
		rY = 1.0f / (float) ResY;
		ratio = ResY * rX;

		u0 = -1.0f;
		u1 = 1.0f;
		v0 = -ratio;
		v1 = ratio;
		dist = 2.5f;

		across = (u1 - u0) * u;
		up = (v1 - v0) * v;
		corner = u0 * u + v0 * v + dist * w;

		this->speed *= 30.0;
		action = 0;
		rotx = 0.0f;
		roty = 0.0f;
		moved = false;
	}
	Ray getRay(const float x, const float y) const;
	void cam_move();
	void cam_w(bool b);
	void cam_a(bool b);
	void cam_s(bool b);
	void cam_d(bool b);
	void cam_rx(float rx);
	void cam_ry(float ry);
};

inline Ray Cam::getRay(const float x, const float y) const
{
	Vec3 dir = corner + (x * rX) * across + (y * rY) * up;
	dir.normalize();

	return Ray(eye, dir, 0.0f, RAY_MAX);
}

inline void Cam::cam_move()
{
	if (!(action || rotx != 0.0f || roty != 0.0f))
	{
		moved = false;
		return;
	}
	moved = true;

	float t = roty * M_PI / 180.0f;
	float p = rotx * M_PI / 180.0f;

	w += u * t + v * p;
	w.normalize();

	u = Vec3::cross(w, Vec3(0.0f, 1.0f, 0.0f));
	u.normalize();

	v = Vec3::cross(u, w);
	v.normalize();

	if (action & 1)
	{
		eye += speed * w;
	}
	if (action & 2)
	{
		eye -= speed * u;
	}
	if (action & 4)
	{
		eye -= speed * w;
	}
	if (action & 8)
	{
		eye += speed * u;
	}
	if (action & 16)
	{

	}
	across = (u1 - u0) * u;
	up = (v1 - v0) * v;

	corner = u0 * u + v0 * v + dist * w;

	rotx = 0.0f;
	roty = 0.0f;
}

inline void Cam::cam_w(bool b)
{
	action = b ? action | 1 : action & ~1;
}

inline void Cam::cam_a(bool b)
{
	action = b ? action | 2 : action & ~2;
}

inline void Cam::cam_s(bool b)
{
	action = b ? action | 4 : action & ~4;
}

inline void Cam::cam_d(bool b)
{
	action = b ? action | 8 : action & ~8;
}

inline void Cam::cam_rx(float rx)
{
	rotx -= rx;
}

inline void Cam::cam_ry(float ry)
{
	roty += ry;
}

#endif

