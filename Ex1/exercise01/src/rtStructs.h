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
{ // 32 Byte
/// Origin of the ray.
	Vec3 origin;
	/// Direction of the ray
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
	{
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
	 * Intersects a ray with the bounding box.
	 * @param r Ray to intersect.
	 * @param intervalMin In/out parameter: Contains the minimum of the
	 * interval to check before the call and the minimum of the interval
	 * in which the ray is inside the check-interval and inside the AABB
	 * after the call.
	 * @param intervalMax In/out parameter: Contains the minimum of the
	 * interval to check before the call and the minimum of the interval
	 * in which the ray is inside the check-interval and inside the AABB
	 * after the call.
	 * @returns True if the ray intersects with the AABB in the given
	 * interval.
	 */
	inline bool intersect(const Ray &r, float &intervalMin,
			float &intervalMax) const
	{
		//compute t values for x axis
		float tmin = (intervalMin- - r.origin.x) / r.dir.x;
		float tmax = (intervalMax - r.origin.x) / r.dir.x;

		if (tmin > tmax){
			std::swap(tmin, tmax);
		}

		//compute t values for y axis
		float tymin = (intervalMin - r.origin.y) / r.dir.y;
		float tymax = (intervalMax - r.origin.y) / r.dir.y;

		if (tymin > tymax){
			std::swap(tymin, tymax);
		}

		if ((tmin > tymax) || (tymin > tmax)){ //miss
			return false;
		}

		if (tymin > tmin){
			tmin = tymin;
		}

		if (tymax < tmax){
			tmax = tymax;
		}

		//compute t values for z axis
		float tzmin = (intervalMin - r.origin.z) / r.dir.z;
		float tzmax = (intervalMax - r.origin.z) / r.dir.z;

		if (tzmin > tzmax){
			std::swap(tzmin, tzmax);
		}
		if ((tmin > tzmax) || (tzmin > tmax)){ //also miss
			return false;
		}

		if (tzmin > tmin){
			tmin = tzmin;
		}
		if (tzmax < tmax){
			tmax = tzmax;
		}

		return true;
	}
};

/**
 * Triangle.
 */
struct Triangle
{
	/// Vertices of the triangle (unordered).
	Vec3 v[3];

	/**
	 * Returns the AABB of the triangle.
	 * @returns AABB containing all vertices of the triangle.
	 */
	inline AABB getAABB() const
	{
		AABB bbox;
		Vec3 min=v[0]; //min(minx,miny,minz)
		Vec3 max=v[0];

		// find greatest x y and z (for max) as well as smallest (for min)
		for(int i=0; i<3;i++){ //min and max coordinates
			for(int j=0;j<3;j++){
				if(v[j][i]<min[i]){
					min[i]=v[j][i];
				}
				if(v[j][i]<max[i]){
					max[i]=v[j][i];
				}
			}
		}

		bbox.bounds[0]=min;
		bbox.bounds[1]=max;

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

		AABB bbox;
		AABB trisBbox[nTris];

		//idea: calculate the AABB for individual triangles and save it in a new array
		//then find the smallest and biggest x, y and z value
		//calculate bounding box out of that



		//Calculate a bounding box for every triangle and save it in an Array
		for ( int i = 0; i < nTris; i++){
			trisBbox[i] = tris[i].getAABB();
		}

		//not sure what to initialize these as, hopefully this gets me the first triangle coordinates from tris
		Vec3 min = Vec3(trisBbox[0].bounds[0][0],trisBbox[1].bounds[0][0],trisBbox[2].bounds[0][0]);
		Vec3 max = Vec3(trisBbox[0].bounds[1][0],trisBbox[1].bounds[1][0],trisBbox[2].bounds[1][0]);

		//Find the smallest and biggest min/max from the bounding boxes
		//dependent on the axis
		for(int i=0;i<nTris; i++){
			for(int j=0;j<3; j++){
				if (trisBbox[i].bounds[0][j] < min[j]){
					min[j] = trisBbox[i].bounds[0][j];
				}
				if (trisBbox[i].bounds[1][j] > max[j]){
					max[j] = trisBbox[i].bounds[1][j];
				}
			}
		}



		bbox.bounds[0] = min;
		bbox.bounds[1] = max; //a bit off from intended? not sure why??

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
	 * @param True if the ray intersects with the triangle in the interval
	 * set in the ray.
	 */
	inline bool intersect(const Ray &ray, HitRec &rec, const int tri_id) const
	{
		//based off of what we learned in the lecture
		//feedback appreciated!
		float kEpsilon = 1e-8;
		Vec3 b=v[1]-v[0];
		Vec3 c=v[2]-v[0];
		Vec3 N=Vec3::cross(b,c);
		float D=v[0]*N;
		float t=(D-ray.origin*N)/(ray.dir*N);

		//miss if parallel
		if(t<kEpsilon||t>ray.tmax){
			return false;
		}

		//determine projections
		int k;
		if(fabs(N[0])>fabs(N[1])){
			if(fabs(N[0])>fabs(N[2])){
				k=0; //project on x
			} else{
				k=2; //project on z
			}
		}else{
			if(fabs(N[1])>fabs(N[2])){
				k=1; //project on y
			}else{
				k=2; //project on z
			}
		}

		int x=(k+1)%3;
		int y=(k+2)%3;

		//hitpoints
		float p[3];
		Vec3 P;

		p[x]=ray.origin[x]+t*ray.dir[x]-v[0][x];
		p[y]=ray.origin[y]+t*ray.dir[y]-v[0][y];

		//check if the points are inside triangle or not
		float b1=b[x]*p[y]-b[y]*p[x];
		float b2=b[x]*c[y]-b[y]*c[x];
		float beta=b1/b2;
		float gamma=(c[y]*p[x]-c[x]*p[y])/(b[x]*c[y]-b[y]*c[x]);
		P=beta*b+gamma*c;
		if((beta+gamma)>1){
			return false;
		}
/* redundant
		if(beta<0){
			return false;
		}

		if(gamma<0){
			return false;
		}
*/
		if(beta>0&&gamma>0){ //if it's a hit we run the loop

				Vec3 dirToVertex = P - ray.origin;
				float distToVertex = dirToVertex.length();
				dirToVertex.normalize();
				if (dirToVertex * ray.dir > kEpsilon) // Vertex approximately in ray direction?
					{
					if (distToVertex < rec.dist) // Distance shorter than previous?
						{
						rec.id = tri_id;
						rec.dist = distToVertex;
						return true;
						}
					}

		}


		return false;

		/* original dummy code:
		for (unsigned int i = 0; i < 3; i++)
				{
					Vec3 dirToVertex = v[i] - ray.origin;
					float distToVertex = dirToVertex.length();
					dirToVertex.normalize();
					if (dirToVertex * ray.dir > 0.9999f) // Vertex approximately in ray direction?
					{
						if (distToVertex < rec.dist) // Distance shorter than previous?
						{
							rec.id = tri_id;
							rec.dist = distToVertex;
							return true;
						}
					}
				}
				return false;*/

	}

	/**
	 * Intersects a ray with the triangle.
	 * @param ray Ray to intersect with the triangle.
	 * @returns True if the ray intersects with the triangle.
	 */
	inline bool intersectShadow(const Ray &ray) const
	{

		//slightly different approach based off my own research

		//same starting values as before
		Vec3 b=v[1]-v[0];
		Vec3 c=v[2]-v[0];
		Vec3 N=Vec3::cross(b,c);
		float length=N.length();

		//check if ray and plane are parallel (a miss)
		float rayDirxN=N*ray.dir;
		if(fabs(rayDirxN)<0.9999f){
			return false;
		}

		//calculate d
		float d=N*v[0];
		//calculate t
		float t=(N*ray.origin + d)/rayDirxN;
		//miss if triangle is behind the ray
		if(t<0){
			return false;
		}

		//calculate intersection point
		Vec3 P=ray.origin+t*ray.dir;

		//check if point is inside or outside of triangle
		Vec3 C;

		Vec3 edge0=v[1]-v[0];
		Vec3 v0=P-v[0];
		C=Vec3::cross(edge0,v0);
		if((N*C)<0){
			return false;
		}

		Vec3 edge1=v[2]-v[1];
		Vec3 v1=P-v[1];
		C=Vec3::cross(edge1,v1);
		if((N*C)<0){
			return false;
		}

		Vec3 edge2=v[0]-v[2];
		Vec3 v2=P-v[2];
		C=Vec3::cross(edge2,v2);
		if((N*C)<0){
			return false;
		}


		return true;
	}

	/**
	 * Returns the normal of the triangle.
	 * @remarks It is undefined if the normal points to the front or the back
	 * side of the triangle because the vertex order is undefined.
	 */
	inline Vec3 getNormal() const
	{
		Vec3 N;
		Vec3 b=v[1]-v[0];
		Vec3 c=v[2]-v[0];
		N=Vec3::cross(b,c);

		return N;
	}

};

#endif

