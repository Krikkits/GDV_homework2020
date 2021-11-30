/**
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#include "bvh.h"

BVH::BVH(const Triangle * const tris, const int nTris) :
		tris(tris), nTris(nTris)
{
	bbox = Triangle::getAABB(tris, nTris);

	nodes = new Node[nTris * 2]; // a bvh has at most 2 * n - 1 nodes

	indices = new int[nTris];
	for (int i = 0; i < nTris; i++)
		indices[i] = i;

	addedNodes = 1; // the root node
	buildBVH(0, 0, nTris, bbox, 0); // recursive construct BVH
}

BVH::~BVH()
{

}

void BVH::buildBVH(int nodeIndex, int triIndex, int numTris, AABB box,
		int depth)
{ // recursive construct BVH
	nodes[nodeIndex].bbox = box;
	nodes[nodeIndex].triIndex = -1;

	if (numTris <= 3 || depth > 63)
	{ // make leaf (normal termination)
		nodes[nodeIndex].left = -1;
		nodes[nodeIndex].right = -1;
		nodes[nodeIndex].triIndex = triIndex;
		nodes[nodeIndex].numTris = numTris;
	}
	else
	{
		float plane;
		int split;
		AABB leftBox, rightBox; // empty boxes (see AABB constructor in rtStructs.h)
		bool foundSplit = false;
		int axis = box.getMaxAxis();

		for (int i = 0; i < 3; i++)
		{ // if split failes in dim(axis), try to split along other dimensions
			plane = (box.bounds[0][axis] + box.bounds[1][axis]) * 0.5f;
			split = sortTris(triIndex, numTris, plane, axis, leftBox, rightBox);
			if (split != 0 && split != numTris)
			{
				foundSplit = true;
				break;
			}
			axis = (axis + 1) % 3;
		}
		if (!foundSplit)
		{ // make leaf (termination because no split could be found)
			nodes[nodeIndex].left = -1;
			nodes[nodeIndex].right = -1;
			nodes[nodeIndex].triIndex = triIndex;
			nodes[nodeIndex].numTris = numTris;
		}
		else
		{ // inner node
			int freeNode = addedNodes;
			addedNodes += 2;

			nodes[nodeIndex].left = freeNode;
			nodes[nodeIndex].right = freeNode + 1; // right = left+1 ... could be optimized

			buildBVH(freeNode, triIndex, split, leftBox, depth + 1); // recursively call buildBVH for left and right sub-tree
			buildBVH(freeNode + 1, triIndex + split, numTris - split, rightBox,
					depth + 1);
		}
	}
}

inline int BVH::sortTris(int triIndex, int numTris, float plane, int axis,
		AABB &leftBox, AABB &rightBox)
{
	int left = 0;
	int right = numTris - 1;

	while (left <= right)
	{ // with unsigned int this wount work! Modify this if you have more than 2^31 triangles ;-)
		AABB triBox = tris[indices[triIndex + left]].getAABB();
		const float tri_mid = (triBox.bounds[0][axis] + triBox.bounds[1][axis])
				* 0.5f;

		if (plane > tri_mid)
		{ // triangle is left
			left++;
			leftBox.extend(triBox);
		}
		else
		{ // triangle is right
			int temp = indices[triIndex + left];
			indices[triIndex + left] = indices[triIndex + right];
			indices[triIndex + right] = temp;
			right--;
			rightBox.extend(triBox);
		}
	}
	return left;
}

HitRec BVH::intersect(const Ray &ray) const
{
	// precompute inverse ray direction for bounding box intersection
	// -> compute only once for bvh travesal
	Vec3 invRayDir(1.0f / ray.dir.x, 1.0f / ray.dir.y, 1.0f / ray.dir.z);

	// index array for ray signs (direction of the ray)
	// -> saving costly if-operations during box intersection and traversal
	unsigned int raySign[3][2];
	raySign[0][0] = invRayDir[0] < 0;
	raySign[1][0] = invRayDir[1] < 0;
	raySign[2][0] = invRayDir[2] < 0;

	raySign[0][1] = invRayDir[0] >= 0;
	raySign[1][1] = invRayDir[1] >= 0;
	raySign[2][1] = invRayDir[2] >= 0;

	HitRec rec;
	float tmin = ray.tmin;
	float tmax = ray.tmax;

	if (!nodes->bbox.intersect(ray, tmin, tmax, invRayDir, raySign))
		return rec;

	TraversalStack stack[64]; // see bvh.h
	int nodeIndex = 0;
	int stackPos = 0;

	// iterative traversal ...
	while (true)
	{
		if (nodes[nodeIndex].triIndex == -1)
		{
			const int child0 = nodes[nodeIndex].left;
			const int child1 = nodes[nodeIndex].right;

			float tmin0 = tmin, tmin1 = tmin, tmax0 = tmax, tmax1 = tmax;
			const bool hit0 = nodes[child0].bbox.intersect(ray, tmin0, tmax0,
					invRayDir, raySign);
			const bool hit1 = nodes[child1].bbox.intersect(ray, tmin1, tmax1,
					invRayDir, raySign);

			// can be further optimized:
			// traverse near child first, either by comparing tmin0 and tmin1 or by
			// the ray direction (raySign) and the split axis (the split axis has to be saved in the node struct)

			if (hit0)
			{
				nodeIndex = child0;
				if (hit1)
				{
					stack[stackPos].tmin = tmin1;
					stack[stackPos].nodeIndex = child1;
					stackPos++;
				}
				continue;
			}
			else if (hit1)
			{
				nodeIndex = child1;
				continue;
			}
		}
		else
		{ // intersect triangles
			for (int i = nodes[nodeIndex].triIndex;
					i < nodes[nodeIndex].triIndex + nodes[nodeIndex].numTris;
					i++)
			{
				const int tri_id = indices[i];
				tris[tri_id].intersect(ray, rec, tri_id);
			}
		}
		while (true)
		{
			if (stackPos <= 0)
				return rec;
			if (stack[--stackPos].tmin <= rec.dist)
			{ // early termination: If we already found an intersection which is nearer than any
			  //                    other box on the stack
				nodeIndex = stack[stackPos].nodeIndex;
				tmax = minf(tmax, rec.dist);
				break;
			}
		}
	}
	return rec;
}

