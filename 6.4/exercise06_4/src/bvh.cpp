#include "bvh.h"
#include <vector>
#include <assert.h>
#include <iostream>

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

bool BVH::findSplitPlane(unsigned int &dimension, float &position,
		const unsigned int triIndex, const unsigned int numTris)
{
	if (numTris <= 1)
		return false;

	AABB volume;
	for (unsigned int t = triIndex; t < triIndex + numTris; t++)
		volume.extend(tris[indices[t]].getAABB());

	dimension = volume.getMaxAxis();
	position = volume.getCenter()[dimension];

	return true;
}

void BVH::buildBVH(int nodeIndex, int triIndex, int numTris, AABB box,
		int depth)
{
	nodes[nodeIndex].bbox = box;

	// Try to split
	bool splitDone = false;
	unsigned int splitDimension;
	float splitPosition;
	if (findSplitPlane(splitDimension, splitPosition, triIndex, numTris))
	{
		// Sort...
		int start = triIndex;
		int end = triIndex + numTris - 1;
		for (; start < end; start++)
		{
			const Triangle &startTriangle = tris[indices[start]];
			if (startTriangle.getAABB().getCenter()[splitDimension]
					>= splitPosition)
			{
				while (tris[indices[end]].getAABB().getCenter()[splitDimension]
						>= splitPosition && end > start)
					end--;
				if (end > start)
				{
					int tmp = indices[start];
					indices[start] = indices[end];
					indices[end] = tmp;
				}
				else
				{
					break; // Don't increment start anymore...
				}
			}
		}

		// start holds the beginning of the second subtree here...
		if (start != triIndex && start != triIndex + numTris)
		{
			// Both subtrees have leaves
			nodes[nodeIndex].triIndex = -1;
			nodes[nodeIndex].numTris = 0;

			// Add the left node
			nodes[nodeIndex].left = addedNodes;
			addedNodes++;
			AABB leftAABB;
			for (int t = triIndex; t < start; t++)
				leftAABB.extend(tris[indices[t]].getAABB());
			buildBVH(addedNodes - 1, triIndex, start - triIndex, leftAABB,
					depth + 1);

			// Add the right node
			nodes[nodeIndex].right = addedNodes;
			addedNodes++;
			AABB rightAABB;
			for (int t = start; t < triIndex + numTris; t++)
				rightAABB.extend(tris[indices[t]].getAABB());
			buildBVH(addedNodes - 1, start, triIndex + numTris - start,
					rightAABB, depth + 1);

			splitDone = true;
		}
	}

	if (!splitDone)
	{
		// If no split, make a leaf node.
		nodes[nodeIndex].left = -1;
		nodes[nodeIndex].right = -1;
		nodes[nodeIndex].triIndex = triIndex;
		nodes[nodeIndex].numTris = numTris;
	}
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

	if (!bbox.intersect(ray, tmin, tmax, invRayDir, raySign))
		return rec;

	std::vector<int> todo;
	todo.reserve(nTris);
	todo.push_back(0);

	while (todo.size() > 0)
	{
		int currentNodeIdx = todo.back();
		todo.pop_back();
		Node &currentNode = nodes[currentNodeIdx];

		float tmpTmin = tmin;
		float tmpTmax = tmax;
		if (currentNode.bbox.intersect(ray, tmpTmin, tmpTmax, invRayDir,
				raySign))
		{
			if (currentNode.left != -1)
			{
				todo.push_back(currentNode.right);
				todo.push_back(currentNode.left);
			}
			else
			{
				for (int i = currentNode.triIndex;
						i < currentNode.triIndex + currentNode.numTris; i++)
				{
					tris[indices[i]].intersect(ray, rec, indices[i]);
				}
			}
		}
	}

	return rec;
}

