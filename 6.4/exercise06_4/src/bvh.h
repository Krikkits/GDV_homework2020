#ifndef BVH_H
#define BVH_H

#include "rtStructs.h"

/**
 * One node of the BVH tree.
 * @remarks a BVH tree node has either child nodes or triangles, never both!
 */
struct Node
{
	/// Includes either all its child nodes or all its vertices.
	AABB bbox;
	/**
	 * Index of the left child node in the BVH::nodes array.
	 * @remarks Should be -1 in leaf nodes.
	 */
	int left;
	/**
	 * Index of the right child node in the BVH::nodes array.
	 * @remarks Should be -1 in leaf nodes.
	 */
	int right;

	/**
	 * A Node contains all the triangles referenced by BVH::indices' elements
	 * triIndex to triIndex + numTris - 1.
	 * @remarks Should be -1 if not a leaf node.
	 */
	int triIndex;
	/**
	 * Number of triangles of this node.
	 * @remarks Should be 0 if not a leaf node.
	 */
	int numTris;
};

struct BVH
{
	/// Bounding box containing the whole scene that is accelerated by the BVH.
	AABB bbox;
	/// Array of all triangles indexed by the BVH (access indirectly via indices!)
	const Triangle * const tris;
	/// Number of triangles.
	const int nTris;

	// this is the bvh!
	/// Array containing all nodes of the BVH.
	Node *nodes;
	/// Amount of nodes that have already been added (is incremented while building the tree).
	int addedNodes;
	/// Used for indirectly accessing the tris. indices[i] contains the index of a triangle in tris.
	int *indices;

	/**
	 * Builds a BVH tree for a set of triangles.
	 * @param tris Array of triangles.
	 * @param nTris Number if triangles in tris.
	 */
	BVH(const Triangle * const tris, const int nTris);
	/**
	 * Default deconstructor.
	 */
	~BVH();

	/**
	 * Finds a good split plane.
	 * @param dimension: Out parameter: Returns the splitted dimension.
	 * @param position: Out parameter: Returns the split position in dimension.
	 * This value should be exclusive in the lower and inclusive in the higher
	 * node.
	 * @param Index of the minimum indices array element that refers to a
	 * triangle belonging to the set that has to be split.
	 * @param Number of triangles that define the volume to be split.
	 * @returns True if more than one triangle -> split makes sense.
	 */
	bool findSplitPlane(unsigned int &dimension, float &position,
			const unsigned int triIndex, const unsigned int numTris);

	/**
	 * Builds a BVH tree out of a set of triangles in tris.
	 * @remarks: Intended to call itself recursively for each subnode.
	 * @param nodeIndex Root node index of the subtree that is processed
	 * during the call.
	 * @param triIndex First element of indices referencing a triangle that
	 * belongs to the current node.
	 * @param numTris Number of triangles belonging to the current node.
	 * @param box Bounding box containing all the current node's triangles.
	 * @param depth Recursion depth = tree depth of the node currently
	 * processed.
	 */
	void buildBVH(int nodeIndex, int triIndex, int numTris, AABB box,
			int depth);

	/**
	 * Finds the closest intersection of a ray and all the triangles in the
	 * BVH efficiently.
	 * @ray Ray to intersect with the BVH's triangles.
	 * @returns Hit record of the closest intersection between triangle and ray.
	 */
	HitRec intersect(const Ray &ray) const;

	// add extra methods if needed (i.e. sortTriangles(...))
};

#endif

