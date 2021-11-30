/**
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

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

/**
 * Contains all that would be importent in a stack frame of a recursive
 * traversal algorithm. Using this own stack frame, traversal can be done
 * iteratively and faster.
 */
struct TraversalStack
{
	/// Index of the node that must be traversed in this stack frame.
	int nodeIndex;
	/// Closest possible intersection distance.
	float tmin;
};

/**
 * Bounded Volume Hierarchy
 */
struct BVH
{
	/// Bounding box containing the whole scene that is accelerated by the BVH.
	AABB bbox;
	/// Array of all triangles indexed by the BVH (access indirectly via indices!)
	const Triangle * const tris;
	/// Number of triangles.
	const int nTris;

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
	 * Sorts an array of triangles based on a split plane into two segments and
	 * calculates their bounding boxes.
	 * @param triIndex Index of the first triangle in indices to concern.
	 * @param numTris Number of triangles following on triIndex in indices to
	 * concern.
	 * @param plane Split plane position.
	 * @param axis Coordinate axis on which the space is split at plane.
	 * @param leftBox Out paramter: Bounding box around all triangles
	 * before the split plane.
	 * @param rightBox Out parameter: Bounding box around all triangles
	 * behind the split plane.
	 * @returns First triangle (element of indices) that is behind the
	 * split plane after sorting.
	 */
	inline int sortTris(int triIndex, int numTris, float plane, int axis,
			AABB &leftBox, AABB &rightBox);

	/**
	 * Finds the closest intersection of a ray and all the triangles in the
	 * BVH efficiently.
	 * @ray Ray to intersect with the BVH's triangles.
	 * @returns Hit record of the closest intersection between triangle and ray.
	 */
	HitRec intersect(const Ray &ray) const;
};

#endif

