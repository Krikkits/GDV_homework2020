#include "bvh.h"
#include <assert.h>
#include <iostream>

using namespace std;

BVH::BVH(const Triangle * const tris, const int nTris) :
		tris(tris), nTris(nTris)
{
	nodes = new Node[nTris * 2]; // a bvh has at most 2 * n - 1 nodes

	indices = new int[nTris];
	for (int i = 0; i < nTris; i++)
		indices[i] = i;

	addedNodes = 1; // the root node
	buildBVH(0, 0, nTris, 0); // recursive construct BVH
}

BVH::~BVH()
{
}

void BVH::findSplitPlane(unsigned int* dimension, float* position, const AABB &box)
{
	// TODO b) Implement the method.
	//0=x , 1=y , 2=z (Axis)
	//ghetto debugging lines inbetween (uncommented)
	//int debug=100;
	unsigned int gay=box.getMaxAxis();
	//std::cout<<"Gay: "<<gay<<endl;
	*dimension=gay; //find the longest axis to split
	//std::cout<<"Dimension: "<<dimension<<endl
	*position=(box.bounds[1][*dimension] - box.bounds[0][*dimension])/2;


}

void BVH::buildBVH(int nodeIndex, int triIndex, int numTris, int depth)
{
	// TODO c) Implement the method by replacing the code below.
	// Rough Steps:
	// - Find the split plane using the method from b)
	// - Sort all triangles belonging to the currently processed node depending on their center's location relative to the split plane.
	// - Recursively call this method to process subnodes.
	// Don't forget about the base case (no split needed)!

	int debug=0;
	unsigned int dim=0;
	unsigned int* dimension=&dim;
	float pos=0;
	float* position=&pos;

	//std::cout<<"Debug: "<<debug++<<endl;
	//root node
	if(nodeIndex==0){
		nodes[nodeIndex].bbox = Triangle::getAABB(tris, numTris);
		nodes[nodeIndex].left=nodeIndex+1;
		nodes[nodeIndex].right=nodeIndex+1;
		BVH::addedNodes=addedNodes;
		nodes[nodeIndex].numTris = 0;
		nodes[nodeIndex].triIndex = -1;
	}
	//std::cout<<"Debug: "<<debug++<<endl;
	BVH::findSplitPlane(dimension, position, nodes[nodeIndex].bbox); //done at every iteration until leaf
	//std::cout<<"Debug: "<<debug++<<endl;
	//calculate next BBox
	Vec3 posVec; //also pivot for sorting (min+position vector?)
	if(*dimension==0){ //if x axis
		posVec = Vec3(*position, 0, 0);
	}
	if(*dimension==1){ //if y axis
		posVec=Vec3(0, *position, 0);
	}
	if(*dimension==2){ //if z axis
		posVec=Vec3(0, 0, *position);
	}
	//std::cout<<"Debug: "<<debug++<<endl;
	//left and right BBox
	AABB leftBBox=AABB(nodes[nodeIndex].bbox.bounds[0], nodes[nodeIndex].bbox.bounds[1]-posVec);
	AABB rightBBox=AABB(nodes[nodeIndex].bbox.bounds[0]-posVec, nodes[nodeIndex].bbox.bounds[1]);
	//std::cout<<"Debug: "<<debug++<<endl;

	if(numTris>2){ //condition to keep calling recursively

		//can't quite figure out how to implement the sort properly
		/*sorting
		float start=nodes[nodeIndex].bbox.bounds[0][dimension];
		float end=nodes[nodeIndex].bbox.bounds[1][dimension];
		float pivot = end; // pivot
		int i = (start - 1); // Index of smaller element

		for (int j = start; j <= end - 1; j++)
		{
			// If current element is smaller than the pivot
			if (BVH::tris[BVH::indices[j][dimension]] < pivot)
			{
				i++; // increment index of smaller element
				std::swap(BVH::tris[BVH::indices[i][dimension]], BVH::tris[BVH::indices[j][dimension]]);
			}
		}
		std::swap(&BVH::tris[BVH::indices[1][dimension]], end); //swap(&arr[i + 1], &arr[high]);
		i++;


		//bvh->tris[bvh->indices[node->triIndex]]

		//Tricompare (const void* a, const void* b)
		//Vec3 diff=a.getAABB().getCenter()-b.getAABB().getCenter();
		//if(diff[*dimension]<0)return -1
		//if(diff[*dimension]==0)return 0
		//f(diff[*dimension]>0)return 1
		*/

		//find out how many triangles are left and right
		int leftNumTris=0;
		int rightNumTris=0;
		for(int i=0;i<numTris-1;i++){
			Vec3 diff=BVH::tris[BVH::indices[i]].getAABB().getCenter() - nodes[nodeIndex].bbox.getCenter();
			if(diff[*dimension]<0){
				leftNumTris++;
			}
			if(diff[*dimension]>=0){
				rightNumTris++;
			}
		}

		//left node

		nodes[nodeIndex].bbox=leftBBox;
		nodes[nodeIndex].left=nodeIndex+1;
		BVH::addedNodes=BVH::addedNodes+1;
		nodes[nodeIndex].numTris = 0;
		nodes[nodeIndex].triIndex = -1;
		buildBVH(nodes[nodeIndex].left,1, leftNumTris, depth++);



		//right node

		nodes[nodeIndex].bbox=rightBBox;
		nodes[nodeIndex].right=nodeIndex+1;
		BVH::addedNodes=BVH::addedNodes+1;
		nodes[nodeIndex].numTris = 0;
		nodes[nodeIndex].triIndex = -1;
		buildBVH(nodes[nodeIndex].right,1, rightNumTris, depth++);
	}


	//leaf nodes here:
	nodes[nodeIndex].bbox=leftBBox;
	nodes[nodeIndex].left=-1;
	nodes[nodeIndex].right=-1;
	BVH::addedNodes=BVH::addedNodes+1;
	nodes[nodeIndex].numTris = numTris;
	nodes[nodeIndex].triIndex = triIndex;

	nodes[nodeIndex].bbox=rightBBox;
	nodes[nodeIndex].right=-1;
	nodes[nodeIndex].left=-1;
	BVH::addedNodes=BVH::addedNodes+1;
	nodes[nodeIndex].numTris = numTris;
	nodes[nodeIndex].triIndex = triIndex;


	//recursive call: buildBVH(nodeIndex+1, triIndex+1, numTris, depth)




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

	// TODO d) instead of intersecting all triangles, traverse the bvh
	// and intersect only the triangles in the leave nodes.
	// Attention: Make sure that hit records and tmin/tmax values are not modified by node intersections!

	float tmin = ray.tmin;
	float tmax = ray.tmax;
	int nodeIndex=0;

	//miss
	if (!nodes->bbox.intersect(ray, tmin, tmax, invRayDir, raySign)){
		return rec;
	}

	//(attempted) recursive method

	if(nodes[nodeIndex].triIndex==-1){
		for (int i = nodes[nodeIndex].triIndex; i < nodes[nodeIndex].triIndex + nodes[nodeIndex].numTris; i++){
			const int tri_id = indices[i];
			tris[tri_id].intersect(ray, rec, tri_id);
			nodeIndex++;
		}
		bool ifHit=nodes[nodeIndex].bbox.intersect(ray, tmin, tmax, invRayDir, raySign); //check if bbox intersects with ray
		if(ifHit==true){
			//recurse through left and right nodes in the same way until leaf is found
			//!!not sure how to do recursive call with this:
			nodeIndex++;
			//this->nodes[nodes->left].intersect(ray);
			//this->nodes[nodes->right].intersect(ray);
		}

	}



//dummy code left in because our intersect isn't correct
		if (nodes[0].triIndex != -1)
		{
			for (int i = nodes[0].triIndex;
					i < nodes[0].triIndex + nodes[0].numTris; i++)
			{
				tris[indices[i]].intersect(ray, rec, i);
			}
		}


	return rec;
}

