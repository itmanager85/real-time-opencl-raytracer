#pragma once
#include "BVH.h"

#include "BVH2.h"

#include <vector>

using namespace std;	

//using namespace FW;	

struct BVH_Node_ {

	AABB aabb;

	int offset_left;
	int offset_right;

	int offset_tris;
	int num_tris;

	BVH_Node_() {
		offset_left = -1;
		offset_right = -1;

		offset_tris = -1;
		num_tris = 0;
	};
} ; 

//struct BVH_Node_Tris {
//float3 v0, v1, v2;
//} ; 

class BVH_Cuda
{
	struct SP {

	} ; 

public:
	vector <BVH_Node_> bvh_nodes;	// 
	vector <int> tri_indices;		// 

public:
	BVH_Cuda(void);
	~BVH_Cuda(void);

	// return offset node in bvh_nodes
	int build ( BVH_Node * node ) { // , int & offLeft, int & offRight, int

		static int num_bvh_nodes = 0;
		static int tri_index = 0;

		int curr_bvh_node_index = num_bvh_nodes;

		BVH_Node_ bvh_node;
		bvh_nodes.push_back( bvh_node );

		bvh_node.aabb = node->bounds;

		if ( !node->isLeaf ) {

			++num_bvh_nodes;
			bvh_node.offset_left = build ( node->pLeft );

			++num_bvh_nodes;
			bvh_node.offset_right = build ( node->pRight );

		} else {
			//tri_indices.resize( tri_indices.size() + node->tris.size() );

			for ( int i=0; i<node->tris.size(); ++i )
				tri_indices.push_back( node->tris[i] );

			bvh_node.offset_tris = tri_index;
			bvh_node.num_tris = node->tris.size();

			tri_index += node->tris.size();
		} 

		bvh_nodes[ curr_bvh_node_index ] = bvh_node;

		return curr_bvh_node_index;
	} ; 

	void build_from_bvh2 ( FW::BVH2 & bvh2 ) {

		//tri_indices.resize( bvh2.m_triIndices.size() );
		tri_indices = bvh2.m_triIndices;

		for ( int i=0; i<tri_indices.size(); ++i )
			tri_indices[i] *= 3;

		build2( bvh2.getRoot() );
	} ; 

int build2 ( FW::BVHNode * node ) { // , int & offLeft, int & offRight, int

		static int num_bvh_nodes = 0;
		static int tri_index = 0;

		int curr_bvh_node_index = num_bvh_nodes;

		BVH_Node_ bvh_node;
		bvh_nodes.push_back( bvh_node );

		//bvh_node.aabb = node->m_bounds;
		bvh_node.aabb.set( node->m_bounds.minf(), node->m_bounds.maxf() );

		if ( !node->isLeaf() ) {

			++num_bvh_nodes;
			bvh_node.offset_left = build2 ( node->getChildNode(0) ); // pLeft

			++num_bvh_nodes;
			bvh_node.offset_right = build2 ( node->getChildNode(1) ); // pRight

		} else {
			//tri_indices.resize( tri_indices.size() + node->tris.size() );

			//for ( int i=0; i<node->tris.size(); ++i )
			//	tri_indices.push_back( node->tris[i] );

			FW::LeafNode * pLeaf = (FW::LeafNode*)node;

			bvh_node.offset_tris = pLeaf->m_lo; // tri_index;
				
			bvh_node.num_tris = pLeaf->getNumTriangles(); // tris.size();

			//tri_index += node->tris.size();
		} 

		bvh_nodes[ curr_bvh_node_index ] = bvh_node;

		return curr_bvh_node_index;
	} ; 
};

