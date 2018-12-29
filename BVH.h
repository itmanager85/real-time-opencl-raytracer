#pragma once

#include "BVH_Node.h"

#include "Mesh.h"

#include "vectors_math.h"

#include <vector>

using namespace std;


class BVH
{
public:

	BVH_Node * root;

	BVH(void);

	void build ( Mesh & mesh ) {
		root = new BVH_Node();
		root->init( mesh );
		root->build_node( mesh );
	} ; 

	~BVH(void);
};

