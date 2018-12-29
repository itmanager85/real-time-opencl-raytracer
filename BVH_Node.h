#pragma once
#include "AABB.h"
#include "Mesh.h"

#include <vector>

using namespace std;

//const int NUM_DIVIDE = 10; // 20; // 
const int NUM_NODE_TRIS = 4; //2; //3;

const int MAX_DEPTH = 5; //10; //20; // 2; // 

//const int MIN_DIVIDE_PART = 2; // 1; 

const float SAH_COST = 0.0f;

class BVH_Node
{
public:
	AABB bounds;

	BVH_Node * pLeft;
	BVH_Node * pRight;

	bool isLeaf;

	vector <int> tris;

public:
	BVH_Node(void);

	void init ( Mesh & mesh );
	void find_bounds ( Mesh & mesh );

	void build_node ( Mesh & mesh, int depth = 0 ) ;
	float compute_sah();
	bool findSplit ( Mesh & mesh, int depth );
	//void find_nodes ( Mesh & mesh, int axis, int div_i, bool flag = false );
	void find_nodes ( Mesh & mesh, int axis, float testSplit, bool flag = false );
	void find_nodes2 ( Mesh & mesh, int axis, int div_i, int div_j );

	~BVH_Node(void);
};

