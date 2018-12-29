#include "BVH_Node.h"

#include <float.h>

BVH_Node::BVH_Node(void)
{
	pLeft = NULL;
	pRight = NULL;

	isLeaf = false ; 
}

	void BVH_Node::init ( Mesh & mesh ) {

		tris.resize( mesh.indices.size()/3 );
		for ( int i=0; i<mesh.indices.size()/3; ++i ) {
			tris[i] = i*3;
		}

		find_bounds ( mesh );
	};


	void BVH_Node::find_bounds ( Mesh & mesh ) {

		bounds.clear();
		for ( int i=0; i<tris.size(); ++i ) {

			int & tri = tris[i];
			bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );

			//printf ( "\n tris[%d] = ( %.2f, %.2f, %.2f ) ", i, tris[i].x, tris[i].y, tris[i].z );
		}
	} ; 

	void BVH_Node::build_node ( Mesh & mesh, int depth ) {

		bool flag_print_info = false; // true; // 

		if ( flag_print_info ) {
			printf ( "\n depth = %d, tris.size() = %d , tris[] by tri = ", depth, tris.size() );

			for ( int i=0; i<tris.size(); ++i)
				printf ( "%d, ", tris[i]/3 );

			bounds.print();
			getchar();
		}

		if ( depth >= MAX_DEPTH || tris.size() <= NUM_NODE_TRIS ) { // 10 // 5 //  || depth >= 5
			isLeaf = true;

			if ( flag_print_info ) {
				printf ( "\n isLeaf = true " );
			}
			return;
		} ; 

		pLeft  = new BVH_Node();
		pRight = new BVH_Node();

		if ( !findSplit( mesh, depth ) ) {
			delete pLeft;
			delete pRight;

			pLeft = NULL; 
			pRight = NULL; 

			isLeaf = true;

			if ( flag_print_info ) {
				printf ( "\n ggg!!!" );
			}

			return ; 
		} ; 

		 //printf ( "\n (2) depth = %d, tris.size() = %d ", depth, tris.size() );
		 //getchar();

		if ( flag_print_info ) {
			printf ( "\n depth = %d, pLeft : ", depth );
		}
		pLeft->build_node( mesh, ++depth );

		if ( flag_print_info ) {
			printf ( "\n depth = %d, pRight : ", (depth-1) );
		}
		pRight->build_node( mesh, depth );

		tris.clear();

		 //printf ( "\n depth = %d, tris.clear() ", depth );
		 //getchar();
	} ; 

	bool BVH_Node::findSplit ( Mesh & mesh, int depth ) {
		
		//float best_sah = exp(80.0f);
		int best_axis = -1, bestSplit; //best_div_i = -1; 
		//int best_div_j = 0;

		float minCost = tris.size() * bounds.getVolume();

		for ( int axis = 0; axis <=2; ++axis ) {
		//for ( int axis = 1; axis <=1; ++axis ) {

			float start = bounds.min.m[axis];
			float stop = bounds.max.m[axis];

			float step = (stop-start)/(32.f/(depth+1.f)); // 1024.f

			for(float testSplit=start+step; testSplit<stop-0.001f; testSplit+=step) {

			//int num_parts_max = (bounds.max.m[axis] - bounds.min.m[axis]) / MIN_DIVIDE_PART;
			//int num_parts = num_parts_max < NUM_DIVIDE ? num_parts_max : NUM_DIVIDE;

			//for ( int i=1; i<num_parts; ++i ) { // NUM_DIVIDE
			//for ( int i=1; i<=1; ++i ) {
				
				//for ( int j=0; j<=i; ++j ) {

					//find_nodes( mesh, axis, i );
					find_nodes( mesh, axis, testSplit );
					//find_nodes2( mesh, axis, i, j );

					//const float traversalCost = 0.0f; // 10.0f; //
					float sah = compute_sah();
					if ( sah < minCost ) {
						minCost = sah;
						best_axis = axis, 
						bestSplit = testSplit;
						//best_div_i = i;
						//best_div_j = j;
					}
				//}
			}
		}

		//printf( "\n best_div_i = %d ", best_div_i );

		if ( best_axis < 0 ) return false; 		

		//printf( "\n best_sah = %f ", best_sah );

		//if ( minCost <= SAH_COST ) return false ; 

		//if ( pLeft->tris.size() == tris.size() && pRight->tris.size() == tris.size() ) return false ; 

		//find_nodes( mesh, best_axis, best_div_i, true );
		find_nodes( mesh, best_axis, bestSplit, true );
		//find_nodes2( mesh, best_axis, best_div_i, best_div_j );

		return true;
	};

	float BVH_Node::compute_sah() {
		//return SAH_COST + pLeft->tris.size() * pLeft->bounds.getVolume() / bounds.getVolume()  + pRight->tris.size() * pRight->bounds.getVolume() / bounds.getVolume();
		//return pLeft->tris.size() * pLeft->bounds.getSurfaceArea()  + pRight->tris.size() * pRight->bounds.getSurfaceArea();
		return pLeft->tris.size() * pLeft->bounds.getSurfaceArea()/bounds.getSurfaceArea()  + pRight->tris.size() * pRight->bounds.getSurfaceArea()/bounds.getSurfaceArea();

		//return SAH_COST + pLeft->tris.size() * ( pLeft->bounds.getVolume() >= 0.1f ? pLeft->bounds.getVolume() : 0.1f)  / bounds.getVolume()
		//	+ pRight->tris.size() * ( pRight->bounds.getVolume() >= 0.1f ? pRight->bounds.getVolume() : 0.1f) / bounds.getVolume();
	} ; 

/*void BVH_Node::find_nodes2 ( Mesh & mesh, int axis, int div_i, int div_j ) {

	pLeft->tris.clear();
	pLeft->bounds.clear();

	pRight->tris.clear();
	pRight->bounds.clear();
	
	float dx = (bounds.max.m[axis] - bounds.min.m[axis]) / NUM_DIVIDE;
	float x1 = bounds.min.m[axis] + div_i * dx; 

	float x1_j = bounds.min.m[axis] + div_i * dx; 

	///*
	for ( int i=0; i<tris.size(); ++i ) {

		int & tri = tris[i];

		//bool b1 = false;

		if ( mesh.vertices [ mesh.indices[tri] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri] ].m[axis] <= x1
			&& mesh.vertices [ mesh.indices[tri+1] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri+1] ].m[axis] <= x1
			&& mesh.vertices [ mesh.indices[tri+2] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri+2] ].m[axis] <= x1 ) {
				//b1 = true;

				pLeft->tris.push_back ( tri );
				pLeft->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );

				continue;
		}

		if ( mesh.vertices [ mesh.indices[tri] ].m[axis] >= x1_j && mesh.vertices [ mesh.indices[tri] ].m[axis] <= bounds.max.m[axis]
			&& mesh.vertices [ mesh.indices[tri+1] ].m[axis] >= x1_j && mesh.vertices [ mesh.indices[tri+1] ].m[axis] <= bounds.max.m[axis]
			&& mesh.vertices [ mesh.indices[tri+2] ].m[axis] >= x1_j && mesh.vertices [ mesh.indices[tri+2] ].m[axis] <= bounds.max.m[axis] ) {

				pRight->tris.push_back ( tri );
				pRight->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );

				continue;
		}		
	

		if ( mesh.vertices [ mesh.indices[tri] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri] ].m[axis] <= x1
			|| mesh.vertices [ mesh.indices[tri+1] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri+1] ].m[axis] <= x1
			|| mesh.vertices [ mesh.indices[tri+2] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri+2] ].m[axis] <= x1 ) {

				pLeft->tris.push_back ( tri );
				pLeft->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );
		}

		if ( mesh.vertices [ mesh.indices[tri] ].m[axis] >= x1_j && mesh.vertices [ mesh.indices[tri] ].m[axis] <= bounds.max.m[axis]
			|| mesh.vertices [ mesh.indices[tri+1] ].m[axis] >= x1_j && mesh.vertices [ mesh.indices[tri+1] ].m[axis] <= bounds.max.m[axis]
			|| mesh.vertices [ mesh.indices[tri+2] ].m[axis] >= x1_j && mesh.vertices [ mesh.indices[tri+2] ].m[axis] <= bounds.max.m[axis] ) {

				pRight->tris.push_back ( tri );
				pRight->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );
		}		
	}

	pLeft->bounds.check_into ( bounds );
	pRight->bounds.check_into ( bounds );

	pLeft->bounds.max.m[axis] = fminf1 ( pLeft->bounds.max.m[axis], x1 );
	pRight->bounds.min.m[axis] = fmaxf1 ( pRight->bounds.min.m[axis], x1_j );

	//pLeft->bounds.max.m[axis] = x1;
	//pRight->bounds.min.m[axis] = x1;
} ;
*/

//void BVH_Node::find_nodes ( Mesh & mesh, int axis, int div_i, bool flag ) {
void BVH_Node::find_nodes ( Mesh & mesh, int axis, float testSplit, bool flag ) {

	pLeft->tris.clear();
	pLeft->bounds.clear();

	pRight->tris.clear();
	pRight->bounds.clear();
	
	//float dx = (bounds.max.m[axis] - bounds.min.m[axis]) / NUM_DIVIDE;
	//float x1 = bounds.min.m[axis] + div_i * dx; 

	///*
	for ( int i=0; i<tris.size(); ++i ) {

		int & tri = tris[i];

		float center = ( mesh.vertices [ mesh.indices[tri] ].m[axis] + mesh.vertices [ mesh.indices[tri+1] ].m[axis] + mesh.vertices [ mesh.indices[tri+2] ].m[axis] ) / 3;

		if (center < testSplit) {
			pLeft->tris.push_back ( tri );
			pLeft->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );
		} else {
			pRight->tris.push_back ( tri );
			pRight->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );
		}

		/*if ( mesh.vertices [ mesh.indices[tri] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri] ].m[axis] <= x1
			|| mesh.vertices [ mesh.indices[tri+1] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri+1] ].m[axis] <= x1
			|| mesh.vertices [ mesh.indices[tri+2] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri+2] ].m[axis] <= x1 ) {

				pLeft->tris.push_back ( tri );
				pLeft->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );

				//if ( flag ) {
				//	printf ("\n tri/3 = %d ", (tri/3) );
				//}
		} else if ( !( mesh.vertices [ mesh.indices[tri] ].m[axis] > x1 && mesh.vertices [ mesh.indices[tri] ].m[axis] <= bounds.max.m[axis]
					&& mesh.vertices [ mesh.indices[tri+1] ].m[axis] > x1 && mesh.vertices [ mesh.indices[tri+1] ].m[axis] <= bounds.max.m[axis]
					&& mesh.vertices [ mesh.indices[tri+2] ].m[axis] > x1 && mesh.vertices [ mesh.indices[tri+2] ].m[axis] <= bounds.max.m[axis] ) ) {

				pLeft->tris.push_back ( tri );
				pLeft->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );
		}

		if ( mesh.vertices [ mesh.indices[tri] ].m[axis] > x1 && mesh.vertices [ mesh.indices[tri] ].m[axis] <= bounds.max.m[axis]
			|| mesh.vertices [ mesh.indices[tri+1] ].m[axis] > x1 && mesh.vertices [ mesh.indices[tri+1] ].m[axis] <= bounds.max.m[axis]
			|| mesh.vertices [ mesh.indices[tri+2] ].m[axis] > x1 && mesh.vertices [ mesh.indices[tri+2] ].m[axis] <= bounds.max.m[axis] ) {

				pRight->tris.push_back ( tri );
				pRight->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );

		} else 	if ( !( mesh.vertices [ mesh.indices[tri] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri] ].m[axis] <= x1
					&& mesh.vertices [ mesh.indices[tri+1] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri+1] ].m[axis] <= x1
					&& mesh.vertices [ mesh.indices[tri+2] ].m[axis] >= bounds.min.m[axis] && mesh.vertices [ mesh.indices[tri+2] ].m[axis] <= x1 ) ) {

				pRight->tris.push_back ( tri );
				pRight->bounds.update3( mesh.vertices [ mesh.indices[tri] ], mesh.vertices [ mesh.indices[tri+1] ], mesh.vertices [ mesh.indices[tri+2] ] );
		}*/
	}

	//pLeft->bounds.check_into ( bounds );
	//pRight->bounds.check_into ( bounds );

	//pLeft->bounds.max.m[axis] = fminf1 ( pLeft->bounds.max.m[axis], x1 );
	//pRight->bounds.min.m[axis] = fmaxf1 ( pRight->bounds.min.m[axis], x1 );

	//pLeft->bounds.max.m[axis] = x1;
	//pRight->bounds.min.m[axis] = x1;
	//*/

	/*
	pLeft->tris.resize(tris.size());
	pRight->tris.resize(tris.size());
	for ( int i=0; i<tris.size(); ++i ) {
		pLeft->tris[i]  = tris[i];
		pRight->tris[i] = tris[i];
	}
	pLeft->bounds = bounds;
	pRight->bounds = bounds;
	*/
} ;


BVH_Node::~BVH_Node(void)
{
	if ( NULL != pLeft ) delete pLeft;
	if ( NULL != pRight ) delete pRight;
}
