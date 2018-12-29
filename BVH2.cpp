#include "BVH2.h"

#include "SplitBVHBuilder.h"

using namespace FW;

//BVH2::BVH2(void)
//{
//}

BVH2::BVH2 (Mesh* mesh ) {
	//m_scene = mesh;
	m_platform.setLeafPreferences ( 1, 8 );

	//m_root = NULL;
	//m_root = SplitBVHBuilder(*this).run();

	//m_scene = NULL;
	m_root = NULL;
	setMesh ( m_scene );
}

void BVH2::setMesh(Mesh* mesh) {
	
	clear ();

	m_scene = mesh;

	if ( NULL != mesh ) 
		m_root = SplitBVHBuilder(*this).run();
};

void BVH2::clear () {
	if ( NULL != m_root ) 
		m_root->deleteSubtree();

	m_root = NULL;
} ; 


BVH2::~BVH2(void)
{
	//if ( NULL != m_root ) m_root->deleteSubtree();
	clear ();
}
