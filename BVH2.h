#pragma once

#include "BVHNode.h"
#include "Platform.h"
#include "Mesh.h"

namespace FW
{

class BVH2
{
public:
	//BVH2(void);
	BVH2                 (Mesh* mesh = NULL);
	~BVH2(void);

	void				setMesh(Mesh* mesh);
	void				clear ();

	Mesh*              getScene                (void) const            { return m_scene; }
    const Platform&     
		getPlatform1             
		(void) 
		const			
	{ return m_platform; }
    BVHNode*            getRoot                 (void) const            { return m_root; }

	vector<S32>&         getTriIndices           (void)                  { return m_triIndices; }
    const vector<S32>&   getTriIndices           (void) const            { return m_triIndices; }

//private:
public:
	Mesh*              m_scene;
    Platform            m_platform;

    BVHNode*            m_root;
    vector<S32>          m_triIndices;
};

}