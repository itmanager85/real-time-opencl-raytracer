#pragma once

#include "Defs.h"
#include "Util.h"

namespace FW
{

class BVHNode
{
public:
	BVHNode(void);

	virtual bool        isLeaf() const = 0;
    virtual S32         getNumChildNodes() const = 0;
    virtual BVHNode*    getChildNode(S32 i) const   = 0;
    virtual S32         getNumTriangles() const { return 0; }

    float       getArea() const     { return m_bounds.area(); }

    AABB        m_bounds;

	void    deleteSubtree();

	~BVHNode(void);
};


class InnerNode : public BVHNode
{
public:
    InnerNode(const AABB& bounds,BVHNode* child0,BVHNode* child1)   { m_bounds=bounds; m_children[0]=child0; m_children[1]=child1; }

    bool        isLeaf() const                  { return false; }
    S32         getNumChildNodes() const        { return 2; }
    BVHNode*    getChildNode(S32 i) const       { if(i>=0 && i<2) return m_children[i]; return NULL; } // FW_ASSERT

    BVHNode*    m_children[2];
};


class LeafNode : public BVHNode
{
public:
    LeafNode(const AABB& bounds,int lo,int hi)  { m_bounds=bounds; m_lo=lo; m_hi=hi; }
    LeafNode(const LeafNode& s)                 { *this = s; }

    bool        isLeaf() const                  { return true; }
    S32         getNumChildNodes() const        { return 0; }
    BVHNode*    getChildNode(S32) const         { return NULL; }

    S32         getNumTriangles() const         { return m_hi-m_lo; }
    S32         m_lo;
    S32         m_hi;
};

}