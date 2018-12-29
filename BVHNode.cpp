#include "BVHNode.h"

namespace FW
{

BVHNode::BVHNode(void)
{
}

void BVHNode::deleteSubtree()
{
    for(int i=0;i<getNumChildNodes();i++)
        getChildNode(i)->deleteSubtree();

    delete this;
}

BVHNode::~BVHNode(void)
{
	//deleteSubtree();
}

}