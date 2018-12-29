#include "BVH.h"

#include <stddef.h>

BVH::BVH(void)
{
	root = NULL;
}


BVH::~BVH(void)
{
	if ( NULL != root ) delete root;
}
