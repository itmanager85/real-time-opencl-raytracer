#pragma once
#include "Defs.h"
#include "Util.h"
#include "BVHNode.h"
#include "Platform.h"
#include "Mesh.h"

#include "BVH2.h"

namespace FW
{

class SplitBVHBuilder
{
private:
    enum
    {
        MaxDepth        = 64,
        MaxSpatialDepth = 48,
        NumSpatialBins  = 128,
    };

    struct Reference
    {
        S32                 triIdx;
        AABB                bounds;

        Reference(void) : triIdx(-1) {}

		//bool operator < (const Reference& rhs) const {
		//}
    };

    struct NodeSpec
    {
        S32                 numRef;
        AABB                bounds;

        NodeSpec(void) : numRef(0) {}
    };

    struct ObjectSplit
    {
        F32                 sah;
        S32                 sortDim;
        S32                 numLeft;
        AABB                leftBounds;
        AABB                rightBounds;

        ObjectSplit(void) : sah(FW_F32_MAX), sortDim(0), numLeft(0) {}
    };

    struct SpatialSplit
    {
        F32                 sah;
        S32                 dim;
        F32                 pos;

        SpatialSplit(void) : sah(FW_F32_MAX), dim(0), pos(0.0f) {}
    };

    struct SpatialBin
    {
        AABB                bounds;
        S32                 enter;
        S32                 exit;
    };

public:
	//SplitBVHBuilder(void);
	//~SplitBVHBuilder(void);

public:
    //                        SplitBVHBuilder     (BVH& bvh, const BVH::BuildParams& params);
	//						SplitBVHBuilder(Mesh& mesh, const Platform& platform);
							SplitBVHBuilder(BVH2& bvh);
                            ~SplitBVHBuilder    (void);

    BVHNode*                run                 (void);

private:
    static bool             sortCompare         (void* data, int idxA, int idxB);
    static void             sortSwap            (void* data, int idxA, int idxB);

    BVHNode*                buildNode           (NodeSpec spec, int level, F32 progressStart, F32 progressEnd);
    BVHNode*                createLeaf          (const NodeSpec& spec);

    ObjectSplit             findObjectSplit     (const NodeSpec& spec, F32 nodeSAH);
    void                    performObjectSplit  (NodeSpec& left, NodeSpec& right, const NodeSpec& spec, const ObjectSplit& split);

    SpatialSplit            findSpatialSplit    (const NodeSpec& spec, F32 nodeSAH);
    void                    performSpatialSplit (NodeSpec& left, NodeSpec& right, const NodeSpec& spec, const SpatialSplit& split);
    void                    splitReference      (Reference& left, Reference& right, const Reference& ref, int dim, F32 pos);

private:
                            SplitBVHBuilder     (const SplitBVHBuilder&); // forbidden
    SplitBVHBuilder&        operator=           (const SplitBVHBuilder&); // forbidden

private:
   BVH2&                    m_bvh;
	//Mesh&					m_Mesh;
    const Platform&         m_platform;
   // const BVH::BuildParams& m_params;

    vector<Reference>        m_refStack;
    F32                     m_minOverlap;
    vector<AABB>             m_rightBounds;
    S32                     m_sortDim;
    SpatialBin              m_bins[3][NumSpatialBins];

    //Timer                   m_progressTimer;
    S32                     m_numDuplicates;

	//vector<S32>				m_triIndices;

	F32						splitAlpha;     // spatial split area threshold
};

}