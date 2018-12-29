#include "SplitBVHBuilder.h"

#include "common.h"

#include "sort.h"

//#include <algorithm>    // std::sort

using namespace FW;

//------------------------------------------------------------------------
SplitBVHBuilder::SplitBVHBuilder(BVH2& bvh)
:   m_bvh           (bvh),
    m_platform      (bvh.getPlatform1()),
    //m_params        (params),
    m_minOverlap    (0.0f),
    m_sortDim       (-1)
{
	splitAlpha      = 1.0e-5f;
}
/*
SplitBVHBuilder::SplitBVHBuilder(Mesh& mesh, const Platform& platform)
:   m_Mesh          (mesh),
    m_platform      (platform),
    //m_params        (params),
    m_minOverlap    (0.0f),
    m_sortDim       (-1)
{
	splitAlpha      = 1.0e-5f;
}
*/
//------------------------------------------------------------------------


SplitBVHBuilder::~SplitBVHBuilder(void)
{
}

//------------------------------------------------------------------------

BVHNode* SplitBVHBuilder::run(void)
{
    // Initialize reference stack and determine root bounds.

	printf ( "\n SplitBVHBuilder::run(void) " );

    const int3* tris = (const int3*) m_bvh.getScene()->getTrisPtr(); //&m_bvh.getScene()->indices[0]; //&m_Mesh.indices[0]; //m_bvh.getScene()->getTriVtxIndexBuffer().getPtr();
    const float4* verts = (const float4*) m_bvh.getScene()->getVertsPtr(); // &m_Mesh.vertices[0]; //m_bvh.getScene()->getVtxPosBuffer().getPtr();

    NodeSpec rootSpec;
    rootSpec.numRef = m_bvh.getScene()->getNumTriangles(); //m_Mesh.indices.size()/3; //->getNumTriangles();
    m_refStack.resize(rootSpec.numRef);

    for (int i = 0; i < rootSpec.numRef; i++)
    {
        m_refStack[i].triIdx = i;
        for (int j = 0; j < 3; j++)
            m_refStack[i].bounds.grow( make_float3( verts[tris[i].m[j]] ) );
        rootSpec.bounds.grow(m_refStack[i].bounds);
    }

    // Initialize rest of the members.

    m_minOverlap = rootSpec.bounds.area() * splitAlpha; // m_params.splitAlpha;
	m_rightBounds.clear();
    m_rightBounds.resize(imax1(rootSpec.numRef, (int)NumSpatialBins) - 1); // reset
    m_numDuplicates = 0;
    //m_progressTimer.start();

    // Build recursively.

    BVHNode* root = buildNode(rootSpec, 0, 0.0f, 1.0f);
    //m_bvh.getTriIndices().compact();

    // Done.

    //if (m_params.enablePrints)
    //    printf("SplitBVHBuilder: progress %.0f%%, duplicates %.0f%%\n",
    //        100.0f, (F32)m_numDuplicates / (F32)m_bvh.getScene()->getNumTriangles() * 100.0f);
    return root;
}

//------------------------------------------------------------------------

bool SplitBVHBuilder::sortCompare(void* data, int idxA, int idxB)
{
    const SplitBVHBuilder* ptr = (const SplitBVHBuilder*)data;
    int dim = ptr->m_sortDim;
    const Reference& ra = ptr->m_refStack[idxA];
    const Reference& rb = ptr->m_refStack[idxB];
    F32 ca = ra.bounds.minf().m[dim] + ra.bounds.maxf().m[dim];
    F32 cb = rb.bounds.minf().m[dim] + rb.bounds.maxf().m[dim];
    return (ca < cb || (ca == cb && ra.triIdx < rb.triIdx));
}

//------------------------------------------------------------------------

void SplitBVHBuilder::sortSwap(void* data, int idxA, int idxB)
{
    SplitBVHBuilder* ptr = (SplitBVHBuilder*)data;
    swap(ptr->m_refStack[idxA], ptr->m_refStack[idxB]);
}


//------------------------------------------------------------------------

BVHNode* SplitBVHBuilder::buildNode(NodeSpec spec, int level, F32 progressStart, F32 progressEnd)
{
    // Display progress.

	/*
    if (m_params.enablePrints && m_progressTimer.getElapsed() >= 1.0f)
    {
        printf("SplitBVHBuilder: progress %.0f%%, duplicates %.0f%%\r",
            progressStart * 100.0f, (F32)m_numDuplicates / (F32)m_bvh.getScene()->getNumTriangles() * 100.0f);
        m_progressTimer.start();
    }*/

    // Remove degenerates.
    {
        int firstRef = m_refStack.size() - spec.numRef;
        for (int i = m_refStack.size() - 1; i >= firstRef; i--)
        {
            float3 size = m_refStack[i].bounds.maxf() - m_refStack[i].bounds.minf();
            if (fminf1(size) < 0.0f || fsumf(size) == fmaxf1(size)) {
                //m_refStack.removeSwap(i);
				m_refStack[i] = m_refStack[m_refStack.size() - 1];
				m_refStack.pop_back();
			}
        }
        spec.numRef = m_refStack.size() - firstRef;
    }

    // Small enough or too deep => create leaf.

    if (spec.numRef <= m_platform.getMinLeafSize() || level >= MaxDepth)
        return createLeaf(spec);

    // Find split candidates.

    F32 area = spec.bounds.area();
    F32 leafSAH = area * m_platform.getTriangleCost(spec.numRef);
    F32 nodeSAH = area * m_platform.getNodeCost(2);
    ObjectSplit object = findObjectSplit(spec, nodeSAH);

    SpatialSplit spatial;
    if (level < MaxSpatialDepth)
    {
        AABB overlap = object.leftBounds;
        overlap.intersect(object.rightBounds);
        if (overlap.area() >= m_minOverlap)
            spatial = findSpatialSplit(spec, nodeSAH);
    }

    // Leaf SAH is the lowest => create leaf.

    F32 minSAH = fminf1(leafSAH, object.sah, spatial.sah);
    if (minSAH == leafSAH && spec.numRef <= m_platform.getMaxLeafSize())
        return createLeaf(spec);

    // Perform split.

    NodeSpec left, right;
    if (minSAH == spatial.sah)
        performSpatialSplit(left, right, spec, spatial);
    if (!left.numRef || !right.numRef)
        performObjectSplit(left, right, spec, object);

    // Create inner node.

    m_numDuplicates += left.numRef + right.numRef - spec.numRef;
    F32 progressMid = lerp(progressStart, progressEnd, (F32)right.numRef / (F32)(left.numRef + right.numRef));
    BVHNode* rightNode = buildNode(right, level + 1, progressStart, progressMid);
    BVHNode* leftNode = buildNode(left, level + 1, progressMid, progressEnd);
    return new InnerNode(spec.bounds, leftNode, rightNode);
}

//------------------------------------------------------------------------


BVHNode* SplitBVHBuilder::createLeaf(const NodeSpec& spec)
{
	vector<S32>& tris = m_bvh.getTriIndices(); //m_triIndices; //m_bvh.getTriIndices();
    for (int i = 0; i < spec.numRef; i++) {
		tris.push_back( m_refStack[m_refStack.size()-1].triIdx ); // add
		m_refStack.pop_back(); //removeLast()
	}
    return new LeafNode(spec.bounds, tris.size() - spec.numRef, tris.size());
}

//------------------------------------------------------------------------

SplitBVHBuilder::ObjectSplit SplitBVHBuilder::findObjectSplit(const NodeSpec& spec, F32 nodeSAH)
{
    ObjectSplit split;
    const Reference* refPtr = &m_refStack[m_refStack.size() - spec.numRef];
    F32 bestTieBreak = FW_F32_MAX;

    // Sort along each dimension.

    for (m_sortDim = 0; m_sortDim < 3; m_sortDim++)
    {
        sort(this, m_refStack.size() - spec.numRef, m_refStack.size(), sortCompare, sortSwap);

        // Sweep right to left and determine bounds.

        AABB rightBounds;
        for (int i = spec.numRef - 1; i > 0; i--)
        {
            rightBounds.grow(refPtr[i].bounds);
            m_rightBounds[i - 1] = rightBounds;
        }

        // Sweep left to right and select lowest SAH.

        AABB leftBounds;
        for (int i = 1; i < spec.numRef; i++)
        {
            leftBounds.grow(refPtr[i - 1].bounds);
            F32 sah = nodeSAH + leftBounds.area() * m_platform.getTriangleCost(i) + m_rightBounds[i - 1].area() * m_platform.getTriangleCost(spec.numRef - i);
            F32 tieBreak = sqr((F32)i) + sqr((F32)(spec.numRef - i));
            if (sah < split.sah || (sah == split.sah && tieBreak < bestTieBreak))
            {
                split.sah = sah;
                split.sortDim = m_sortDim;
                split.numLeft = i;
                split.leftBounds = leftBounds;
                split.rightBounds = m_rightBounds[i - 1];
                bestTieBreak = tieBreak;
            }
        }
    }
    return split;
}

//------------------------------------------------------------------------

void SplitBVHBuilder::performObjectSplit(NodeSpec& left, NodeSpec& right, const NodeSpec& spec, const ObjectSplit& split)
{
    m_sortDim = split.sortDim;
    sort(this, m_refStack.size() - spec.numRef, m_refStack.size(), sortCompare, sortSwap);
	 //std::sort ( m_refStack.end()-spec.numRef, m_refStack.end(), sortCompare );

    left.numRef = split.numLeft;
    left.bounds = split.leftBounds;
    right.numRef = spec.numRef - split.numLeft;
    right.bounds = split.rightBounds;
}

//------------------------------------------------------------------------

SplitBVHBuilder::SpatialSplit SplitBVHBuilder::findSpatialSplit(const NodeSpec& spec, F32 nodeSAH)
{
    // Initialize bins.

    float3 origin = spec.bounds.minf();
    float3 binSize = (spec.bounds.maxf() - origin) * (1.0f / (F32)NumSpatialBins);
    float3 invBinSize = 1.0f / binSize;

    for (int dim = 0; dim < 3; dim++)
    {
        for (int i = 0; i < NumSpatialBins; i++)
        {
            SpatialBin& bin = m_bins[dim][i];
            bin.bounds = AABB();
            bin.enter = 0;
            bin.exit = 0;
        }
    }

    // Chop references into bins.

    for (int refIdx = m_refStack.size() - spec.numRef; refIdx < m_refStack.size(); refIdx++)
    {
        const Reference& ref = m_refStack[refIdx];
        int3 firstBin = clamp(make_int3((ref.bounds.minf() - origin) * invBinSize), 0, NumSpatialBins - 1);
        int3 lastBin = clamp(make_int3((ref.bounds.maxf() - origin) * invBinSize), firstBin, NumSpatialBins - 1);

        for (int dim = 0; dim < 3; dim++)
        {
            Reference currRef = ref;
            for (int i = firstBin.m[dim]; i < lastBin.m[dim]; i++)
            {
                Reference leftRef, rightRef;
                splitReference(leftRef, rightRef, currRef, dim, origin.m[dim] + binSize.m[dim] * (F32)(i + 1));
                m_bins[dim][i].bounds.grow(leftRef.bounds);
                currRef = rightRef;
            }
            m_bins[dim][lastBin.m[dim]].bounds.grow(currRef.bounds);
            m_bins[dim][firstBin.m[dim]].enter++;
            m_bins[dim][lastBin.m[dim]].exit++;
        }
    }

    // Select best split plane.

    SpatialSplit split;
    for (int dim = 0; dim < 3; dim++)
    {
        // Sweep right to left and determine bounds.

        AABB rightBounds;
        for (int i = NumSpatialBins - 1; i > 0; i--)
        {
            rightBounds.grow(m_bins[dim][i].bounds);
            m_rightBounds[i - 1] = rightBounds;
        }

        // Sweep left to right and select lowest SAH.

        AABB leftBounds;
        int leftNum = 0;
        int rightNum = spec.numRef;

        for (int i = 1; i < NumSpatialBins; i++)
        {
            leftBounds.grow(m_bins[dim][i - 1].bounds);
            leftNum += m_bins[dim][i - 1].enter;
            rightNum -= m_bins[dim][i - 1].exit;

            F32 sah = nodeSAH + leftBounds.area() * m_platform.getTriangleCost(leftNum) + m_rightBounds[i - 1].area() * m_platform.getTriangleCost(rightNum);
            if (sah < split.sah)
            {
                split.sah = sah;
                split.dim = dim;
                split.pos = origin.m[dim] + binSize.m[dim] * (F32)i;
            }
        }
    }
    return split;
}

//------------------------------------------------------------------------

void SplitBVHBuilder::performSpatialSplit(NodeSpec& left, NodeSpec& right, const NodeSpec& spec, const SpatialSplit& split)
{
    // Categorize references and compute bounds.
    //
    // Left-hand side:      [leftStart, leftEnd[
    // Uncategorized/split: [leftEnd, rightStart[
    // Right-hand side:     [rightStart, refs.getSize()[

    vector<Reference>& refs = m_refStack;
    int leftStart = refs.size() - spec.numRef;
    int leftEnd = leftStart;
    int rightStart = refs.size();
    left.bounds = right.bounds = AABB();

    for (int i = leftEnd; i < rightStart; i++)
    {
        // Entirely on the left-hand side?

        if (refs[i].bounds.maxf().m[split.dim] <= split.pos)
        {
            left.bounds.grow(refs[i].bounds);
            swap1(refs[i], refs[leftEnd++]);
        }

        // Entirely on the right-hand side?

        else if (refs[i].bounds.minf().m[split.dim] >= split.pos)
        {
            right.bounds.grow(refs[i].bounds);
            swap1(refs[i--], refs[--rightStart]);
        }
    }

    // Duplicate or unsplit references intersecting both sides.

    while (leftEnd < rightStart)
    {
        // Split reference.

        Reference lref, rref;
        splitReference(lref, rref, refs[leftEnd], split.dim, split.pos);

        // Compute SAH for duplicate/unsplit candidates.

        AABB lub = left.bounds;  // Unsplit to left:     new left-hand bounds.
        AABB rub = right.bounds; // Unsplit to right:    new right-hand bounds.
        AABB ldb = left.bounds;  // Duplicate:           new left-hand bounds.
        AABB rdb = right.bounds; // Duplicate:           new right-hand bounds.
        lub.grow(refs[leftEnd].bounds);
        rub.grow(refs[leftEnd].bounds);
        ldb.grow(lref.bounds);
        rdb.grow(rref.bounds);

        F32 lac = m_platform.getTriangleCost(leftEnd - leftStart);
        F32 rac = m_platform.getTriangleCost(refs.size() - rightStart);
        F32 lbc = m_platform.getTriangleCost(leftEnd - leftStart + 1);
        F32 rbc = m_platform.getTriangleCost(refs.size() - rightStart + 1);

        F32 unsplitLeftSAH = lub.area() * lbc + right.bounds.area() * rac;
        F32 unsplitRightSAH = left.bounds.area() * lac + rub.area() * rbc;
        F32 duplicateSAH = ldb.area() * lbc + rdb.area() * rbc;
        F32 minSAH = fminf1(unsplitLeftSAH, unsplitRightSAH, duplicateSAH);

        // Unsplit to left?

        if (minSAH == unsplitLeftSAH)
        {
            left.bounds = lub;
            leftEnd++;
        }

        // Unsplit to right?

        else if (minSAH == unsplitRightSAH)
        {
            right.bounds = rub;
            swap(refs[leftEnd], refs[--rightStart]);
        }

        // Duplicate?

        else
        {
            left.bounds = ldb;
            right.bounds = rdb;
            refs[leftEnd++] = lref;
            refs.push_back(rref);
        }
    }

    left.numRef = leftEnd - leftStart;
    right.numRef = refs.size() - rightStart;
}

//------------------------------------------------------------------------

void SplitBVHBuilder::splitReference(Reference& left, Reference& right, const Reference& ref, int dim, F32 pos)
{
    // Initialize references.

    left.triIdx = right.triIdx = ref.triIdx;
    left.bounds = right.bounds = AABB();


    // Loop over vertices/edges.

    const int3* tris = (const int3*) m_bvh.getScene()->getTrisPtr(); //&m_Mesh.indices[0]; //m_bvh.getScene()->getTriVtxIndexBuffer().getPtr();
    const float4* verts = (const float4*) m_bvh.getScene()->getVertsPtr(); //&m_Mesh.vertices[0]; //m_bvh.getScene()->getVtxPosBuffer().getPtr();
    const int3& inds = tris[ref.triIdx];
    const float4* v1 = &verts[inds.z];

    for (int i = 0; i < 3; i++)
    {
        const float4* v0 = v1;
        v1 = &verts[inds.m[i]];
        F32 v0p = v0->get(dim);
        F32 v1p = v1->get(dim);

        // Insert vertex to the boxes it belongs to.

        if (v0p <= pos)
            left.bounds.grow( make_float3(v0->x, v0->y, v0->z) ); // *v0
        if (v0p >= pos)
            right.bounds.grow( make_float3(v0->x, v0->y, v0->z) ); // *v0

        // Edge intersects the plane => insert intersection to both boxes.

        if ((v0p < pos && v1p > pos) || (v0p > pos && v1p < pos))
        {
            float4 t = lerp(*v0, *v1, clamp((pos - v0p) / (v1p - v0p), 0.0f, 1.0f));
            left.bounds.grow( make_float3(t.x, t.y, t.z) ); // t
            right.bounds.grow( make_float3(t.x, t.y, t.z) );
        }
    }

    // Intersect with original bounds.

    left.bounds.maxf().m[dim] = pos;
    right.bounds.minf().m[dim] = pos;
    left.bounds.intersect(ref.bounds);
    right.bounds.intersect(ref.bounds);
}

//------------------------------------------------------------------------