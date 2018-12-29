#pragma once

#include "Defs.h"

namespace FW
{

class LeafNode;
class BVHNode;
class Platform
{
public:
    //Platform()                                                                                                          { m_name=String("Default"); m_SAHNodeCost = 1.f; m_SAHTriangleCost = 1.f; m_nodeBatchSize = 1; m_triBatchSize = 1; m_minLeafSize=1; m_maxLeafSize=0x7FFFFFF; }
    //Platform(const String& name,float nodeCost=1.f, float triCost=1.f, S32 nodeBatchSize=1, S32 triBatchSize=1) { m_name=name; m_SAHNodeCost = nodeCost; m_SAHTriangleCost = triCost; m_nodeBatchSize = nodeBatchSize; m_triBatchSize = triBatchSize; m_minLeafSize=1; m_maxLeafSize=0x7FFFFFF; }

	//Platform()                                                                               { m_SAHNodeCost = 1.f; m_SAHTriangleCost = 1.f; m_nodeBatchSize = 1; m_triBatchSize = 1; m_minLeafSize=1; m_maxLeafSize=0x7FFFFFF; }
    Platform(float nodeCost=1.f, float triCost=1.f, S32 nodeBatchSize=1, S32 triBatchSize=1) { m_SAHNodeCost = nodeCost; m_SAHTriangleCost = triCost; m_nodeBatchSize = nodeBatchSize; m_triBatchSize = triBatchSize; m_minLeafSize=1; m_maxLeafSize=0x7FFFFFF; }

    //const String&   getName() const                     { return m_name; }

    // SAH weights
    float getSAHTriangleCost() const                    { return m_SAHTriangleCost; }
    float getSAHNodeCost() const                        { return m_SAHNodeCost; }

    // SAH costs, raw and batched
    float getCost(int numChildNodes,int numTris) const  { return getNodeCost(numChildNodes) + getTriangleCost(numTris); }
    float getTriangleCost(S32 n) const                  { return roundToTriangleBatchSize(n) * m_SAHTriangleCost; }
    float getNodeCost(S32 n) const                      { return roundToNodeBatchSize(n) * m_SAHNodeCost; }

    // batch processing (how many ops at the price of one)
    S32   getTriangleBatchSize() const                  { return m_triBatchSize; }
    S32   getNodeBatchSize() const                      { return m_nodeBatchSize; }
    void  setTriangleBatchSize(S32 triBatchSize)        { m_triBatchSize = triBatchSize; }
    void  setNodeBatchSize(S32 nodeBatchSize)           { m_nodeBatchSize= nodeBatchSize; }
    S32   roundToTriangleBatchSize(S32 n) const         { return ((n+m_triBatchSize-1)/m_triBatchSize)*m_triBatchSize; }
    S32   roundToNodeBatchSize(S32 n) const             { return ((n+m_nodeBatchSize-1)/m_nodeBatchSize)*m_nodeBatchSize; }

    // leaf preferences
    void  setLeafPreferences(S32 minSize,S32 maxSize)   { m_minLeafSize=minSize; m_maxLeafSize=maxSize; }
    S32   getMinLeafSize() const                        { return m_minLeafSize; }
    S32   getMaxLeafSize() const                        { return m_maxLeafSize; }

    //U32   computeHash() const                           { return hashBits(hash<String>(m_name), floatToBits(m_SAHNodeCost), floatToBits(m_SAHTriangleCost), hashBits(m_triBatchSize, m_nodeBatchSize, m_minLeafSize, m_maxLeafSize)); }

private:
    //String  m_name;
    float   m_SAHNodeCost;
    float   m_SAHTriangleCost;
    S32     m_triBatchSize;
    S32     m_nodeBatchSize;
    S32     m_minLeafSize;
    S32     m_maxLeafSize;
};


} //
