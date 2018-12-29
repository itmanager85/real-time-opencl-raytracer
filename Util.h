#pragma once
#include "Defs.h"

#include "common.h"
#include "vectors_math.h"

namespace FW
{
//------------------------------------------------------------------------

class AABB
{
public:
	AABB        (void) : m_mn(FW_F32_MAX, FW_F32_MAX, FW_F32_MAX), m_mx(-FW_F32_MAX, -FW_F32_MAX, -FW_F32_MAX) {}
	AABB        (const float3& mn, const float3& mx) : m_mn(mn), m_mx(mx) {}

	void            grow        (const float3& pt)   { m_mn = fminf1(m_mn, pt); m_mx = fmaxf1(m_mx, pt); }
	void            grow        (const AABB& aabb)  { grow(aabb.m_mn); grow(aabb.m_mx); }
    void            intersect   (const AABB& aabb)  { m_mn = fmaxf1(m_mn, aabb.m_mn); m_mx = fminf1(m_mx, aabb.m_mx); }

	float           volume      (void) const        { if(!valid()) return 0.0f; return (m_mx.x-m_mn.x) * (m_mx.y-m_mn.y) * (m_mx.z-m_mn.z); }
    float           area        (void) const        { if(!valid()) return 0.0f; float3 d = m_mx - m_mn; return (d.x*d.y + d.y*d.z + d.z*d.x)*2.0f; }
    bool            valid       (void) const        { return m_mn.x<=m_mx.x && m_mn.y<=m_mx.y && m_mn.z<=m_mx.z; }
    float3           midPoint    (void) const        { return (m_mn+m_mx)*0.5f; }
    const float3&    minf         (void) const        { return m_mn; }
    const float3&    maxf         (void) const        { return m_mx; }
    float3&          minf         (void)              { return m_mn; }
    float3&          maxf         (void)              { return m_mx; }

private:
    float3           m_mn; // Vec3f
    float3           m_mx; // Vec3f
};

/*
class Util
{
public:
	Util(void);
	~Util(void);
};*/

}
