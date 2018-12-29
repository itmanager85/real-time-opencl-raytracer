#pragma once

#include "vectors_math.h"

#include <vector>

//#include <string>

//#include <iostream>

using namespace std;

#define DEBUG_ON 1 
//#define DEBUG_OFF 2 

//#define DEBUG_STATE DEBUG_OFF

// Obj loader ------------------------------------
/*struct TriangleFace 
{
	int v[3]; // vertex indices
};

struct TriangleMesh
{
	vector<float3> verts;
	vector<TriangleFace> faces;
	float3 bounding_box[2];
};*/


#include <vector>

using namespace std;


template <class T> void swap1(T& a, T& b) { T t = a; a = b; b = t; };
//template <class T> T min(T a, T b, T c) { return min(min(a, b), c); };

template <class T> T sqr(const T& a) { return a * a; }


struct TriangleFace 
{
	//union {
		int v[3]; // vertex indices

		int vn[3]; // for smothing triangle reflections

		int temp ; 
	//} ; 
};

struct TriangleMesh
{
	vector<float3> verts;

	vector<float3> norms; // for smothing triangle reflections

	vector<TriangleFace> faces;
	float3 bounding_box[2];

	//int getSize1() {
	//	return sizeof ( float3 ) * verts.size() + sizeof ( TriangleFace ) * faces.size() + sizeof ( bounding_box ) ; //sizeof ( float3 ) * 2;  
	//} ; 

	vector<float4> getVerts () {

		vector<float4> v ; 
		v.resize ( verts.size() )  ;

		for ( int i=0; i<verts.size(); ++i ) {
			v[i] = make_float4 ( verts[i].x, verts[i].y, verts[i].z, 1.0 );
		} ; 

		return  v ; 
	} ; 

};

/*
struct TriangleMeshCuda {

	int num_verts ; 
	float4 * verts;

	int num_faces ; 
	TriangleFace * faces;

	//float3 bounding_box[2];
} ; */

static int count = 0; 

static TriangleMesh mesh;

static TriangleMesh sphere;

static vector<float4> triangles; 


struct spec { 

	struct Ray {

		float3 o;
		float3 dir;
		float3 inv_dir;

		/*union {
			float3 o;
			//float3 ori;
			struct {
				float ox,oy,oz ;
			} ; 
		} ; 

		union {
			float3 dir;
			struct {
				float dx,dy,dz ; 
			} ; 
		} ; 

		union {
			float3 inv_dir;
			struct {
				float dxr,dyr,dzr ; 
			} ; 
		} ; */

		//__host__ __device__ 
		Ray(){};

		//__host__ __device__ 
		Ray( float3 &o, const float3 &d ) {	
			Ray::o = o;
			dir = normalize(d);
			inv_dir = make_float3(1.0/dir.x,1.0/dir.y,1.0/dir.z);
		}
		void print () {
			//printf ( "\n ray :: o = ( %.2f, %.2f, %.2f ), dir = ( %.2f, %.2f, %.2f ), inv_dir = ( %.2f, %.2f, %.2f ) ", ox,oy,oz, dx,dy,dz, dxr,dyr,dzr ) ;  
		} ;
	} ; 

	struct HitRecord
	{
		//__host__ __device__ 
		HitRecord() { t = UINT_MAX; hit_index = -1; color = make_float3(0,0,0); }

		//__host__ __device__ 
		void resetT(){ t = UINT_MAX; hit_index = -1; }

		float t;
		float3 color;
		float3 normal;
		int hit_index; 
		
	};


	// convert floating point rgb color to 8-bit integer
	static int rgbToInt(float r, float g, float b)
	{
		r = clamp(r, 0.0f, 255.0f);
		g = clamp(g, 0.0f, 255.0f);
		b = clamp(b, 0.0f, 255.0f);
		return (int(r)<<16) | (int(g)<<8) | int(b); // notice switch red and blue to counter the GL_BGRA
	}

	// intersection code --------------------------------------------------
	static int RayBoxIntersection(const float3 &BBMin, const float3 &BBMax, const float3 &RayOrg, const float3 &RayDirInv, float &tmin, float &tmax)
	{
		float l1   = (BBMin.x - RayOrg.x) * RayDirInv.x;
		float l2   = (BBMax.x - RayOrg.x) * RayDirInv.x;
		tmin = fminf1(l1,l2);
		tmax = fmaxf1(l1,l2);

		l1   = (BBMin.y - RayOrg.y) * RayDirInv.y;
		l2   = (BBMax.y - RayOrg.y) * RayDirInv.y;
		tmin = fmaxf1(fminf1(l1,l2), tmin);
		tmax = fminf1(fmaxf1(l1,l2), tmax);

		l1   = (BBMin.z - RayOrg.z) * RayDirInv.z;
		l2   = (BBMax.z - RayOrg.z) * RayDirInv.z;
		tmin = fmaxf1(fminf1(l1,l2), tmin);
		tmax = fminf1(fmaxf1(l1,l2), tmax);

		return ((tmax >= tmin) && (tmax >= 0.0f));
	}

	// the classic ray triangle intersection: http://www.cs.virginia.edu/~gfx/Courses/2003/ImageSynthesis/papers/Acceleration/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
	static float RayTriangleIntersection(const Ray &r, 
											 const float3 &v0, 
											 const float3 &edge1, 
											 const float3 &edge2)
	{  

		float3 tvec = r.o- v0;  
		float3 pvec = cross(r.dir, edge2);  
		float  det  = dot(edge1, pvec);  

		//det = __fdividef(1.0f, det);  
		det = 1.0f / det;  

		float u = dot(tvec, pvec) * det;  

		if (u < 0.0f || u > 1.0f)  
			return -1.0f;  

		float3 qvec = cross(tvec, edge1);  

		float v = dot(r.dir, qvec) * det;  

		if (v < 0.0f || (u + v) > 1.0f)  
			return -1.0f;  

		return dot(edge2, qvec) * det;  
	}  

	static int RaySphereIntersection(const Ray  &ray, const float3 sphere_center, const float sphere_radius, float &t)
	{
		float b, c, d;

		float3 sr = ray.o - sphere_center;
		b =  dot(sr,ray.dir);
		c = dot(sr,sr) - (sphere_radius*sphere_radius);
		d = b*b - c;
		if (d > 0) 
		{
			float e = sqrt(d);
			float t0 = -b-e;
			if(t0 < 0)
				t = -b+e;
			else
				t = fminf1(-b-e,-b+e);
			return 1;
		}
		return 0;
	}

} ;

