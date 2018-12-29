#pragma once

//#ifndef _VEC_MATH
//#define _VEC_MATH

#include "vectors_math.h"

#include <math.h>

/*float fmax ( float a, float b ) {
	return a >= b ? a : b;
} ; 

float fmin ( float a, float b ) {
	return a <= b ? a : b;
} ; */

float rsqrtf(float x) {
    return 1.0f / sqrtf(x);
}

	float2 make_float2 ( float x, float y ) {
		return float2 ( x,y );
	} ; 

	float3 operator / (float3 &a, float3 &b) {
		return make_float3( a.x / b.x, a.y / b.y, a.z / b.z );
	} ; 

	
	float4 fminf1 ( const float4 &v0, const float4 &v1 ) {
		return make_float4 ( fminf1 ( v0.x, v1.x ), fminf1 ( v0.y, v1.y ), fminf1 ( v0.z, v1.z ), 1.0f );
	} ; 

	float4 fmaxf1 ( const float4 &v0, const float4 &v1 ) {
		return make_float4 ( fmaxf1 ( v0.x, v1.x ), fmaxf1 ( v0.y, v1.y ), fmaxf1 ( v0.z, v1.z ), 1.0f );
	} ; 

	float3 fminf1 ( const float3 &v0, const float3 &v1 ) {
		return make_float3 ( fminf1 ( v0.x, v1.x ), fminf1 ( v0.y, v1.y ), fminf1 ( v0.z, v1.z ) );
	} ; 

	float3 fmaxf1 ( const float3 &v0, const float3 &v1 ) {
		return make_float3 ( fmaxf1 ( v0.x, v1.x ), fmaxf1 ( v0.y, v1.y ), fmaxf1 ( v0.z, v1.z ) );
	} ; 

	float4 operator-(float4 &a, float4 &b)
	{
		return make_float4( a.x - b.x, a.y - b.y, a.z - b.z, 1.0f );
	}


	float3 make_float3 ( float x, float y, float z ) {
		return float3 ( x,y,z );
	} ; 

	float3 operator-(float3 &a) {
		return make_float3(-a.x, -a.y, -a.z);
	} ; 

	float3 operator+(float3 a, float3 b) {
		return make_float3(a.x + b.x, a.y + b.y, a.z + b.z );
	}

	void operator*=(float3 &a, float b) {
		a.x *= b; a.y *= b; a.z *= b;
	}

	float3 operator*(float b, float3 a) {
		return make_float3(b * a.x, b * a.y, b * a.z );
	}

	float dot(float3 a, float3 b) {
		return a.x * b.x + a.y * b.y + a.z * b.z;
	} ; 

	float3 cross(float3 a, float3 b) { 
		return make_float3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x); 
	}

	float3 normalize(float3 v) {
		float invLen = rsqrtf(dot(v, v));
		return invLen * v;
	}

	float4 make_float4 ( float x, float y, float z, float w ) {
		return float4 ( x,y,z,w );
	} ; 

	float4 operator-(float4 &a) {
		return make_float4(-a.x, -a.y, -a.z,-a.w);
	} ; 

	float4 operator+(float4 a, float4 b) {
		return make_float4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w );
	}

	void operator*=(float4 &a, float b) {
		a.x *= b; a.y *= b; a.z *= b; a.w *= b;
	}

	float4 operator*(float b, float4 a) {
		return make_float4(b * a.x, b * a.y, b * a.z, b * a.w );
	}

	float dot(float4 a, float4 b) {
		return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
	} ; 


	float4 normalize(float4 v) {
		const float invLen = rsqrtf(dot(v, v));
		return invLen * v;
	}


	float imin1(const int a, const int b) {
		return a < b ? a : b;
	}

	float imax1(const int a, const int b) {
		return a > b ? a : b;
	}


	float fminf1(const float a, const float b)
	{
	  return a < b ? a : b;
	}

	float fmaxf1(const float a, const float b)
	{
	  return a > b ? a : b;
	}

	float clamp(const float f, const float a, const float b) {
		return fmaxf1(a, fminf1(f, b));
	}

float3 operator-(float3 a, float3 b)
{
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}

float3 operator+(float3 a, float b)
{
    return make_float3(a.x + b, a.y + b, a.z + b);
}

void operator+=(float3 &a, float b)
{
    a.x += b; a.y += b; a.z += b;
}

	float3 operator*(float3 a, float b)
{
    return make_float3(a.x * b, a.y * b, a.z * b);
}


	float length(float3 v)
{
    return sqrtf(dot(v, v));
}

	void operator+=(float3 &a, float3 b)
{
    a.x += b.x; a.y += b.y; a.z += b.z;
}

	float3 operator/(float3 a, float b)
{
    return make_float3(a.x / b, a.y / b, a.z / b);
}

	void operator/=(float3 &a, float3 b)
{
    a.x /= b.x; a.y /= b.y; a.z /= b.z;
}

	void operator/=(float3 &a, float b)
{
    a.x /= b; a.y /= b; a.z /= b;
}

	float3 get_float3( float4 &v ) {
		return make_float3( v.x, v.y, v.z );
	} ;

	float3 operator/(float b, float3 a) {
		return make_float3( b/a.x, b/a.y, b/a.z );
	};

	float4 operator*(float4 a, float b) {
		return make_float4(b * a.x, b * a.y, b * a.z, b * a.w );
	} ; 

	float4 lerp(const float4& a, const float4& b, const float t) { return (float4)(a * ((float)1 - t) + b * t); };

	float fminf1(float a, float b, float c) { return fminf1(fminf1(a, b), c); };

	float fmaxf1(float a, float b, float c) { return fmaxf1(fmaxf1(a, b), c); };

	int3 make_int3 ( int x, int y, int z ) {
		return int3( x, y, z );
	} ; 

	int3 operator-(int3 &a, int3 &b) {
		return int3(a.x - b.x, a.y - b.y, a.z - b.z );
	} ; 

	int3 operator*(int3 a, float3 b) {
		return int3(a.x * b.x, a.y * b.y, a.z * b.z );
	} ; 

	float3 operator*(float3 a, float3 b) {
		return float3(a.x * b.x, a.y * b.y, a.z * b.z );
	} ; 

	int3 min3i ( const int3 &v0, const int a ) {
		return int3( min(v0.x,a), min(v0.y,a), min(v0.z,a) );
	} ; 

	int3 max3i ( const int3 &v0, const int a ) {
		return int3( max(v0.x,a), max(v0.y,a), max(v0.z,a) );
	} ; 

	int3 clamp(const int3 i3, const int a, const int b) {
		return max3i(min3i(i3, b), a);
	} ; 

	int3 max3i ( const int3 &v0, const int3 &a ) {
		return int3( max(v0.x,a.x), max(v0.y,a.y), max(v0.z,a.z) );
	} ; 

	int3 clamp(const int3 i3, const int3 &a, const int b) {
		return max3i(min3i(i3, b), a);
	} ; 

	int3 make_int3 ( const float3 & a ) {
		return int3( a.x, a.y, a.z );
	} ; 

	float lerp(const float a, const float b, const float t) { return (float)(a * ((float)1 - t) + b * t); };

	float fminf1 ( const float3& a ) {
		return fminf1 ( a.x, a.y, a.z );
	}; 

	float fmaxf1 ( const float3& a ) {
		return fmaxf1 ( a.x, a.y, a.z );
	}; 

	float fsumf ( const float3& a ) {
		return a.x + a.y + a.z;
	}; 

	float3 make_float3 ( const float4& v ) {
		return float3 ( v.x, v.y, v.z );
	} ; 

	float4 make_float4 ( const float3& v ) {
		return float4 ( v.x, v.y, v.z, 1.0f );
	} ; 

	float4& operator+=(float4 &a, float4 &b) {
		a = a + b;
		return a;
	}


//#endif  /* _VEC_MATH */