#pragma once

#include <math.h>

#include "vectors_math.h"
#include "stdio.h"

class AABB
{
public:

float4 min, max; 

	AABB () {
		clear() ; 
	} ;

	void clear() {
		min = make_float4 ( exp( 80.0f ), exp( 80.0f ), exp( 80.0f ), 1.0f );
		max = make_float4 ( -exp( 80.0f ), -exp( 80.0f ), -exp( 80.0f ), 1.0f );
	} ; 

	void set ( const float3 & m_min, const float3 & m_max ) {
		min = make_float4 (m_min);
		max = make_float4 (m_max);
	} ; 

	void update ( float4 &v0 ) {
		min = fminf1 ( min, v0 ); // fmin4
		max = fmaxf1 ( max, v0 ); // fmax4
	} ; 

	void update3 ( float4 &v0, float4 &v1, float4 &v2 ) {
		update ( v0 );
		update ( v1 );
		update ( v2 );
	} ; 

	void check_into ( AABB & box ) { // сомнительно
		min = fmaxf1 ( min, box.min ); // fmax4
		max = fminf1 ( max, box.max ); // fmin4
	} ; 

	float getVolume() {
		return ( max.x - min.x ) * ( max.y - min.y ) * ( max.z - min.z );
	} ; 

	float getSurfaceArea() {
		float side1 = max.x - min.x;
		float side2 = max.y - min.y;
		float side3 = max.z - min.z;
		// The current box has a cost of (No of triangles)*surfaceArea
		return (side1*side2 + side2*side3 + side3*side1);
	} ; 

	void print() {
		printf ( "\n min = ( %.2f, %.2f, %.2f ) ", min.x, min.y, min.z );
		printf ( "\n max = ( %.2f, %.2f, %.2f ) \n", max.x, max.y, max.z );
	} ;

	//AABB(void);
	~AABB(void);
};

