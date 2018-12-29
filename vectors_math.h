#pragma once

const int AXIS_X = 0, AXIS_Y = 1, AXIS_Z = 2;

float rsqrtf(float x);

//struct float3;

struct int2 {
	union {
		struct { int x, y; };
		int m[2];
	};

	int2() {
		x = 0; y = 0;
	};

	int2(int x_, int y_) {
		x = x_; y = y_;
	};
};

struct int3 {
	union {
		struct { int x,y,z ; };
		int m[3];
	};

	int3 ( ) {
		x = 0; y = 0; z = 0; 
	} ;

	int3 ( int x_, int y_, int z_ ) {
		x = x_; y = y_; z = z_; 
	} ; 

	//int3 ( float3 &b ) {
	//	x = b.x; y = b.y; z = b.z; 
	//} ; 

	const int&        get         (int idx) const             { /*FW_ASSERT(idx >= 0 && idx < L);*/ return m[idx]; }
    int&              get         (int idx)                   { /*FW_ASSERT(idx >= 0 && idx < L);*/ return m[idx]; }
} ;

struct int4 {
	union {
		struct { int x, y, z, w; };
		int m[4];
	};

	int4() {
		x = 0; y = 0; z = 0; w = 0;
	};

	int4(int x_, int y_, int z_, int w_) {
		x = x_; y = y_; z = z_; w = w_;
	};
};

struct float2 {
	union {
		struct { float x,y ; };
		float m[2];
	};

	float2 ( ) {
		x = 0; y = 0;
	} ;

	float2 ( float x_, float y_ ) {
		x = x_; y = y_;
	} ; 
} ; 

struct float3 {
	union {
		struct { float x,y,z ; };
		float m[3];
	};
	
	float3 ( ) {
		x = 0; y = 0; z = 0; 
	} ;

	float3 ( float x_, float y_, float z_ ) {
		x = x_; y = y_; z = z_; 
	} ; 

	const float&        get         (int idx) const             { /*FW_ASSERT(idx >= 0 && idx < L);*/ return m[idx]; }
    float&              get         (int idx)                   { /*FW_ASSERT(idx >= 0 && idx < L);*/ return m[idx]; }
} ; 

struct float4 {

	union {
		struct { float x,y,z,w ; };
		float m[4];
	};

	float4 ( ) {
		x = 0; y = 0; z = 0; w = 0; 
	} ;

	float4 ( float x_, float y_, float z_, float w_ ) {
		x = x_; y = y_; z = z_; w = w_; 
	} ; 

	const float&        get         (int idx) const             { /*FW_ASSERT(idx >= 0 && idx < L);*/ return m[idx]; }
    float&              get         (int idx)                   { /*FW_ASSERT(idx >= 0 && idx < L);*/ return m[idx]; }
} ; 

	float2 make_float2 ( float x, float y );

	float3 operator/(float3 &a, float3 &b);

	float4 fminf1 ( const float4 &v0, const float4 &v1 );
	float4 fmaxf1 ( const float4 &v0, const float4 &v1 );

	float3 fminf1 ( const float3 &v0, const float3 &v1 );
	float3 fmaxf1 ( const float3 &v0, const float3 &v1 );

	float4 operator-(float4 &a, float4 &b);

	float3 make_float3 ( float x, float y, float z );
	float3 operator-(float3 &a);
	float3 operator+(float3 a, float3 b);
	void operator*=(float3 &a, float b);
	float3 operator*(float b, float3 a);
	float dot(float3 a, float3 b);
	float3 cross(float3 a, float3 b);
	float3 normalize(float3 v);

	float4 make_float4 ( float x, float y, float z, float w );
	float4 operator-(float4 &a);
	float4 operator+(float4 a, float4 b);
	void operator*=(float4 &a, float b);
	float4 operator*(float b, float4 a);
	float dot(float4 a, float4 b);
	float4 normalize(float4 v);

	float imin1(const int a, const int b);
	float imax1(const int a, const int b);

	float fminf1(const float a, const float b);
	float fmaxf1(const float a, const float b);
	float clamp(const float f, const float a, const float b);

float3 operator-(float3 a, float3 b);
float3 operator+(float3 a, float b);
void operator+=(float3 &a, float b);
	float3 operator*(float3 a, float b);
	float length(float3 v);
	void operator+=(float3 &a, float3 b);
	float3 operator/(float3 a, float b);
	void operator/=(float3 &a, float3 b);
	void operator/=(float3 &a, float b);

	float3 get_float3( float4 &v );

	float3 operator/(float b, float3 a);
	float4 operator*(float4 a, float b);
	float4 lerp(const float4& a, const float4& b, const float t);

	float fminf1(float a, float b, float c);
	float fmaxf1(float a, float b, float c);

	int3 make_int3 ( int x, int y, int z );
	int3 operator-(int3 &a, int3 &b);
	int3 operator*(int3 a, float3 b);

	float3 operator*(float3 a, float3 b);

	int3 min3i ( const int3 &v0, const int a );
	int3 max3i ( const int3 &v0, const int a );

	int3 clamp(const int3 i3, const int a, const int b);

	int3 max3i ( const int3 &v0, const int3 &a );

	int3 clamp(const int3 i3, const int3 &a, const int b);

	int3 make_int3 ( const float3 & a );

	float lerp(const float a, const float b, const float t);

	float fminf1 ( const float3& a );
	float fmaxf1 ( const float3& a );
	float fsumf ( const float3& a );

	float3 make_float3 ( const float4& v );

	float4 make_float4 ( const float3& v );

	float4& operator+=(float4 &a, float4& b);