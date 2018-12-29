#pragma once

#include "vectors_math.h"

class Matrix4x4
{
public:

	union
	{
		float m[16];

		struct {
			float	m00, m01, m02, m03,
					m10, m11, m12, m13,
					m20, m21, m22, m23,
					m30, m31, m32, m33;
		};
	};

public:
	Matrix4x4();
	~Matrix4x4();

	void set( float	_m00, float	_m01, float	_m02, float	_m03,
		float _m10, float _m11, float _m12, float _m13,
		float _m20, float _m21, float _m22, float _m23,
		float _m30, float _m31, float _m32, float _m33
	);

	void rotateX(float angle);
	void rotateY(float angle);
	void rotateZ(float angle);

	void rotateXYZ(float rot_x, float rot_y, float rot_z);

	void translate(float tx, float ty, float tz);

	void multiply(Matrix4x4& mat);

	void transponse();
};

static float4 multiply(float4 &v, Matrix4x4 & mat) {

	//float4 v( vert.x, vert.y, vert.z, 1.0f);

	float4 result;

	for (int i = 0; i < 4; i++) {
		float s = 0.0f;
		for (int j = 0; j < 4; j++) {
			s += v.m[j] * mat.m[j * 4 + i];
		}
		result.m[i] = s;
	}

	//return float3(result.x, result.y, result.z);
	return result;
}; 

