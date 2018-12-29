#include "Matrix4x4.h"

#include <math.h>

Matrix4x4::Matrix4x4()
{
	set(1, 0, 0, 0,
		0, 1, 0, 0, 
		0, 0, 1, 0, 
		0, 0, 0, 1 );
}


Matrix4x4::~Matrix4x4()
{
}

void Matrix4x4::set( float _m00, float _m01, float _m02, float _m03,
	float _m10, float _m11, float _m12, float _m13,
	float _m20, float _m21, float _m22, float _m23,
	float _m30, float _m31, float _m32, float _m33	) {

	m00 = _m00; m01 = _m01; m02 = _m02; m03 = _m03;
	m10 = _m10; m11 = _m11; m12 = _m12; m13 = _m13;
	m20 = _m20; m21 = _m21; m22 = _m22; m23 = _m23;
	m30 = _m30; m31 = _m31; m32 = _m32; m33 = _m33;
}

void Matrix4x4::rotateX(float angle) {

	float sx = sin(angle), cx = cos(angle);

	Matrix4x4 rx;
	rx.set(1, 0, 0, 0,
		0, cx, sx, 0,
		0, -sx, cx, 0,
		0, 0, 0, 1);

	multiply(rx);
}

void Matrix4x4::rotateY(float angle) {
	float sy = sin(angle), cy = cos(angle);

	Matrix4x4 ry;
	ry.set(cy, 0, -sy, 0,
		0, 1, 0, 0,
		sy, 0, cy, 0,
		0, 0, 0, 1);

	multiply(ry);
}

void Matrix4x4::rotateZ(float angle) {
	float sz = sin(angle), cz = cos(angle);

	Matrix4x4 rz;
	rz.set(cz, sz, 0, 0,
		-sz, cz, 0, 0,
		0, 0, 1, 0,
		0, 0, 0, 1);

	multiply(rz);
}

void Matrix4x4::rotateXYZ(float rot_x, float rot_y, float rot_z) {
	rotateX(rot_x);
	rotateY(rot_y);
	rotateZ(rot_z);
}

void Matrix4x4::translate(float tx, float ty, float tz) {

	Matrix4x4 rt;

	rt.set(1, 0, 0, 0,
		0, 1, 0, 0,
		0, 0, 1, 0,
		tx, ty, tz, 1);

	multiply(rt);
}

void Matrix4x4::multiply(Matrix4x4& mat) {

	Matrix4x4 result;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {

			result.m[i * 4 + j] = 0;

			for (int k = 0; k < 4; ++k) {
				//m[i * 4 + j] * mat.m[]

				result.m[i * 4 + j] += m[i * 4 + k] * mat.m[k * 4 + j];
			}
		}
	}

	*this = result;
}

void Matrix4x4::transponse() {

	Matrix4x4 result;

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			//result.m[i * 4 + j] = m[j*4+i];
			result.m[j * 4 + i] = m[i * 4 + j];
		}
	}

	*this = result;
}