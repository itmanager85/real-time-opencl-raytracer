#pragma once

#include "vectors_math.h"

class Camera
{
public:
	float3 eye, center, up;

	float cam_radius, cam_alpha, cam_beta;

	float3 camera_right;
	float3 camera_up;
	float3 camera_direction;

public:
	Camera();

	void add_radius(float dr);
	void add_rotate( float da, float db );

	~Camera();

private:
	void update_eye();

	void update_full();
};

