#include "Camera.h"

#include <corecrt_math_defines.h>
#include <math.h>

Camera::Camera()
{
	//eye = float3( 0.0f, 0.0f, 1.0f );
	center = float3( 0.0f, 0.0f, 0.0f );
	//up = float3( 0.0f, 1.0f, 0.0f );

	cam_radius = 200.0f;
	//update_eye();

	cam_alpha = 0.0f;
	cam_beta = 0.0f;

	add_rotate(45*(3+2) * M_PI / 180.0f, 45 * M_PI / 180.0f);
}

void Camera::add_radius(float dr) {
	cam_radius += dr;
	update_eye();
}

void Camera::add_rotate(float da, float db) {

	cam_alpha += da;
	cam_beta += db;

	//int u;
	if (cam_beta < 0.0f) {
		cam_beta += 2.0f * M_PI;

	} else if (cam_beta > 2.0f * M_PI) {
		cam_beta -= 2.0f * M_PI;
	}

	if ((cam_beta > M_PI / 2.0f) && (cam_beta < 3.0f * M_PI / 2.0f)) {
		up = float3(0.0f, -1.0f, 0.0f); //u = -1;
	} else {
		up = float3(0.0f, 1.0f, 0.0f); //u = 1;
	}

	update_eye();

	//gluLookAt(x1, y1, z1, x2, y2, z2, 0, u, 0); // eye, center, up
}

void Camera::update_eye() {

	eye.x = center.x + cam_radius * cos( cam_beta ) * cos( cam_alpha );
	eye.y = center.y + cam_radius * sin( cam_beta );
	eye.z = center.z + cam_radius * cos( cam_beta ) * sin( cam_alpha );

	//gluLookAt(x1, y1, z1, x2, y2, z2, 0, u, 0); // eye, center, up

	update_full();
}

void Camera::update_full() {

	camera_direction = normalize ( center - eye );

	camera_right = normalize ( cross( camera_direction, up ) );

	camera_up = normalize(-cross(camera_direction, camera_right));
}

Camera::~Camera()
{
}
