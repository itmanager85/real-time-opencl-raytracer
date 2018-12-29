#pragma once

//#include "vectors_math.h"
#include "common.h"

#include <vector>

#include "ColladaLoader.h"

using namespace std;

/*struct TriFace {
	union {
		struct { int i0, i1, i2; };
		int m[3];
	};

} ; */

struct Material {
	int4 technique; // int // PHONG or COOK_TORRANCE 

	float4 emission;
	float4 ambient;
	float4 diffuse;
	float4 specular;
	float4 shininess; // float

	float4 reflective;
	float4 reflectivity; // float
	float4 transparent;
	float4 transparency; // float
	float4 glossiness; // float

	Material() {
		// by default
		set(Effect::PHONG, float4(0, 0, 0, 1.0f), float4(0.0f, 0.0f, 0.0f, 1.0f), float4(1, 1, 1, 1.0f), 
			float4(1, 1, 1, 1.0f), 2, float4(1, 1, 1, 1.0f), 1, float4(1, 1, 1, 1), 0, 1.0f);
	}

	Material(Effect & effect) {
		set(effect);
	}

	void set(int tech, float4 emi, float4 amb, float4 diff, float4 spec, float shini, float4 refl, float refl_ty, float4 transp, float transp_cy, float gloss) {

		technique.x = tech;

		emission = emi;
		ambient = amb;
		diffuse = diff;
		specular = spec;
		shininess.x = shini;

		reflective = refl;
		reflectivity.x = refl_ty;
		transparent = transp;
		transparency.x = transp_cy;

		glossiness.x = gloss;
	}

	void set(Effect & effect) {
		set(effect.technique, effect.emission, effect.ambient, effect.diffuse, effect.specular, effect.shininess,
			effect.reflective, effect.reflectivity, effect.transparent, effect.transparency, effect.glossiness);
	}
};

class Mesh
{
public:
	vector <int>	indices;
	vector <float4> vertices;

	// for smothing triangle reflections
	vector <int>	normals_indices;
	vector <float4> normals; 

	vector <Material> materials;

	vector <int> triangle_index_to_material_index;

	float3 scene_aabbox_min;
	float3 scene_aabbox_max;

public:
	Mesh(void);

	void init( ColladaLoader & collada_loader );

	void init ( TriangleMesh & triMesh );

	void add_tri ( float4& v0, float4& v1, float4& v2 );

	int* getTrisPtr() { return &indices[0]; };
	float4* getVertsPtr() { return &vertices[0]; };

	int getNumTriangles() { return indices.size()/3; };

	~Mesh(void);
};

