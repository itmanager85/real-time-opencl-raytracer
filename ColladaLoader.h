#pragma once

#include <vector>

//#include <iostream>
#include <string>

#include <iostream>

#include <unordered_map>

#include "./pugixml/src/pugixml.hpp"

//#include "vectors_math.h"
#include "common.h"

#include "Matrix4x4.h"

using namespace std;

static int stof_array(string s, int num_floats, float * pf) {

	if (s.length() < 1 || num_floats < 1 || NULL == pf) return 0;

	std::string::size_type sz = 0;     // alias of size_t

	for (int i = 0; i < num_floats; ++i) {
		pf[i] = std::stof(s = s.substr(sz), &sz);
	}

	return num_floats;
}

static void printf_ints(int num_ints, int * pi) {

	std::cout << "( ";

	for (int i = 0; i < num_ints - 1; ++i) {
		std::cout << pi[i] << ", ";
	}
	std::cout << pi[num_ints - 1];

	std::cout << " ) ";
}

static void printf_floats(int num_floats, float * pf) {

	std::cout << "( ";

	for (int i = 0; i < num_floats - 1; ++i) {
		std::cout << pf[i] << ", ";
	}
	std::cout << pf[num_floats - 1];

	std::cout << " ) ";
}

const static int ATTRS_NAMES_LENGTH = 10;
const static string ATTRIBUTES_NAMES[ATTRS_NAMES_LENGTH] = { "emission", "ambient", "diffuse", "specular", "shininess", "reflective", "reflectivity", "transparent", "transparency", "glossiness" };

const static int ATTRIBUTES_NUM_FLOATS[ATTRS_NAMES_LENGTH] = { 4, 4, 4, 4, 1, 4, 1, 4, 1, 1 };

const static string ATTRIBUTES_SUB_NAMES[ATTRS_NAMES_LENGTH] = { "color", "color", "color", "color", "float", "color", "float", "color", "float", "float" };

struct Effect {

	const static int PHONG = 0x001;
	const static int COOK_TORRANCE = 0x002;

	int technique; // PHONG or COOK_TORRANCE

	float4 emission;
	float4 ambient;
	float4 diffuse;
	float4 specular;
	float shininess;

	float4 reflective;
	float reflectivity;
	float4 transparent;
	float transparency;
	float glossiness;

	Effect() {
	};

	void set( int tech, float4 emi, float4 amb, float4 diff, float4 spec, float shini, float4 refl, float refl_ty, float4 transp, float transp_cy, float gloss) {

		technique = tech;

		emission = emi;
		ambient = amb;
		diffuse = diff;
		specular = spec;
		shininess = shini;

		reflective = refl;
		reflectivity = refl_ty;
		transparent = transp;
		transparency = transp_cy;

		glossiness = gloss;
	}

	void printf() {
		std::cout << "technique = " << technique << std::endl;

		std::cout << "emission = ";
		printf_floats(4, &emission.m[0]);
		std::cout << std::endl;

		std::cout << "ambient = ";
		printf_floats(4, &ambient.m[0]);
		std::cout << std::endl;

		std::cout << "diffuse = ";
		printf_floats(4, &diffuse.m[0]);
		std::cout << std::endl;

		std::cout << "specular = ";
		printf_floats(4, &specular.m[0]);
		std::cout << std::endl;

		std::cout << "shininess = " << shininess << std::endl;
		//std::cout << "shininess = ";
		//printf_floats(4, &shininess.m[0]);
		//std::cout << std::endl;

		std::cout << "reflective = ";
		printf_floats(4, &reflective.m[0]);
		std::cout << std::endl;

		std::cout << "reflectivity = " << reflectivity << std::endl;
		//std::cout << "reflectivity = ";
		//printf_floats(4, &reflectivity.m[0]);
		//std::cout << std::endl;

		std::cout << "transparent = ";
		printf_floats(4, &transparent.m[0]);
		std::cout << std::endl;

		std::cout << "transparency = " << transparency << std::endl;
		//std::cout << "transparency = ";
		//printf_floats(4, &transparency.m[0]);
		//std::cout << std::endl;

		std::cout << "glossiness = " << glossiness << std::endl;
	}
};

struct PolygonTriangle {

	int effect_index;

	int3 vertex_indices;
	int3 normal_indices;
	int3 uv0_indices;

	void printf() {
		std::cout << "effect_index = " << effect_index << std::endl;

		std::cout << "vertex_indices = ";
		printf_ints(3, &vertex_indices.m[0]);
		std::cout << std::endl;

		std::cout << "normal_indices = ";
		printf_ints(3, &normal_indices.m[0]);
		std::cout << std::endl;

		std::cout << "uv0_indices = ";
		printf_ints(3, &uv0_indices.m[0]);
		std::cout << std::endl;
	}
};

struct Geometry {
	vector <float3> float_array_positions;
	vector <float3> float_array_normals;
	vector <float2> float_array_uv0;

	vector <PolygonTriangle> polygons;
};

typedef unordered_map<string, int> Mymap;

struct SceneNode {
	int geometry_index;
	Matrix4x4 matrix;
};

class ColladaLoader
{
private:
	unordered_map<string, int> library_effect_name_to_index;
	unordered_map<string, int> library_geometry_name_to_index;

public:
	vector <Effect> library_effects;
	vector <Geometry> library_geometries;
	vector <SceneNode> library_visual_scenes;
		
public:
	
	ColladaLoader();

	bool load(const char* filename);

	~ColladaLoader();

	float3 get_vertex(float3 &v, int geometry_index);
	float3 get_normal(float3 &v, int geometry_index);

private:
	//bool load_material(FILE * file);

	bool load_effects( pugi::xml_node node_library_effects );
	bool load_effect(pugi::xml_node node_effect, int count);

	bool load_geometries(pugi::xml_node node_collada );
	bool load_geometry(pugi::xml_node node_collada, int count);

	bool load_polygons(pugi::xml_node node_polygons, Geometry &geometry);
	bool load_polygon_triangle(string s, PolygonTriangle & poly_tri);

	string get_input_source_attr(const string semantic, pugi::xml_node node_polygons);

	pugi::xml_node find_mesh_float_array(const string mesh_source_id, pugi::xml_node node_mesh);

	bool load_vertices(pugi::xml_node node_mesh_float_array, Geometry &geometry);
	bool load_normals(pugi::xml_node node_mesh_float_array, Geometry &geometry);
	bool load_tex_coords(pugi::xml_node node_mesh_float_array, Geometry &geometry);

	bool load_float_array(pugi::xml_node node_mesh_float_array, float * pf, int num_floats);

	string find_mesh_vertices_input_source_str(const string mesh_vertices_id, pugi::xml_node node_mesh);

	bool load_visual_scenes(pugi::xml_node node_collada);
	bool load_scene_node(pugi::xml_node scene_node, int count);

	Matrix4x4 load_scene_node_matrix(pugi::xml_node scene_node);

	unordered_map<string, pugi::xml_node> load_sid_to_rotate_map(pugi::xml_node scene_node);

	void compute_geometry_to_scene_index();

private:
	unordered_map <int, int> geometry_to_scene_index;

	/*
	vector<vector<float3>> arrays_of_positions;
	vector<vector<float3>> arrays_of_normals;
	vector<vector<float2>> arrays_of_uv0;

	unordered_map<string, int> arrays_of_positions_name_to_index;
	unordered_map<string, int> arrays_of_normals_name_to_index;
	unordered_map<string, int> arrays_of_uv0_name_to_index;
	*/
//public:
	//static int stof_array(string s, int num_floats, float * pf);
	//static void printf_floats(int num_floats, float * pf);
};

