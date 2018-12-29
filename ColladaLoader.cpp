#include "ColladaLoader.h"

//#include "./pugixml/src/pugixml.hpp"

#include <iostream>

#include <corecrt_math_defines.h>

ColladaLoader::ColladaLoader()
{
}

bool ColladaLoader::load(const char* filename) {

	//filename = "data/collada/tree.xml";

	pugi::xml_document doc;

	pugi::xml_parse_result result = doc.load_file(filename);

	//std::cout << "Load result: " << result.description() << ", mesh name: " << doc.child("mesh").attribute("name").value() << std::endl;

	std::cout << "Load result: " << result.description() << ", collada: " << doc.child("COLLADA").attribute("xmlns").value() << std::endl; // GOOD !!

	//std::cout << "effect_id = " << doc.child("COLLADA").child("library_effects").child("effect").attribute("id").value() << std::endl; // GOOD !!
	
	pugi::xml_node node_COLLADA = doc.child("COLLADA");

	pugi::xml_node node_library_effects = node_COLLADA.child("library_effects");

	//std::cout << "effect_id = " << library_effects.first_child().attribute("id").value() << std::endl; 
	//std::cout << "effect_id = " << library_effects.child("effect").attribute("id").value() << std::endl; 

	/*
	pugi::xml_node effect = library_effects.child("effect");
	std::cout << "effect_id = " << effect.attribute("id").value() << std::endl;

	effect = effect.next_sibling("effect");
	std::cout << "effect_id = " << effect.attribute("id").value() << std::endl;

	effect = effect.next_sibling("effect");
	std::cout << "effect_id3 = " << effect.attribute("id").value() << std::endl;
	*/

	//pugi::xml_node effect = library_effects_node.child("effect");

	load_effects(node_library_effects);

	//std::cout << "geometry_id = " << node_COLLADA.child("library_geometries").child("geometry").attribute("id").as_string() << std::endl;
	//std::cout << "polygons_material = " << node_COLLADA.child("library_geometries").child("geometry").child("mesh").child("polygons").attribute("material").as_string() << std::endl;
	//std::cout << "node_polygons = " << std::endl;

	//Geometry geometry;
	//load_polygons (node_COLLADA.child("library_geometries").child("geometry").child("mesh").child("polygons"), geometry );

	load_geometries( node_COLLADA );

	load_visual_scenes(node_COLLADA);

	compute_geometry_to_scene_index(); // 

//	getchar();

	/*
	FILE * file = fopen(filename, "r");
	if (!file) {
	
		printf("Impossible to open the file %s !\n", filename );
		return false;
	}

	while (true) {
	
		char line_header[256];

		if (EOF == fscanf(file, "%s", line_header))
			break;

		if (0 == strcmp(line_header, "<library_effects>")) {

			while (true) {

				char line_header[256];

				if (EOF == fscanf(file, "%s", line_header))
					break;

				if (0 == strcmp(line_header, "</library_effects>")) break;

				else if (0 == strcmp(line_header, "<effect")) {

				}
				
			}
		
		}
	}
	
	fclose( file );
	*/

	return true; 
};


ColladaLoader::~ColladaLoader()
{
}

//bool ColladaLoader::load_material( FILE * file ) {
//}

bool ColladaLoader::load_effect(pugi::xml_node node_effect, int count) {

	library_effect_name_to_index[node_effect.attribute("name").value()] = count;

	pugi::xml_node node_technique = node_effect.child("profile_COMMON").child("technique");

	//std::cout << "technique_sid = " << node_technique.attribute("sid").as_string() << std::endl;
	//getchar();

	Effect effect;

	pugi::xml_node node_technique_current;

	if (node_technique_current = node_technique.child("cook-torrance")) {
		effect.technique = Effect::COOK_TORRANCE;
	}
	else if (node_technique_current = node_technique.child("phong")) {
		effect.technique = Effect::PHONG;
	}
	else {
		std::cout << "effect_technique_not_found = " << node_effect.attribute("name").value() << std::endl;

		return false;
	}

	//std::cout << "effect_technique = " << node_technique_current.name() << std::endl;
	//std::cout << "effect_technique = " << effect.technique << std::endl;

	//pugi::xml_text txt = node_technique_current.child("emission").child("color").text();
	//std::cout << "emission_color = " << txt.as_string() << std::endl;
	//getchar();

	float4 colors[ATTRS_NAMES_LENGTH];
	
	//for (int i = 0; i < 1; ++i) {
	for (int i = 0; i < ATTRS_NAMES_LENGTH; ++i) {

		//string txt_floats = "0.1 0.2 0.3 1.0";
		string txt_floats = node_technique_current.child( ATTRIBUTES_NAMES[i].c_str() ).child( ATTRIBUTES_SUB_NAMES[i].c_str() ).text().as_string();

//		std::cout << ATTRIBUTES_NAMES[i] << "_" << ATTRIBUTES_SUB_NAMES[i] << " = " << txt_floats << std::endl;
		 
		stof_array(txt_floats, ATTRIBUTES_NUM_FLOATS[i], &colors[i].m[0]);

		//printf_floats(4, &colors[i].m[0]);

		//getchar();
	}

	effect.set( effect.technique, colors[0], colors[1], colors[2], colors[3], colors[4].x, colors[5], colors[6].x, colors[7], colors[8].x, colors[9].x);
	
	//effect.emission = float4(0.1f, 0.2f, 0.3f, 1.0f);
	
	//float4 test = float4(0.1f, 0.2f, 0.3f, 1.0f);
	//effect.set(effect.technique, test, test, test, test, test, test, test, test, test);

//	effect.printf();

//	getchar();

	library_effects.push_back( effect );

	return true;
}

bool ColladaLoader::load_effects(pugi::xml_node node_library_effects) {
	
	if (!node_library_effects) return false;
	
	int count = 0;
	for (pugi::xml_node node_effect = node_library_effects.child("effect"); node_effect; node_effect = node_effect.next_sibling("effect"), ++count) {
		//std::cout << "effect_id" << count << " = " << effect.attribute("id").value() << std::endl;
		load_effect(node_effect, count);
	}

	//for (int i = 0; i < count; ++i) {
	//std::cout << "effect_name #" << count << " = " << effect.attribute("id").value() << std::endl;
	//}

	/*
	for (Mymap::const_iterator it = library_effect_name_to_index.begin(); it != library_effect_name_to_index.end(); ++it) {
		std::cout << " [" << it->first << ", " << it->second << "]";
		std::cout << std::endl;
	}
	*/

	return true;
}

bool ColladaLoader::load_geometries(pugi::xml_node node_collada) {

	if (!node_collada) return false;

	pugi::xml_node node_library_geometries = node_collada.child("library_geometries");

	if (!node_library_geometries) return false;

	int count = 0;
	for (pugi::xml_node node_geometry = node_library_geometries.child("geometry"); node_geometry; node_geometry = node_geometry.next_sibling("geometry"), ++count) {
		//std::cout << "effect_id" << count << " = " << effect.attribute("id").value() << std::endl;
		load_geometry(node_geometry, count);
	}

	return true;
}

bool ColladaLoader::load_geometry(pugi::xml_node node_geometry, int count) {

	pugi::xml_node node_polygons = node_geometry.child("mesh").child("polygons");

	library_geometry_name_to_index[node_geometry.attribute("id").as_string()] = count;

	Geometry geometry;

	load_polygons(node_polygons, geometry);

	pugi::xml_node node_mesh = node_polygons.parent();

	//cout << "imput_source = " << node_polygons.child("input").attribute("source").as_string() << endl;

	string node_name_vertex = get_input_source_attr("VERTEX", node_polygons).substr(1);
	string node_name_normal = get_input_source_attr("NORMAL", node_polygons).substr(1);
	string node_name_uv0 = get_input_source_attr("TEXCOORD", node_polygons).substr(1);

//	cout << "imput_VERTEX_source = " << node_name_vertex << endl;
//	cout << "imput_NORMAL_source = " << node_name_normal << endl;
//	cout << "imput_TEXCOORD_source = " << node_name_uv0 << endl;

	string node_name_vertex_2 = find_mesh_vertices_input_source_str(node_name_vertex, node_mesh);
//	cout << "imput_VERTEX_2_source = " << node_name_vertex_2 << endl;

	pugi::xml_node node_mesh_float_array_vertex2 = find_mesh_float_array(node_name_vertex_2, node_mesh);
	pugi::xml_node node_mesh_float_array_normal = find_mesh_float_array(node_name_normal, node_mesh);
	pugi::xml_node node_mesh_float_array_uv0 = find_mesh_float_array(node_name_uv0, node_mesh);

//	cout << "node_mesh_float_array_VERTEX_2_source_id= " << node_mesh_float_array_vertex2.attribute("id").as_string() << endl;
//	cout << "node_mesh_float_array_NORMAL_source_id= " << node_mesh_float_array_normal.attribute("id").as_string() << endl;
//	cout << "node_mesh_float_array_TEXCOORD_source_id= " << node_mesh_float_array_uv0.attribute("id").as_string() << endl;

	load_vertices(node_mesh_float_array_vertex2, geometry);
	load_normals(node_mesh_float_array_normal, geometry);
	load_tex_coords(node_mesh_float_array_uv0, geometry);

	//printf_floats(geometry.float_array_uv0.size() * 2, &geometry.float_array_uv0[0].m[0]);

//	getchar();

	//pugi::xml_node node_uv0 = find_node_array(node_mesh, );

	//load_texture_coords();

	library_geometries.push_back(geometry);

	return true;
}

bool ColladaLoader::load_polygons(pugi::xml_node node_polygons, Geometry &geometry) {
	//cout << "node_polygons_p = " << node_polygons.child("p").text().as_string() << std::endl;

	int num_polygons;
	sscanf(node_polygons.attribute("count").value(), "%d", &num_polygons);	

	//std::cout << "num_polygons = " << num_polygons << std::endl;

	string material = node_polygons.attribute("material").as_string();

//	std::cout << "material = " << material << std::endl;

	int effect_index = library_effect_name_to_index[ material ];

//	std::cout << "effect_index = " << effect_index << std::endl;

//	getchar();

	geometry.polygons.resize(num_polygons);

	pugi::xml_node node_polygons_p = node_polygons.child("p");
	for (int i = 0; i < num_polygons; ++i) {

		geometry.polygons[i].effect_index = effect_index;

		//cout << "node_polygons_p = " << node_polygons_p.text().as_string() << std::endl;

		load_polygon_triangle( node_polygons_p.text().as_string(), geometry.polygons[i]);

		//geometry.polygons[i].printf();

		node_polygons_p = node_polygons_p.next_sibling("p");
		
		//getchar();
	}

	return true;
}



bool ColladaLoader::load_vertices(pugi::xml_node node_mesh_float_array, Geometry &geometry) {

	int num_floats = node_mesh_float_array.attribute("count").as_int();
//	cout << "vertices_float_array_count = " << num_floats << endl;

	int num_verts = num_floats / 3;
	
	geometry.float_array_positions.resize(num_verts);

	load_float_array(node_mesh_float_array, &geometry.float_array_positions[0].m[0], num_floats);

	return true;
}

bool ColladaLoader::load_normals(pugi::xml_node node_mesh_float_array, Geometry &geometry) {

	int num_floats = node_mesh_float_array.attribute("count").as_int();
//	cout << "normals_float_array_count = " << num_floats << endl;

	int num_verts = num_floats / 3;

	geometry.float_array_normals.resize(num_verts);

	load_float_array(node_mesh_float_array, &geometry.float_array_normals[0].m[0], num_floats);

	return true;
}

bool ColladaLoader::load_tex_coords(pugi::xml_node node_mesh_float_array, Geometry &geometry) {

	int num_floats = node_mesh_float_array.attribute("count").as_int();
//	cout << "uv0_float_array_count = " << num_floats << endl;

	int num_verts = num_floats / 2;

	geometry.float_array_uv0.resize(num_verts);

	load_float_array(node_mesh_float_array, &geometry.float_array_uv0[0].m[0], num_floats);

	return true;
}

bool ColladaLoader::load_float_array(pugi::xml_node node_mesh_float_array, float * pf, int num_floats) {
	stof_array(node_mesh_float_array.text().as_string(), num_floats, pf);
	return true;
}

string ColladaLoader::find_mesh_vertices_input_source_str(const string mesh_vertices_id, pugi::xml_node node_mesh) {

	if (!node_mesh) return "";

	int count = 0;
	for (pugi::xml_node node_mesh_vertices = node_mesh.child("vertices"); node_mesh_vertices; node_mesh_vertices = node_mesh_vertices.next_sibling("vertices"), ++count) {
		//std::cout << "effect_id" << count << " = " << effect.attribute("id").value() << std::endl;

		//cout << "mesh_source_id = " << node_mesh_source.attribute("id").as_string() << endl;
		//getchar();

		if (mesh_vertices_id == node_mesh_vertices.attribute("id").as_string()) {

			//cout << "found, mesh_source_id = " << node_mesh_source.attribute("id").as_string() << endl;
			//getchar();

			return string( node_mesh_vertices.child("input").attribute("source").as_string() ).substr(1);
		}
	}

	return "";
}

pugi::xml_node ColladaLoader::find_mesh_float_array(const string mesh_source_id, pugi::xml_node node_mesh ) {

	if (!node_mesh) return pugi::xml_node();

	int count = 0;
	for (pugi::xml_node node_mesh_source = node_mesh.child("source"); node_mesh_source; node_mesh_source = node_mesh_source.next_sibling("source"), ++count) {
		//std::cout << "effect_id" << count << " = " << effect.attribute("id").value() << std::endl;

		//cout << "mesh_source_id = " << node_mesh_source.attribute("id").as_string() << endl;
		//getchar();

		if ( mesh_source_id == node_mesh_source.attribute( "id" ).as_string() ) {

			//cout << "found, mesh_source_id = " << node_mesh_source.attribute("id").as_string() << endl;
			//getchar();

			return node_mesh_source.child("float_array");
		}
	}

	return pugi::xml_node();
}

string ColladaLoader::get_input_source_attr( const string semantic, pugi::xml_node node_polygons ) {

	//const static string 

	pugi::xml_node node_polygons_input = node_polygons.child("input");

	for(int i = 0; i < 3; ++i) {
		if (semantic == node_polygons_input.attribute("semantic").as_string()) {
			return  node_polygons_input.attribute("source").as_string();
		}

		node_polygons_input = node_polygons_input.next_sibling("input");
	}

	return "";
}

bool ColladaLoader::load_polygon_triangle( string s, PolygonTriangle & poly_tri) {
	
	sscanf( s.c_str(), "%d %d %d %d %d %d %d %d %d", 
		&poly_tri.vertex_indices.x, &poly_tri.normal_indices.x, &poly_tri.uv0_indices.x,
		&poly_tri.vertex_indices.y, &poly_tri.normal_indices.y, &poly_tri.uv0_indices.y,
		&poly_tri.vertex_indices.z, &poly_tri.normal_indices.z, &poly_tri.uv0_indices.z );

	return true;
}

bool ColladaLoader::load_visual_scenes(pugi::xml_node node_collada) {
	if (!node_collada) return false;

	pugi::xml_node scene_nodes = node_collada.child("library_visual_scenes").child("visual_scene");

	if (!scene_nodes) return false;

	int count = 0;
	for (pugi::xml_node scene_node = scene_nodes.child("node"); scene_node; scene_node = scene_node.next_sibling("node"), ++count) {
//		std::cout << "count = " << count << ", scene_node _id = " << scene_node.attribute("id").as_string() << std::endl;
		load_scene_node(scene_node, count);
	}

	return true;
}

bool ColladaLoader::load_scene_node(pugi::xml_node scene_node, int count) {

	if (!scene_node) return false;

	string geometry_id = string( scene_node.child("instance_geometry").attribute("url").as_string() ).substr(1);
	int geometry_id_index = library_geometry_name_to_index[geometry_id];

//	std::cout << "geometry_id" << geometry_id << std::endl;
//	std::cout << "geometry_id_index = " << geometry_id_index << std::endl;

	SceneNode node;

	node.geometry_index = geometry_id_index;

	node.matrix = load_scene_node_matrix(scene_node);

	library_visual_scenes.push_back(node);

//	getchar();

	return true;
}

Matrix4x4 ColladaLoader::load_scene_node_matrix(pugi::xml_node scene_node) {

	Matrix4x4 matrix;

	if (!scene_node) return matrix;	

	pugi::xml_node node_matrix = scene_node.child("matrix");

	if (node_matrix) {
		string node_matrix_str = scene_node.child("matrix").text().as_string();
		stof_array(node_matrix_str, 4*4, &matrix.m[0] );

		//std::cout << "node_matrix = ";
		//printf_floats( 4*4, &matrix.m[0] );
		//getchar();

		matrix.transponse();

		return matrix;
	}

	unordered_map<string, pugi::xml_node> sid_to_node_rotate_map = load_sid_to_rotate_map (scene_node);

	const static int NUM_ROTATE_SIDS = 6;
	//const static string ROTATE_SID_ARRAY[NUM_ROTATE_SIDS] = { "jointOrientX", "jointOrientY", "jointOrientZ", "rotateX", "rotateY", "rotateZ" };

	//const static string ROTATE_SID_ARRAY[NUM_ROTATE_SIDS] = { "jointOrientZ", "jointOrientY", "jointOrientX", "rotateZ", "rotateY", "rotateX" };
	//const static string ROTATE_SID_ARRAY[NUM_ROTATE_SIDS] = { "jointOrientX", "jointOrientY", "jointOrientZ", "rotateZ", "rotateY", "rotateX" };

	//const static string ROTATE_SID_ARRAY[NUM_ROTATE_SIDS] = { "jointOrientX", "jointOrientY", "jointOrientZ", "rotateY", "rotateZ", "rotateX" }; // 
	//const static string ROTATE_SID_ARRAY[NUM_ROTATE_SIDS] = { "jointOrientX", "jointOrientY", "jointOrientZ", "rotateY", "rotateX", "rotateZ" };

	const static string ROTATE_SID_ARRAY[NUM_ROTATE_SIDS] = { "jointOrientX", "jointOrientY", "jointOrientZ", "rotateX", "rotateZ", "rotateY" };

	for (int i = 0; i < NUM_ROTATE_SIDS; i++) {
		
		pugi::xml_node rotate_node = sid_to_node_rotate_map[ ROTATE_SID_ARRAY[i] ];

		//pugi::xml_node rotate_node = sid_to_node_rotate_map["jointOrientX"];

		if (!rotate_node) {
			//std::cout << "not_found : ROTATE_SID_ARRAY[" << i << "] = " << ROTATE_SID_ARRAY[i] << std::endl;
			//getchar();
			continue;
		}

//		std::cout << "ROTATE_SID_ARRAY[" << i << "] = " << ROTATE_SID_ARRAY[i] << std::endl;
//		std::cout << "node_rotate_sid = " << rotate_node.attribute("sid").as_string() << std::endl;
		
		string angle_str = string(rotate_node.text().as_string()).substr(6);

		float angle;
		stof_array(angle_str, 1, &angle);

//		std::cout << "angle" << angle << std::endl;

		float rad_angle = angle * M_PI / 180.0f;

		switch ( i%3 ) {

			case 0: matrix.rotateX(rad_angle);
				break;

			case 1: matrix.rotateY(rad_angle);
				break;

			case 2: matrix.rotateZ(rad_angle);
				//matrix.rotateZ(45 * M_PI / 180);
				break;

			default: 
				break;
		}
	}

	string translate_str = string(scene_node.child("translate").text().as_string());

	float translate[3];

	stof_array(translate_str, 3, &translate[0]);

//	std::cout << "translate = ";
//	printf_floats( 3, &translate[0]);

	matrix.translate(translate[0], translate[1], translate[2]);

	return matrix;
}

unordered_map<string, pugi::xml_node> ColladaLoader::load_sid_to_rotate_map(pugi::xml_node scene_node) {

	unordered_map<string, pugi::xml_node> sid_to_node_rotate_map;

	if (!scene_node) return sid_to_node_rotate_map;

	int count = 0;
	for (pugi::xml_node rotate_node = scene_node.child("rotate"); rotate_node; rotate_node = rotate_node.next_sibling("rotate"), ++count) {
		//std::cout << "count = " << count << ", scene_node _id = " << scene_node.attribute("id").as_string() << std::endl;

		sid_to_node_rotate_map[rotate_node.attribute("sid").as_string()] = rotate_node;
	}

	return sid_to_node_rotate_map;
}

float3 ColladaLoader::get_vertex(float3 &v, int geometry_index) {
	float4 result = multiply(float4(v.x, v.y, v.z, 1.0f), library_visual_scenes[geometry_to_scene_index[geometry_index]].matrix);
	return float3(result.x, result.y, result.z);
}

float3 ColladaLoader::get_normal(float3 &n, int geometry_index) {
	float4 result = multiply(float4(n.x, n.y, n.z, 0.0f), library_visual_scenes[geometry_to_scene_index[geometry_index]].matrix);
	return float3(result.x, result.y, result.z);
}

void ColladaLoader::compute_geometry_to_scene_index() {
	for (int i = 0; i < library_geometries.size(); ++i) {

		//std::cout << "library_visual_scenes.size() = " << library_visual_scenes.size() << endl;

		//std::cout << "geometry_index = " << library_visual_scenes[i].geometry_index << endl;
		geometry_to_scene_index[library_visual_scenes[i].geometry_index] = i;
		//getchar();
	}
}