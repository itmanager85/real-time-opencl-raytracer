#include "Mesh.h"


Mesh::Mesh(void)
{
	//scene_aabbox_min = float3(0, 0, 0);
	//scene_aabbox_max = float3(0, 0, 0);
}

void Mesh::init(ColladaLoader & collada_loader) {

	//scene_aabbox_min = float3(0, 0, 0);
	//scene_aabbox_max = float3(0, 0, 0);

	for (int i = 0; i < collada_loader.library_effects.size(); ++i) {
		materials.push_back( Material(collada_loader.library_effects[i]) );
	}

	//vertices.resize( );

	int triangle_index_count = 0;
	int vertex_count = 0;
	int normal_count = 0;

	for (int i = 0; i < collada_loader.library_geometries.size(); ++i) {

		Geometry &geometry = collada_loader.library_geometries[i];

		for (int j = 0; j < geometry.polygons.size(); ++j) {

			PolygonTriangle &p = geometry.polygons[j];

			indices.push_back(vertex_count + p.vertex_indices.x);
			indices.push_back(vertex_count + p.vertex_indices.y);
			indices.push_back(vertex_count + p.vertex_indices.z);

			normals_indices.push_back(normal_count + p.normal_indices.x);
			normals_indices.push_back(normal_count + p.normal_indices.y);
			normals_indices.push_back(normal_count + p.normal_indices.z);

			//triangle_index_to_material[ triangle_index_count ] = p.effect_index;
			triangle_index_to_material_index.push_back(p.effect_index);

			++triangle_index_count;
		}

		for (int j = 0; j < geometry.float_array_positions.size(); ++j) {


			float3 v = collada_loader.get_vertex( geometry.float_array_positions[j], i);

			//float3 & v = geometry.float_array_positions[j];
			//float4 v4 = collada_loader.get_vert( float4(v.x, v.y, v.z, 1.0f), i);

			if (vertex_count < 1) {
				scene_aabbox_min = v;
				scene_aabbox_max = v;
			} else {
				scene_aabbox_min = fminf1(scene_aabbox_min, v);
				scene_aabbox_max = fmaxf1(scene_aabbox_max, v);
			}

			vertices.push_back(float4(v.x, v.y, v.z, 1.0f));

			++vertex_count;
		}

		for (int j = 0; j < geometry.float_array_normals.size(); ++j) {

			//float3 & n = geometry.float_array_normals[j];
			float3 n = collada_loader.get_normal(geometry.float_array_normals[j], i);

			normals.push_back(float4(n.x, n.y, n.z, 1.0f));

			++normal_count;
		}
	}
}

void Mesh::init ( TriangleMesh & triMesh ) {

	indices.resize( triMesh.faces.size() * 3 );
	vertices.resize( triMesh.verts.size() );

	// for smothing triangle reflections
	normals_indices.resize(triMesh.faces.size() * 3);
	normals.resize( triMesh.norms.size() ); 

	for ( int i=0; i<triMesh.faces.size(); ++i ) {
		indices[ i*3 + 0 ] = triMesh.faces[i].v[0]-1;
		indices[ i*3 + 1 ] = triMesh.faces[i].v[1]-1;
		indices[ i*3 + 2 ] = triMesh.faces[i].v[2]-1;
	} ; 

	//float3 mid3 = ( triMesh.bounding_box[0] + triMesh.bounding_box[1] ) / 2;
	//float4 mid4 = make_float4 ( mid3.x, mid3.y, mid3.z, 1.0f );

	for ( int i=0; i<triMesh.verts.size(); ++i ) {
		vertices[i].x = triMesh.verts[i].x;
		vertices[i].y = triMesh.verts[i].y;
		vertices[i].z = triMesh.verts[i].z;

	//	vertices[i] += mid4;

		if (i < 1) {
			scene_aabbox_min = triMesh.verts[i];
			scene_aabbox_max = triMesh.verts[i];
		}
		else {
			scene_aabbox_min = fminf1(scene_aabbox_min, triMesh.verts[i]);
			scene_aabbox_max = fmaxf1(scene_aabbox_max, triMesh.verts[i]);
		}
	} ; 

	// for smothing triangle reflections
	for (int i = 0; i<triMesh.faces.size(); ++i) {
		normals_indices[i * 3 + 0] = triMesh.faces[i].vn[0] - 1;
		normals_indices[i * 3 + 1] = triMesh.faces[i].vn[1] - 1;
		normals_indices[i * 3 + 2] = triMesh.faces[i].vn[2] - 1;
	};

	
	for (int i = 0; i<triMesh.norms.size(); ++i) { 
		normals[i].x = triMesh.norms[i].x;
		normals[i].y = triMesh.norms[i].y;
		normals[i].z = triMesh.norms[i].z;

		normals[i] = normalize( normals[i] );
	};
} 

void Mesh::add_tri ( float4& v0, float4& v1, float4& v2 ) {

	indices.push_back( vertices.size() + 0 );
	indices.push_back( vertices.size() + 1 );
	indices.push_back( vertices.size() + 2 );

	vertices.push_back( v0 );
	vertices.push_back( v1 );
	vertices.push_back( v2 );
} 

Mesh::~Mesh(void)
{
}
