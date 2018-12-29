/*
 * Copyright 1993-2010 NVIDIA Corporation.  All rights reserved.
 *
 * Please refer to the NVIDIA end user license agreement (EULA) associated
 * with this source code for terms and conditions that govern your use of
 * this software. Any use, reproduction, disclosure, or distribution of
 * this software and related documentation outside the terms of the EULA
 * is strictly prohibited.
 *
 */

 #define RAY_TRACE_DEPTH 3 // 3 // 5 //  10 // 2 // 
 #define AMBIENT_LIGHT

#define maxSteps 500
#define tstep 0.01f

struct Material_pbr {
    float3 albedo;
    float3 f0;
    float  roughness;
};

//float2 reflect(const float2 i, const float2 n) { return i - 2.0f * n * dot(n,i); }
float3 reflect(const float3 i, const float3 n) { return i - 2.0f * n * dot(n,i); }

float3 get_normal_at_tri_point( const float3 pNew, const float4 p0, const float4 p1, const float4 p2, const float4 vn0, const float4 vn1, const float4 vn2 ) {

	// {	
	//   l0 * p0.x + l1 * p1.x + l2 * p2.x = pNew.x ; 
	//   l0 * p0.y + l1 * p1.y + l2 * p2.y = pNew.y ; 
	//   l0 * p0.z + l1 * p1.z + l2 * p2.z = pNew.z ; 
	// }

	const float Det = p0.x * (p1.y*p2.z - p2.y*p1.z) - p1.x * (p0.y*p2.z - p2.y*p0.z) + p2.x * ( p0.y*p1.z - p1.y*p0.z );

	const float Det_l0 = pNew.x * (p1.y*p2.z - p2.y*p1.z) - p1.x * (pNew.y*p2.z - p2.y*pNew.z) + p2.x * ( pNew.y*p1.z - p1.y*pNew.z );

	const float Det_l1 = p0.x * (pNew.y*p2.z - p2.y*pNew.z) - pNew.x * (p0.y*p2.z - p2.y*p0.z) + p2.x * ( p0.y*pNew.z - pNew.y*p0.z );

	const float Det_l2 = p0.x * (p1.y*pNew.z - pNew.y*p1.z) - p1.x * (p0.y*pNew.z - pNew.y*p0.z) + pNew.x * ( p0.y*p1.z - p1.y*p0.z );

	const float3 l = (float3) ( Det_l0 / Det, Det_l1 / Det, Det_l2 / Det );

	//const float4 normNew = l.x * vn0 + l.y * vn1 + l.z * vn2;
	const float3 normNew = l.x * vn0.xyz + l.y * vn1.xyz + l.z * vn2.xyz;

	//const float3 normNew = l.x * normalize(vn0.xyz) + l.y * normalize(vn1.xyz) + l.z * normalize(vn2.xyz);

	return normNew;//.xyz;

	
}

// intersect ray with a box
// http://www.siggraph.org/education/materials/HyperGraph/raytrace/rtinter3.htm

int intersectBox(float4 r_o, float4 r_d, float4 boxmin, float4 boxmax, float *tnear, float *tfar)
{
    // compute intersection of ray with all six bbox planes
    float4 invR = (float4)(1.0f,1.0f,1.0f,1.0f) / r_d;
    float4 tbot = invR * (boxmin - r_o);
    float4 ttop = invR * (boxmax - r_o);

    // re-order intersections to find smallest and largest on each axis
    float4 tmin = min(ttop, tbot);
    float4 tmax = max(ttop, tbot);

    // find the largest tmin and the smallest tmax
    float largest_tmin = max(max(tmin.x, tmin.y), max(tmin.x, tmin.z));
    float smallest_tmax = min(min(tmax.x, tmax.y), min(tmax.x, tmax.z));

	*tnear = largest_tmin;
	*tfar = smallest_tmax;

	return smallest_tmax > largest_tmin;
}

uint rgbaFloatToInt(float4 rgba)
{
    rgba.x = clamp(rgba.x,0.0f,1.0f);  
    rgba.y = clamp(rgba.y,0.0f,1.0f);  
    rgba.z = clamp(rgba.z,0.0f,1.0f);  
    rgba.w = clamp(rgba.w,0.0f,1.0f);  
    return ((uint)(rgba.w*255.0f)<<24) | ((uint)(rgba.z*255.0f)<<16) | ((uint)(rgba.y*255.0f)<<8) | (uint)(rgba.x*255.0f);
}

__kernel void
d_render(__global uint *d_output, 
         uint imageW, uint imageH,
         float density, float brightness,
         float transferOffset, float transferScale,
         __constant float* invViewMatrix
 #ifdef IMAGE_SUPPORT
          ,__read_only image3d_t volume,
          __read_only image2d_t transferFunc,
          sampler_t volumeSampler,
          sampler_t transferFuncSampler
 #endif
         )

{	
    uint x = get_global_id(0);
    uint y = get_global_id(1);

    float u = (x / (float) imageW)*2.0f-1.0f;
    float v = (y / (float) imageH)*2.0f-1.0f;

    //float tstep = 0.01f;
    float4 boxMin = (float4)(-1.0f, -1.0f, -1.0f,1.0f);
    float4 boxMax = (float4)(1.0f, 1.0f, 1.0f,1.0f);

    // calculate eye ray in world space
    float4 eyeRay_o;
    float4 eyeRay_d;

    eyeRay_o = (float4)(invViewMatrix[3], invViewMatrix[7], invViewMatrix[11], 1.0f);   

    float4 temp = normalize(((float4)(u, v, -2.0f,0.0f)));
    eyeRay_d.x = dot(temp, ((float4)(invViewMatrix[0],invViewMatrix[1],invViewMatrix[2],invViewMatrix[3])));
    eyeRay_d.y = dot(temp, ((float4)(invViewMatrix[4],invViewMatrix[5],invViewMatrix[6],invViewMatrix[7])));
    eyeRay_d.z = dot(temp, ((float4)(invViewMatrix[8],invViewMatrix[9],invViewMatrix[10],invViewMatrix[11])));
    eyeRay_d.w = 0.0f;

    // find intersection with box
	float tnear, tfar;
	int hit = intersectBox(eyeRay_o, eyeRay_d, boxMin, boxMax, &tnear, &tfar);
    if (!hit) {
        if ((x < imageW) && (y < imageH)) {
            // write output color
            uint i =(y * imageW) + x;
            d_output[i] = 0;
        }
        return;
    }
	if (tnear < 0.0f) tnear = 0.0f;     // clamp to near plane

    // march along ray from back to front, accumulating color
    temp = (float4)(0.0f,0.0f,0.0f,0.0f);
    float t = tfar;

    for(uint i=0; i<maxSteps; i++) {		
        float4 pos = eyeRay_o + eyeRay_d*t;
        pos = pos*0.5f+0.5f;    // map position to [0, 1] coordinates

        // read from 3D texture        
#ifdef IMAGE_SUPPORT        
        float4 sample = read_imagef(volume, volumeSampler, pos);
        
        // lookup in transfer function texture
        float2 transfer_pos = (float2)((sample.x-transferOffset)*transferScale, 0.5f);
        float4 col = read_imagef(transferFunc, transferFuncSampler, transfer_pos);
#else
        float4 col = (float4)(pos.x,pos.y,pos.z,.25f);
#endif


        // accumulate result
        float a = col.w*density;
        temp = mix(temp, col, (float4)(a, a, a, a));

        t -= tstep;
        if (t < tnear) break;
    }
    temp *= brightness;

    if ((x < imageW) && (y < imageH)) {
        // write output color
        uint i =(y * imageW) + x;
        d_output[i] = rgbaFloatToInt(temp);
    }
}


/*float fminf1(float a, float b)
{
  return a <= b ? a : b;
}

float fmaxf1(float a, float b)
{
  return a >= b ? a : b;
}*/

// convert floating point rgb color to 8-bit integer
int rgbToInt(float r, float g, float b)
{
	r = clamp(r, 0.0f, 255.0f);
	g = clamp(g, 0.0f, 255.0f);
	b = clamp(b, 0.0f, 255.0f);
	//return (int(r)<<16) | (int(g)<<8) | int(b); // notice switch red and blue to counter the GL_BGRA

	//return (convert_uint(r) << 16) | (convert_uint(g)<<8) | (convert_uint(b));
	return (convert_uint(b) << 16) | (convert_uint(g)<<8) | (convert_uint(r));
}

// Ray structure
struct Ray
{	
	float3 ori;
	float3 dir;
	float3 inv_dir;
};

	void RayInit( struct Ray *r, float3 o, float3 d )
	{
		r->ori = o;
		r->dir = d;
		r->dir = normalize(r->dir);
		r->inv_dir = (float3)(1.0/r->dir.x,1.0/r->dir.y,1.0/r->dir.z);
	} 

struct HitRecord
{
	float t;

	float3 color;
	float3 normal;

	int hit_index; 
};

	void HitRecordInit( struct HitRecord *ht ) {
		ht->t = UINT_MAX;
		ht->hit_index = -1; 
		ht->color = (float3)(0,0,0);
	} ; 

	void HitRecord_resetT( struct HitRecord *ht ) {
		ht->t = UINT_MAX;
		ht->hit_index = -1; 
	} ; 

// intersection code --------------------------------------------------

int RayBoxIntersection(const float3 BBMin, const float3 BBMax, const float3 RayOrg, const float3 RayDirInv, float *tmin, float *tmax)
{
	float l1   = (BBMin.x - RayOrg.x) * RayDirInv.x;
	float l2   = (BBMax.x - RayOrg.x) * RayDirInv.x;
	*tmin = fmin(l1,l2);
	*tmax = fmax(l1,l2);

	l1   = (BBMin.y - RayOrg.y) * RayDirInv.y;
	l2   = (BBMax.y - RayOrg.y) * RayDirInv.y;
	*tmin = fmax(fmin(l1,l2), *tmin);
	*tmax = fmin(fmax(l1,l2), *tmax);

	l1   = (BBMin.z - RayOrg.z) * RayDirInv.z;
	l2   = (BBMax.z - RayOrg.z) * RayDirInv.z;
	*tmin = fmax(fmin(l1,l2), *tmin);
	*tmax = fmin(fmax(l1,l2), *tmax);

	return ((*tmax >= *tmin) && (*tmax >= 0.0f));
}

// the classic ray triangle intersection: http://www.cs.virginia.edu/~gfx/Courses/2003/ImageSynthesis/papers/Acceleration/Fast%20MinimumStorage%20RayTriangle%20Intersection.pdf
float RayTriangleIntersection(const struct Ray r, 
										 const float3 v0, 
										 const float3 edge1, 
										 const float3 edge2)
{  

	float3 tvec = r.ori- v0;  
	float3 pvec = cross(r.dir, edge2);  
	float  det  = dot(edge1, pvec);  

	det = 1.0f / det;  

	float u = dot(tvec, pvec) * det;  

	if (u < 0.0f || u > 1.0f)  
		return -1.0f;  

	float3 qvec = cross(tvec, edge1);  

	float v = dot(r.dir, qvec) * det;  

	if (v < 0.0f || (u + v) > 1.0f)  
		return -1.0f;  

	return dot(edge2, qvec) * det;  
}  

int RaySphereIntersection(const struct Ray ray, const float3 sphere_center, const float sphere_radius, float *t)
{
	float b, c, d;

	float3 sr = ray.ori - sphere_center;
	b =  dot(sr,ray.dir);
	c = dot(sr,sr) - (sphere_radius*sphere_radius);
	d = b*b - c;

	if (d > 0) 
	{
		float e = sqrt(d);
		float t0 = -b-e;
		if(t0 < 0)
			*t = -b+e;
		else
			*t = min(-b-e,-b+e);
		return 1;
	}
	return 0;
}

struct Params {

	float3 a; 
	float3 b; 
	float3 c; 
   float3 campos;
   float3 light_pos;
   float3 light_color;
   float3 scene_aabb_min; 
   float3 scene_aabb_max;
};

typedef struct Params InitParams1;

struct Params2 {

	float4 a; 
	float4 b; 
	float4 c; 
   float4 campos;
   float4 light_pos;
   float4 light_color;
   float4 scene_aabb_min; 
   float4 scene_aabb_max;
};

typedef struct Params2 InitParams2;

//const static int PHONG = 0x001;
//const static int COOK_TORRANCE = 0x002;

#define PHONG 0x001
#define COOK_TORRANCE 0x002

struct Mat {

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
};

typedef struct Mat Material;

__kernel void raytracer( __global uint *out_data,
						   const uint w,
						   const uint h,
						   __global float4 * triangles,
						   const int number_of_triangles,
						   __global float * params
						   //__global struct Params * params
						   //__global int * p0,
						   //__global int *p_res,
						   //__global float * p_res2,
						   //__read_only image2d_t tris
						   )
{

	//const sampler_t  samp = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
	 //const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_LINEAR | CLK_ADDRESS_REPEAT;
	 //float4 lValue = read_imagef(tris, sampler, (int2)(1,0));

	//p_res[0] = 2 ; 

	///*
	int i = 0 ; 

	   const float3 a = (float3) ( params[i], params[++i], params[++i] ); 
	   const float3 b = (float3) ( params[++i], params[++i], params[++i] ); 
	   const float3 c = (float3) ( params[++i], params[++i], params[++i] ); 

	   const float3 campos = (float3) ( params[++i], params[++i], params[++i] ); 

	   const float3 light_pos = (float3) ( params[++i], params[++i], params[++i] ); 
	   const float3 light_color = (float3) ( params[++i], params[++i], params[++i] ); 

	   const float3 scene_aabb_min = (float3) ( params[++i], params[++i], params[++i] ); 
	   const float3 scene_aabb_max = (float3) ( params[++i], params[++i], params[++i] ); 
   //*/
   /*
   	   const float3 a = params->a; 
	   const float3 b = params->b;
	   const float3 c = params->c;

	   const float3 campos = params->campos; 

	   const float3 light_pos = params->light_pos; 
	   const float3 light_color = params->light_color; 

	   const float3 scene_aabb_min = params->scene_aabb_min; 
	   const float3 scene_aabb_max = params->scene_aabb_max; 
	*/

	unsigned int x = get_global_id(0) ; // blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int y = get_global_id(1) ; // blockIdx.y*blockDim.y + threadIdx.y;
	
	if ( y * w + x >= w*h ) return ; 
	
	//out_data[ y * w + x ] = 0x005555; // (y * w + x) &0xFFFFFFFF;
	//return ; 

	//int2 pos = (int2)(get_global_id(0), get_global_id(1));

	float xf = (x-0.5)/((float)w);
	float yf = (y-0.5)/((float)h);

	int ray_depth = 0;
	bool continue_path = true;

	float3 t1 = c+(a*xf);
	float3 t2 = b*yf;
	float3 image_pos = t1 + t2;

	struct Ray r;
	RayInit ( &r, image_pos,image_pos-campos );

	struct HitRecord hit_r;
	HitRecordInit ( &hit_r );

	float t_min,t_max;
	continue_path = RayBoxIntersection(scene_aabb_min, scene_aabb_max, r.ori, r.inv_dir, &t_min, &t_max );

	hit_r.color = (float3)(0,0,0);

	// hack to display the light source we simple make a ray sphere intersection and 
	// compare the depth with the found t value from the triangles
	float sphere_t;
	bool sphere_hit = RaySphereIntersection(r,light_pos,2.0, &sphere_t);

	if(sphere_hit && sphere_t > 0.001)
	{
		if(!continue_path)
		{
			hit_r.color = light_color;
		}
		sphere_hit = true;
	}
	
	//float4 v0 = triangles[10];
	
	for(int i = 0; i < number_of_triangles; i++) {
			float4 v0 = triangles[i*3];
			float4 e1 = triangles[i*3+1];
			float4 e2 = triangles[i*3+2];

			float t = RayTriangleIntersection(r, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

			if(t < hit_r.t && t > 0.001)
			{
				hit_r.t = t; 
				hit_r.hit_index = i;
			}
	}
			
		if(sphere_hit && sphere_t < hit_r.t)
		{
			hit_r.color += light_color;
			continue_path = false;
			//break;
		}
		
		if(hit_r.hit_index >= 0)
		{
			ray_depth++;

			float4 e1 = triangles[ hit_r.hit_index*3+1 ];
			float4 e2 = triangles[ hit_r.hit_index*3+2 ];
			
			// create the normal
			//float4 e1 = tex1Dfetch(triangle_texture,hit_r.hit_index*3+1);
			//float4 e2 = tex1Dfetch(triangle_texture,hit_r.hit_index*3+2);

			hit_r.normal = cross((float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));
			hit_r.normal = normalize(hit_r.normal);

			// calculate simple diffuse light
			float3 hitpoint = r.ori + r.dir *hit_r.t;
			float3 L = light_pos - hitpoint;
			float dist_to_light = length(L);
			
			L = normalize(L);
			float diffuse_light = fmax( dot(L,hit_r.normal), 0.0f);
			diffuse_light = fmin( (diffuse_light),1.0f);
			//calculate simple specular light
			float3 H = L + (-r.dir);
			H = normalize(H);
			float specular_light = pow(fmax(dot(H,hit_r.normal),0.0f), 25.0f);

			diffuse_light  *=  16.0f/dist_to_light;
			specular_light *=  16.0f/dist_to_light;

			clamp(diffuse_light, 0.0f, 1.0f);
			clamp(specular_light, 0.0f, 1.0f);

			hit_r.color += light_color * diffuse_light + (float3)(1.0,1.0,1.0)*specular_light*0.2f + (float3)(0.2,0.2,0.2);
		}
		else
		{
			continue_path = false;
			hit_r.color += (float3)(0.5,0.5,0.95*yf+0.3);
		}

	/*while(continue_path && ray_depth < 1) // 4
	{
		
		// search through the triangles and find the nearest hit point
		for(int i = 0; i < number_of_triangles; i++)
		{
			//float4 v0 = read_imagef( tris, samp, (int2)(i*3,0) ).xyzw;
			//float4 e1 = read_imagef( tris, samp, (int2)(i*3+1,0) ).xyzw;
			//float4 e2 = read_imagef( tris, samp, (int2)(i*3+2,0) ).xyzw;

			float4 v0 = triangles[i*3];
			float4 e1 = triangles[i*3+1];
			float4 e2 = triangles[i*3+2];

			//float4 v0 = tex1Dfetch(triangle_texture,i*3);
			//float4 e1 = tex1Dfetch(triangle_texture,i*3+1);
			//float4 e2 = tex1Dfetch(triangle_texture,i*3+2);

			float t = RayTriangleIntersection(r, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

			if(t < hit_r.t && t > 0.001)
			{
				hit_r.t = t; 
				hit_r.hit_index = i;
			}
		}

		if(sphere_hit && sphere_t < hit_r.t)
		{
			hit_r.color += light_color;
			continue_path = false;
			break;
		}

		if(hit_r.hit_index >= 0)
		{
			ray_depth++;

			float4 e1 = triangles[ hit_r.hit_index*3+1 ];
			float4 e2 = triangles[ hit_r.hit_index*3+2 ];
			
			// create the normal
			//float4 e1 = tex1Dfetch(triangle_texture,hit_r.hit_index*3+1);
			//float4 e2 = tex1Dfetch(triangle_texture,hit_r.hit_index*3+2);

			hit_r.normal = cross((float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));
			hit_r.normal = normalize(hit_r.normal);

			// calculate simple diffuse light
			float3 hitpoint = r.ori + r.dir *hit_r.t;
			float3 L = light_pos - hitpoint;
			float dist_to_light = length(L);
			
			L = normalize(L);
			float diffuse_light = max( dot(L,hit_r.normal), 0.0);
			diffuse_light = min( (diffuse_light),1.0);
			//calculate simple specular light
			float3 H = L + (-r.dir);
			H = normalize(H);
			float specular_light = pow(fmaxf1(dot(H,hit_r.normal),0.0),25.0f);

			diffuse_light  *=  16.0/dist_to_light;
			specular_light *=  16.0/dist_to_light;

			clamp(diffuse_light, 0.0f, 1.0f);
			clamp(specular_light, 0.0f, 1.0f);

			hit_r.color += light_color * diffuse_light + (float3)(1.0,1.0,1.0)*specular_light*0.2 + (float3)(0.2,0.2,0.2);
		}
		else
		{
			continue_path = false;
			hit_r.color += (float3)(0.5,0.5,0.95*yf+0.3);
		}
	}*/

	if(ray_depth >= 1 || sphere_hit)
	{
		ray_depth = max(ray_depth,1);
		hit_r.color /= ray_depth; // normalize the colors
	}
	else
	{
		hit_r.color = (float3)(0.5,0.5,yf+0.3);
	}

	int val = rgbToInt(hit_r.color.x*255,hit_r.color.y*255,hit_r.color.z*255);
	out_data[y * w + x] = val;
}

struct AABB {
	float3 min, max; 
} ; 

//float2 ray_box ( struct Ray * r, struct AABB * box ) {
float2 ray_box ( struct Ray * r, float4 min, float4 max ) {

	float3 t0 = ( (float3)(min.x, min.y, min.z) - r->ori) / r->dir; // * r->inv_dir; // 
    float3 t1 = ( (float3)(max.x, max.y, max.z) - r->ori) / r->dir; // * r->inv_dir; // 

	float3 tmin = fmin(t0,t1);
    float3 tmax = fmax(t0,t1);

    float tmin1 = fmax( fmax( tmin.x, tmin.y ), tmin.z );
    float tmax1 = fmin( fmin( tmax.x, tmax.y ), tmax.z );

	return (float2)( tmin1, tmax1 );
} ; 

struct bvh_node {
	float4 min, max;

	int offset_left;
	int offset_right;

	int offset_tris;
	int num_tris;
} ; 

#define STACK_SIZE  65  // Size of the traversal stack in local memory.

#define EntrypointSentinel 0x76543210   // Bottom-most stack entry, indicating the end of traversal.

#define TMIN 0.001f

//#define SWAP ( t, a, b ) ( t = a; a = b; b = t; )

#define SWAP_TYPE(type, a, b) \
    { \
        type __swap_temp; \
        __swap_temp = (b); \
        (b) = (a); \
        (a) = __swap_temp; \
    }


//void swap1(bool * a, bool  & b) { }//bool  t = a; a = b; b = t; }
//void swap1(int & a, int  & b) { int  t = a; a = b; b = t; }

///*
// return triangle index in tris
int traverse_bvh ( struct Ray * ray, float * tHit, __global float4 * mesh_vertices , __global int * mesh_indices, __global float4 * bvh_nodes, 
			__global int * bvh_tris_indices, bool needClosestHit, int num_tris, int num_bvh_tris, int num_bvh_nodes, __global int * temp ) {//, __global int4 * temp ) {
//int traverse_bvh ( struct Ray * ray, __global struct bvh_node * bvh_nodes ) {

	//temp[0].x = 11;
	//temp[0].y = 12;

	/*
	temp[1].x = as_int(ray->ori.x);
	temp[1].y = as_int(ray->ori.y);
	temp[1].z = as_int(ray->ori.z);
	temp[1].w = as_int(*tHit);

	temp[2].x = as_int(ray->dir.x);
	temp[2].y = as_int(ray->dir.y);
	temp[2].z = as_int(ray->dir.z);

	return -1; 
	*/

	/*
	for(int i = 0; i < num_bvh_nodes; i++) {
		temp[i].x = as_int(bvh_nodes[i*3+2].x);
		temp[i].y = as_int(bvh_nodes[i*3+2].y);

		temp[i].z = as_int(bvh_nodes[i*3+2].z);
		temp[i].w = as_int(bvh_nodes[i*3+2].w);
	}

	return -1;
	*/

	/* works !!!
	int hit_index = -1;
	for(int i = 0; i < num_tris; i++) {
			//float4 v0 = triangles[i*3];
			//float4 e1 = triangles[i*3+1];
			//float4 e2 = triangles[i*3+2];

			float4 v0 = mesh_vertices[ mesh_indices[i*3+0] ];
			float4 e1 = mesh_vertices[ mesh_indices[i*3+1] ] - v0;
			float4 e2 = mesh_vertices[ mesh_indices[i*3+2] ] - v0;

			float t = RayTriangleIntersection( *ray, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

			if(t < *tHit && t > TMIN) // 0.001f
			{
				*tHit = t; 
				hit_index = i*3;	
				if (!needClosestHit) return hit_index;			
			}
	}

	return hit_index;
	*/

	/*
	int hit_index = -1;
	for(int i = 0; i < num_bvh_tris; i++) {
			//float4 v0 = triangles[i*3];
			//float4 e1 = triangles[i*3+1];
			//float4 e2 = triangles[i*3+2];

			int tri1 = bvh_tris_indices [ i ];

			float4 v0 = mesh_vertices[ mesh_indices[tri1] ];
			float4 e1 = mesh_vertices[ mesh_indices[tri1+1] ] - v0;
			float4 e2 = mesh_vertices[ mesh_indices[tri1+2] ] - v0;

			float t = RayTriangleIntersection( *ray, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

			if(t < *tHit && t > 0.001)
			{
				*tHit = t; 
				hit_index = i;
				if (!needClosestHit) return hit_index;
			}
	}
	
	return hit_index;
	*/

	/*
	int hit_index = -1;
	for(int i = 0; i < num_bvh_nodes; i++) {
			//float4 v0 = triangles[i*3];
			//float4 e1 = triangles[i*3+1];
			//float4 e2 = triangles[i*3+2];

		if ( as_int(bvh_nodes[i*3+2].x) < 0 ) {

			int bvh_node_tris_index = as_int(bvh_nodes[i*3+2].z);

			for(int j = 0; j < as_int(bvh_nodes[i*3+2].w); j++) {

				int tri1 = bvh_tris_indices [ bvh_node_tris_index + j];

				float4 v0 = mesh_vertices[ mesh_indices[tri1] ];
				float4 e1 = mesh_vertices[ mesh_indices[tri1+1] ] - v0;
				float4 e2 = mesh_vertices[ mesh_indices[tri1+2] ] - v0;

				float t = RayTriangleIntersection( *ray, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

				if(t < *tHit && t > 0.001)
				{
					*tHit = t; 
					hit_index = tri1; // i;
				}
			}
		}
	}
	
	return hit_index;
	*/

	//char*   stackPtr;               // Current position in traversal stack.
	//int     nodeAddr;               // Non-negative: current internal node, negative: second postponed leaf.

	int traversalStack[STACK_SIZE];
	int stack_count = 1;

	traversalStack[0] = 0; // Bottom-most entry.

	//traversalStack[0] = EntrypointSentinel; // Bottom-most entry.
	//stackPtr = (char*)&traversalStack[0];
	//nodeAddr = 0;   // Start from the root.

	//float   tmin;                   // t-value from which the ray starts. Usually 0.
	//float   hitT;                   // t-value of the closest intersection.

//	bool flag = false ; 
	int tri_index = -1;

	int flag2 = 0;
	//while ( stackPtr >= 0 ) {
	while ( stack_count > 0 ) 
	{	
	//	temp[0].x = flag2+1;
	//	temp[flag2+1].x = stack_count;
		//temp[0].x = num_bvh_nodes;
		//temp[flag2*2+1].x = stack_count;
		//temp[flag2*2+1].x = 15; // flag2*2+1;
		//temp[8*2+1].x = stack_count;

		//return 0; 
	
		//if ( stack_count > 1 ) return -1; //0; //
		//if ( stack_count > 100 ) return 0; // -1; //
		//return 0; // -1; //
		//else {
		//	++stack_count;

			/*if ( as_int(bvh_nodes[0*3+2].x) >= 0 ) {

				int offset_left = as_int(bvh_nodes[0*3+2].x);
				int offset_right = as_int(bvh_nodes[0*3+2].y);

				++stack_count;
			};*/

			//continue;
		//}
		//if ( stack_count > 10 ) return -1; 

		//struct bvh_node & node = *((struct bvh_node *)(&as_bvh_nodes[nodeIndex]));

		int nodeIndex = traversalStack[stack_count - 1];

		//if ( flag2 < 1 ) temp[flag2*2+1].y = nodeIndex;
	//	temp[flag2+1].y = nodeIndex;

		//if ( nodeIndex >= num_bvh_nodes ) return -1; //0; 
		
		//if (false) { 
		//if (true) {
		//if ( bvh_nodes[nodeIndex*3*4+2*4] >= 0 ) {
		//if ( as_int(bvh_nodes[nodeIndex*3*4+2*4]) >= 0 ) {

		int offset_left = as_int(bvh_nodes[nodeIndex*3+2].x);
		if ( offset_left >= 0 ) { // is node
			///*			
			int offset_right = as_int(bvh_nodes[nodeIndex*3+2].y);

			if ( offset_right < 0 ) return -1;

	//		temp[flag2+1].z = offset_left;

			if ( flag2+1 < 0 ) {//if ( flag2+1 > 0xEFFFFFFF ) { //64 ) { // AMD BUG !!!
				//temp[flag2+1].w = 1; //offset_right;
				//temp[0].x = 0; //offset_right;
				temp[0] = 0; //offset_right; 
			}

			//temp[flag2*2+1].z = offset_left;
			//temp[flag2*2+1].w = offset_right;

			//if ( offset_left >= num_bvh_nodes || offset_right >= num_bvh_nodes ) return -1; 
			if ( offset_left  >= num_bvh_nodes ) return -1; // 0
			if ( offset_right >= num_bvh_nodes ) return -1; // 0

			//BVH_Node_ child0 = bvh_cuda.bvh_nodes [ node.offset_left ];
			//BVH_Node_ child1 = bvh_cuda.bvh_nodes [ node.offset_right ];

			//struct bvh_node child0

			float2 tspan0 = ray_box ( ray, bvh_nodes[offset_left*3+0], bvh_nodes[offset_left*3+1] );
			float2 tspan1 = ray_box ( ray, bvh_nodes[offset_right*3+0], bvh_nodes[offset_right*3+1] );

			bool intersect0 = (tspan0.x <= tspan0.y ) && ( tspan0.y >= TMIN) && (tspan0.x <= *tHit);
			bool intersect1 = (tspan1.x <= tspan1.y ) && ( tspan1.y >= TMIN) && (tspan1.x <= *tHit);

			//temp[flag2*2+1+1].x = 500; //as_int(intersect0);
			//temp[flag2*2+1+1].y = 600; // as_int(intersect1);

			//temp[flag2*2+1+1].z = 1.15f; // as_int(tspan0.x);
			//temp[flag2*2+1+1].w = 2.25f; // as_int(tspan1.x);
			//*/

			/*
			if (true) { //if (true) { //if(intersect0 || intersect1) {
			
				for(int i = 0; i < num_tris; i++) {
					//float4 v0 = triangles[i*3];
					//float4 e1 = triangles[i*3+1];
					//float4 e2 = triangles[i*3+2];

					float4 v0 = mesh_vertices[ mesh_indices[i*3] ];
					float4 e1 = mesh_vertices[ mesh_indices[i*3+1] ] - v0;
					float4 e2 = mesh_vertices[ mesh_indices[i*3+2] ] - v0;
					
					float t = RayTriangleIntersection( *ray, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

					if(t < *tHit && t > 0.001)
					{
						*tHit = t; 
						return tri_index = i*3; // return -1; // 
					}
				}
			}
			*/

			//++stack_count;
			///*
			if(intersect0 && intersect1) {

				if( tspan0.x > tspan1.x ) {
					//SWAP_TYPE(float2,tspan0,tspan1); // swap1
					//SWAP_TYPE(int, offset_left,offset_right); //swap1(child0,child1); // swap1
					int t = offset_left; offset_left = offset_right; offset_right = t;
				}

				//nodeIndex = offset_left;
				//traversalStack [ stack_count - 1 ] = offset_right;
				
				traversalStack [ stack_count - 1 ] = offset_right; 
				
				if ( stack_count >= STACK_SIZE ) return -1;
				traversalStack [ stack_count ] = offset_left; 
				++stack_count ;
				

			} else {
			
				if(intersect0) {
					//nodeIndex = offset_left;
					traversalStack [ stack_count - 1 ] = offset_left;
				} else if (intersect1)  {
					//nodeIndex = offset_right;
					traversalStack [ stack_count - 1 ] = offset_right;
				} else {
					--stack_count;
				};	
				
			}
			//*/

			//return 0; //-1; //

		} else { // is leaf

		//	temp[flag2+1].z = -1;
		//	temp[flag2+1].w = -1;

			//return 0; //-1; //
			/*
			for(int i = 0; i < num_tris; i++) {
					//float4 v0 = triangles[i*3];
					//float4 e1 = triangles[i*3+1];
					//float4 e2 = triangles[i*3+2];

					float4 v0 = mesh_vertices[ mesh_indices[i*3] ];
					float4 e1 = mesh_vertices[ mesh_indices[i*3+1] ] - v0;
					float4 e2 = mesh_vertices[ mesh_indices[i*3+2] ] - v0;

					float t = RayTriangleIntersection( *ray, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

					if(t < *tHit && t > 0.001)
					{
						*tHit = t; 
						//return 
						tri_index = i*3;
						//return 0; 
					}
			}
			*/

			///*
			int node_offset_tris = as_int(bvh_nodes[nodeIndex*3+2].z);
			int node_num_tris = as_int(bvh_nodes[nodeIndex*3+2].w);

			for ( int i=0; i<node_num_tris; ++i ) {

				int tri1 = bvh_tris_indices[ node_offset_tris + i ]; // pNode->tris[i];
				
				float4 v0 = mesh_vertices [ mesh_indices[ tri1 + 0 ] ];
				float4 e1 = mesh_vertices [ mesh_indices[ tri1 + 1 ] ] - v0;
				float4 e2 = mesh_vertices [ mesh_indices[ tri1 + 2 ] ] - v0;

				float t = RayTriangleIntersection(*ray, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

				if(t < *tHit && t > TMIN) {
			
					*tHit = t; 

					//tri.v0 = v0;
					//tri.v1 = e1;
					//tri.v2 = e2;

					if(!needClosestHit) return tri1;
					//return -1;

					//flag = true ; 
					//return 
					tri_index = tri1;
			//		temp[flag2+1].w = tri_index;
				} 
			} // bvh_node_tris
			//*/

			--stack_count;

		} // is leaf

		flag2++;

	} // while ( stack_count > 0 )

	//if ( !flag )	return -1;
	return tri_index;

	//return 0;
	//return -1;
}
//*/

float4 cook_torrance (
			float3 LightDirection1,
			float3 specularColor,
			float3 diffuseColor,
			float roughnessValue, // 0 : smooth, 1: rough
			float F0, // fresnel reflectance at normal incidence
			float k, // fraction of diffuse reflection (specular reflection = 1 - k)
			float3 vnormal,
			float3 vpos
		);

float4 cook_torrance_steps3d (
			float3 l,
			float3 h,
			float3 v,
			float3 n,
			
			float r0, 
			float roughness,

			//uniform sampler2D	lookupMap;

			float4 diffColor1, float4 specColor1
		);

float4 phong_steps3d ( float3 l, float3 v, float3 n,
						float4 diffColor, float4 specColor );

float3 CookTorrance_GGX(float3 n, float3 l, float3 v, struct Material_pbr m);

__kernel void raytracer_bvh( __global __write_only uint *out_data,
						   const uint w,
						   const uint h,
						   __global float4 * triangles,
						   const int number_of_triangles,
						   //__global float * params,
						   //__global InitParams1 * params1,
						   __global InitParams2 * params2,

						   __global float4 * mesh_vertices , 
						   __global int * mesh_indices, 

						   __global float4 * bvh_nodes, 
						   __global int * bvh_tris_indices, 
						   
						   const int num_bvh_tris, 
						   const int num_bvh_nodes,

						   __global int * temp,  //__global int4 * temp 

						   __global float4 * mesh_normals , 
						   __global int * mesh_normals_indices,

						   __global Material * mesh_materials,
						   __global int * mesh_triangle_index_to_material_index
						   
						   
						   //__global struct Params * params,
						   //__global int * p0,
						   //__global int *p_res,
						   //__global float * p_res2,
						   //__read_only image2d_t tris
						   )
{
	//const sampler_t  samp = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
	 //const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE | CLK_FILTER_LINEAR | CLK_ADDRESS_REPEAT;
	 //float4 lValue = read_imagef(tris, sampler, (int2)(1,0));

	//p_res[0] = 2 ; 

	/*
	//__global float * params = (__global float *)params1;

	//__global float3 * params3 = (__global float3 *)params;
	//__global float3 * paramsTemp = (__global float3 *)(params+3);
	//__global float4 * params4 = (__global float4 *)(params+0);

	//__global InitParams1 * paramsIPar = (__global InitParams1 *)params;
	//__global InitParams2 * paramsIPar2 = (__global InitParams2 *)params;
	int i = 0 ; 

	   //const float3 a = (float3) ( ((__global float3 *)params)[0].x, params[++i], params[++i] ); 
	   //const float3 a = (float3) ( params3[0].x, params3[0].y, params3[0].z ); 
	   //const float3 a = params3[0];
	   //const float3 a = paramsIPar[0].a;
	   //const float3 a = paramsIPar->a;
	   //const float3 a = (float3) ( paramsIPar2->a.x, paramsIPar2->a.y, paramsIPar2->a.z ); 
	   //++i; ++i;

	   const float3 a = (float3) ( params[i], params[++i], params[++i] ); // [i]
	   const float3 b = (float3) ( params[++i], params[++i], params[++i] ); 
	   //const float3 b = params3[1];
	   //const float3 b = paramsTemp[0];
	   //const float3 b = (float3) ( params4[0].w, params4[1].x, params4[1].y ); 
	   //const float3 b = (float3) ( params3[1].x, params3[1].y, params3[1].z ); 
	   //const float3 b = paramsIPar[0].b;
	   //const float3 b = paramsIPar->b;
	   //const float3 b = (float3) ( paramsIPar2->a.w, paramsIPar2->b.x, paramsIPar2->b.y ); 
	   //++i; ++i; ++i;

	   const float3 c = (float3) ( params[++i], params[++i], params[++i] ); 
	   //const float3 c = (float3) ( params4[1].z, params4[1].w, params4[2].x ); 
	   //++i; ++i; ++i;

	   const float3 campos = (float3) ( params[++i], params[++i], params[++i] ); 

	   const float3 light_pos = (float3) ( params[++i], params[++i], params[++i] ); 
	   const float3 light_color = (float3) ( params[++i], params[++i], params[++i] ); 

	   const float3 scene_aabb_min = (float3) ( params[++i], params[++i], params[++i] ); 
	   const float3 scene_aabb_max = (float3) ( params[++i], params[++i], params[++i] ); 
   */

   /*
	__global float4 * params = (__global float4 *)params2;
	int i = 0 ; 

	   const float3 a = params[i].xyz;
	   const float3 b = params[++i].xyz;
	   const float3 c = params[++i].xyz;

	   const float3 campos = params[++i].xyz;

	   const float3 light_pos = params[++i].xyz;
	   const float3 light_color = params[++i].xyz;

	   const float3 scene_aabb_min = params[++i].xyz;
	   const float3 scene_aabb_max = params[++i].xyz;
   */
   ///*
   	   const float3 a = params2->a.xyz; 
	   const float3 b = params2->b.xyz; 
	   const float3 c = params2->c.xyz; 

	   const float3 campos = params2->campos.xyz; 

	   const float3 light_pos = params2->light_pos.xyz; 
	   const float3 light_color = params2->light_color.xyz; 

	   const float3 scene_aabb_min = params2->scene_aabb_min.xyz; 
	   const float3 scene_aabb_max = params2->scene_aabb_max.xyz; 
	//*/

	unsigned int x = get_global_id(0) ; // blockIdx.x*blockDim.x + threadIdx.x;
	unsigned int y = get_global_id(1) ; // blockIdx.y*blockDim.y + threadIdx.y;

	//unsigned int x = 400;
	//unsigned int y = 400;
	
	if ( y * w + x >= w*h ) return ; 
	
	//out_data[ y * w + x ] = 0x005555; // (y * w + x) &0xFFFFFFFF;
	//return ; 

	//int2 pos = (int2)(get_global_id(0), get_global_id(1));

	float xf = (x-0.5)/((float)w);
	float yf = (y-0.5)/((float)h);

	int ray_depth = 0;
	bool continue_path = true;

	float3 t1 = c+(a*xf);
	float3 t2 = b*yf;
	float3 image_pos = t1 + t2;

	/*
	float aspect_ratio = (float)(w) / (float)(h);

	float xf = ( (x / (float)w) - 0.5f ) * aspect_ratio;
	float yf = ( (y / (float)h) - 0.5f );
	float3 image_pos = xf * a + //camera_right +
                        yf * b + // camera_up +
                        campos + c; // camera_direction;
	*/

	struct Ray r;
	RayInit ( &r, image_pos, image_pos-campos );

	struct HitRecord hit_r;
	HitRecordInit ( &hit_r );

	float t_min,t_max;
	continue_path = RayBoxIntersection(scene_aabb_min, scene_aabb_max, r.ori, r.inv_dir, &t_min, &t_max );

	hit_r.color = (float3)(0,0,0);
	//hit_r.color = (float3)(1.0f, 1.0f, 1.0f);

	/*
	// hack to display the light source we simple make a ray sphere intersection and 
	// compare the depth with the found t value from the triangles
	float sphere_t;
	bool sphere_hit = RaySphereIntersection(r,light_pos,2.0, &sphere_t);

	if(sphere_hit && sphere_t > 0.001)
	{
		if(!continue_path)
		{
			hit_r.color = light_color;
		} else {
			//hit_r.color = (float3)(0,0,0);
		}
		sphere_hit = true;
	}
	*/
	
	//float4 v0 = triangles[10];
	
	//float t = 10000000; // 
	//int tri_index = 
	//if ( tri_index >= 0 )	
	//	hit_r.hit_index = tri_index;

	float shadow_coef = 1.0f;

	float shadow_coef_2 = 1.0f;

	float shadow_coef_2_resault = 0.0f;

	//float3 vPrev;

	//float3 rez_colors[4]; // rez_colors[RAY_TRACE_DEPTH]

	while( continue_path && ray_depth < RAY_TRACE_DEPTH ) // 4 // 10 // 2
	{
		hit_r.hit_index = traverse_bvh ( &r, &hit_r.t, mesh_vertices , mesh_indices, bvh_nodes, bvh_tris_indices, true, number_of_triangles, num_bvh_tris, num_bvh_nodes, temp ); // tHit

		/*
		for(int i = 0; i < number_of_triangles; i++) {
				//float4 v0 = triangles[i*3];
				//float4 e1 = triangles[i*3+1];
				//float4 e2 = triangles[i*3+2];

				float4 v0 = mesh_vertices[ mesh_indices[i*3] ];
				float4 e1 = mesh_vertices[ mesh_indices[i*3+1] ] - v0;
				float4 e2 = mesh_vertices[ mesh_indices[i*3+2] ] - v0;

				float t = RayTriangleIntersection(r, (float3)(v0.x,v0.y,v0.z), (float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));

				if(t < hit_r.t && t > 0.001)
				{
					hit_r.t = t; 
					hit_r.hit_index = i;
				}
		}
		*/

		shadow_coef_2 = 1.0f;
			
		/*
		sphere_hit = RaySphereIntersection(r,light_pos,2.0, &sphere_t);

		if(sphere_hit && sphere_t > 0.001 && sphere_t < hit_r.t)
		{
		//	hit_r.color += light_color * shadow_coef;

			//hit_r.color *= light_color * shadow_coef; //_2 ?

			if (0 == ray_depth) hit_r.color = (float3)(1.0f, 1.0f, 1.0f);

			//hit_r.color *= light_color;

			continue_path = false;
			//ray_depth++;

			//shadow_coef_2_resault += shadow_coef_2;

			break;
		}
		*/
		
		if(hit_r.hit_index >= 0 )
		{
			ray_depth++;

			/*
			//float4 e1 = triangles[ hit_r.hit_index*3+1 ];
			//float4 e2 = triangles[ hit_r.hit_index*3+2 ];

			// *3
			float4 v0 = mesh_vertices [ mesh_indices[ hit_r.hit_index + 0 ] ];
			float4 e1 = mesh_vertices [ mesh_indices[ hit_r.hit_index + 1 ] ] - v0;
			float4 e2 = mesh_vertices [ mesh_indices[ hit_r.hit_index + 2 ] ] - v0;
			
			// create the normal
			//float4 e1 = tex1Dfetch(triangle_texture,hit_r.hit_index*3+1);
			//float4 e2 = tex1Dfetch(triangle_texture,hit_r.hit_index*3+2);

			hit_r.normal = cross((float3)(e1.x,e1.y,e1.z), (float3)(e2.x,e2.y,e2.z));
			hit_r.normal = normalize(hit_r.normal);
			*/

			///*
			const float4 v0 = mesh_vertices [ mesh_indices[ hit_r.hit_index + 0 ] ];
			const float4 v1 = mesh_vertices [ mesh_indices[ hit_r.hit_index + 1 ] ];
			const float4 v2 = mesh_vertices [ mesh_indices[ hit_r.hit_index + 2 ] ];

			const float4 vn0 = mesh_normals [ mesh_normals_indices[ hit_r.hit_index + 0 ] ];
			const float4 vn1 = mesh_normals [ mesh_normals_indices[ hit_r.hit_index + 1 ] ];
			const float4 vn2 = mesh_normals [ mesh_normals_indices[ hit_r.hit_index + 2 ] ];

			const float3 vNew = r.ori + r.dir * ( hit_r.t - 0.001f );
			
			hit_r.normal = get_normal_at_tri_point ( vNew, v0, v1, v2, vn0, vn1, vn2 ); // (vn0.xyz + vn1.xyz + vn2.xyz)/3; //(float3)(1,1,1); //
			hit_r.normal = normalize( hit_r.normal );
			//*/

			/*
			float3 viewPos = r.ori;
			float3 n = hit_r.normal;
			float3 l = normalize ( light_pos - viewPos );
			float3 e = normalize ( -r.dir ); 			
			//float3 e = ray_depth <= 1 ? normalize ( -r.dir ) : normalize ( vNew - vPrev ); 			
			
			float3 r1 = normalize(-reflect(l, n));

			float3 Idiff = light_color * fmax(dot(n, l), 0.0f);
			Idiff = clamp(Idiff, 0.0f, 1.0f);

			float3 Ispec = pow(fmax(dot(r1, e), 0.0f), 10.0f ); 
			Ispec = clamp(Ispec, 0.0f, 1.0f);
 
			float3 rez_color = Idiff + Ispec*0.2;
			*/

			/*
			cook_torrance(	float3 LightDirection1,
							float3 specularColor,
							float3 diffuseColor,
							float roughnessValue, // 0 : smooth, 1: rough
							float F0, // fresnel reflectance at normal incidence
							float k, // fraction of diffuse reflection (specular reflection = 1 - k)
							float3 vnormal,
							float3 vpos )
			*/

			//float3 rez_color = cook_torrance( (float3)(0,-1,-1), (float3)(1.0f, 1.0f, 1.0f), (float3)(1.0f, 1.0f, 1.0f), 0.1f, 0.8f, 0.7f, hit_r.normal, vNew ).xyz;

			//vec3	p = vec3 ( gl_ModelViewMatrix * gl_Vertex );			// transformed point to world space
			//l = normalize ( vec3 ( lightPos ) - p );					// vector to light source
			//v = normalize ( vec3 ( eyePos )   - p );					// vector to the eye
			//h = normalize ( l + v );
			//n = normalize ( gl_NormalMatrix * gl_Normal );							// transformed n

			float3 l1 = normalize ( light_pos - vNew ); // vector to light source
			float3 v = normalize ( r.ori - vNew ); // vector to the "eye"" 
			float3 h = normalize ( l1 + v );
			float3 n = normalize ( hit_r.normal ); // normal to point on polygon

			/*
			float4 cook_torrance_step3d (
											float3 l,
											float3 h,
											float3 v,
											float3 n,
			
											float r0, 
											float roughness
										)
			*/

			Material mat = mesh_materials[ mesh_triangle_index_to_material_index[ hit_r.hit_index/3 ] ];

			const float diff_koef = 1.25f;

			//float3 rez_color = cook_torrance_steps3d( l1, h, v, n, 0.5, 0.15, mat.diffuse * diff_koef, mat.specular ).xyz;
			//float3 rez_color = phong_steps3d(l1, v, n, mat.diffuse * diff_koef, mat.specular).xyz;

			struct Material_pbr m;
			m.albedo = mat.diffuse.xyz; // (float3)( 0.01f, 0.01f, 0.01f ); // 
			//m.f0 = (float3)( 255, 219, 145 ) * (1/255.0f);
			m.f0 = (float3)( 40, 40, 40 ) * (1/255.0f);
			//m.f0 = (float3)( 127, 127, 127 ) * (1/255.0f);
			//m.f0 =  mat.diffuse.xyz;
			m.roughness = 0.5f;

			float3 rez_color = CookTorrance_GGX(n, l1, v, m) * 3.0f; // 1.0f; // 
			
			/*
			float3 rez_color = (float3)(1.0f, 1.0f, 1.0f);
			
			if ( PHONG == mat.technique.x ) {
				rez_color = phong_steps3d(l1, v, n, mat.diffuse * diff_koef, mat.specular).xyz;
			} else if ( COOK_TORRANCE == mat.technique.x ) { 
				rez_color = cook_torrance_steps3d( l1, h, v, n, 0.5, 0.15, mat.diffuse * diff_koef, mat.specular ).xyz;
			} else {
			}
			*/

			float3 ambient_color = (float3)(0.3f,0.3f,0.3f) * mat.diffuse.xyz; // (float3)(0.18f,0.18f,0.18f) // (float3)(0.1f,0.1f,0.1f); // (float3)(0.0f,0.0f,0.0f); // 
			rez_color += ambient_color;
			
			// calculate simple diffuse light
			//float3 hitpoint = r.ori + r.dir *hit_r.t;
			float3 hitpoint = vNew;

			
			float3 L = light_pos - hitpoint;
			//float dist_to_light = length(L);
			
			L = normalize(L);
			/*
			float diffuse_light = fmax( dot(L,hit_r.normal), 0.0f);
			diffuse_light = fmin( (diffuse_light),1.0f);

			//calculate simple specular light
			float3 H = L + (-r.dir);
			H = normalize(H);
			float specular_light = pow(fmax(dot(H,hit_r.normal),0.0f), 25.0f);

			diffuse_light  *=  16.0f/dist_to_light;
			specular_light *=  16.0f/dist_to_light;

			clamp(diffuse_light, 0.0f, 1.0f);
			clamp(specular_light, 0.0f, 1.0f);

			//hit_r.color += light_color * diffuse_light + (float3)(1.0,1.0,1.0)*specular_light*0.2 + (float3)(0.2,0.2,0.2);
			float3 rez_color = light_color * diffuse_light + (float3)(1.0,1.0,1.0)*specular_light*0.2 ;//+ (float3)(0.2,0.2,0.2);
			*/
			
			//shadow_coef_2 = 1.0f;

			///*
			//if ( ray_depth < 2 ) 
			{
				struct Ray ray_shadow;

				// create a shadow ray
				RayInit ( &ray_shadow, hitpoint + L * 0.001f, L ); //Ray shadow_ray(hitpoint, L);
				HitRecord_resetT( &hit_r ); //hit_r.resetT();
				hit_r.hit_index = traverse_bvh ( &ray_shadow, &hit_r.t, mesh_vertices , mesh_indices, bvh_nodes, bvh_tris_indices, false, number_of_triangles, num_bvh_tris, num_bvh_nodes, temp ); // tHit
				if ( hit_r.hit_index >= 0 && hit_r.t > 0.025f) {
					//hit_r.color *= 0.25f;
					//rez_color *= 0.25f;
//					shadow_coef *= 0.25f;

					shadow_coef_2 = 0.25f; // 0.55f; // 0.95f; //

					//rez_color = ambient_color ; // (float3)(1.0f, 1.0f, 1.0f); // (float3)(0.5f, 0.5f, 0.5f); // 

					//hit_r.color += (mat.diffuse.xyz * diff_koef + (float3)(0.2,0.2,0.2)) * shadow_coef_2; //* refl_koef; 

					//break;
				} else {

					//hit_r.color += (rez_color + (float3)(0.2,0.2,0.2)) * shadow_coef_2; //* refl_koef; 
				}
			}
			//*/

			//float refl_koef = 1.0f;
			//if ( ray_depth > 1 ) refl_koef = 0.9f;

			//hit_r.color += rez_color * shadow_coef; //+ (float3)(0.2,0.2,0.2); // (float3)(0.1,0.1,0.1); // 
		//	hit_r.color += (rez_color + (float3)(0.2,0.2,0.2)) * shadow_coef; //* refl_koef; 
			//hit_r.color += (rez_color + (float3)(0.1,0.1,0.1)) * shadow_coef; 

	//		hit_r.color += (rez_color + (float3)(0.2,0.2,0.2)) * shadow_coef_2; //* refl_koef; 

			//if (1 == ray_depth) {
				//hit_r.color = (float3)(1.0f, 1.0f, 1.0f);
			//}

			//hit_r.color *= (rez_color + (float3)(0.2,0.2,0.2)); //* shadow_coef_2;

			//hit_r.color *= rez_color + (float3)(0.2,0.2,0.2) * ray_depth;// * shadow_coef_2;
			//hit_r.color *= rez_color ; //* ray_depth ; //* shadow_coef_2;

			//hit_r.color *= rez_color ;//* (ray_depth+0)*1.5f ; //* shadow_coef_2;

			hit_r.color += rez_color;

			//rez_colors[ray_depth-1] = rez_color;

			shadow_coef_2_resault += shadow_coef_2;

			//hit_r.color += rez_color * shadow_coef_2; //
			//hit_r.color += (rez_color + (float3)(0.1,0.1,0.1)) * shadow_coef_2; 
				
			///*		
			{ // reflect
				//vPrev = vNew;
				float3 refl = reflect( r.dir, hit_r.normal ); //normalize( reflect(r.dir,hit_r.normal) );
				RayInit ( &r, hitpoint + refl * 0.001f, refl ); //r = Ray(hitpoint, reflect(r.dir,hit_r.normal));
				HitRecord_resetT( &hit_r ); //hit_r.resetT();
			}
			//*/
		}
		else
		{
			continue_path = false;
			//hit_r.color += (float3)(0.5,0.5,0.95*yf+0.3) * shadow_coef;
//			hit_r.color += (float3)(0.5,0.5,yf+0.3) * shadow_coef;

			//ray_depth++;

			//hit_r.color *= (float3)(0.9f, 0.9f, 0.9f);

			//shadow_coef_2 = 0.0f; //0.25f;
		}

		//shadow_coef_2_resault += shadow_coef_2;

		//ray_depth++;
	}

	if(ray_depth >= 1 ) //|| sphere_hit)
	{
		//ray_depth = max(ray_depth,1);
	//	hit_r.color /= ray_depth; // normalize the colors

		/*
		hit_r.color = (float3)(1.0f,1.0f,1.0f);

		for(int i=0; i<ray_depth ; ++i){
			hit_r.color *= rez_colors[i];
		}
		*/

		//if ( ray_depth >= 2 )
		//	hit_r.color *= pow(2.0f, (ray_depth-1) ) ; // 2*(ray_depth-1); // 

		hit_r.color /= ray_depth ; // 2*(ray_depth-1);

		shadow_coef_2_resault /= ray_depth;
		hit_r.color *= shadow_coef_2_resault;
	}
	else
	{
		hit_r.color = (float3)(0.0f,0.0f,0.0f); //(float3)(0.5,0.5,yf+0.3);
	}

	int val = rgbToInt(hit_r.color.x*255,hit_r.color.y*255,hit_r.color.z*255);
	out_data[y * w + x] = val;
}

///*
// NOT IN USE - SEE function "cook_torrance_steps3d"
float4 cook_torrance (
			float3 LightDirection1,
			float3 specularColor,
			float3 diffuseColor,
			float roughnessValue, // 0 : smooth, 1: rough
			float F0, // fresnel reflectance at normal incidence
			float k, // fraction of diffuse reflection (specular reflection = 1 - k)
			float3 vnormal,
			float3 vpos
		) 
{
	float3 lightDirection = normalize(-LightDirection1); // to light
	float3 normal = normalize(vnormal);

	float NdotL = dot(normal, lightDirection);

	float specular = 0.0;

	if(NdotL > 0.0) {
	
		float3 eyeDir = normalize(-vpos); // to eye
		float3 halfVector = normalize(lightDirection + eyeDir);
		float NdotH = max(dot(normal, halfVector), 0.0f); 
		float NdotV = max(dot(normal, eyeDir), 0.0f);
		float VdotH = max(dot(eyeDir, halfVector), 0.0f);
		float mSquared = roughnessValue * roughnessValue;
			        
		// geometric attenuation
		// Blinn's model
		float NH2 = 2.0f * NdotH;
		float g1 = (NH2 * NdotV) / VdotH;
		float g2 = (NH2 * NdotL) / VdotH;
		float geoAtt = min(1.0f, min(g1, g2));
			     
		// roughness (microfacet distribution function)
		// Beckmann distribution
		float r1 = 1.0f / ( M_PI_F * mSquared * pow(NdotH, 4.0f));

		float r2 = (NdotH * NdotH - 1.0f) / (mSquared * NdotH * NdotH);
		float roughness = r1 * exp(r2);
			        
		// Fresnel
		// Schlick's approximation
		float fresnel = pow(1.0f - VdotH, 5.0f);
		fresnel *= (1.0f - F0);
		fresnel += F0;
			        
		specular = (fresnel * geoAtt * roughness) / (NdotV * NdotL * M_PI_F);
	}

	float3 finalValue = NdotL * ((k * diffuseColor) + (specularColor * specular * (1.0f - k)));
	float4 gl_FragColor = (float4)(finalValue, 1.0f);

	return gl_FragColor;
}
//*/

// NOT IN USE - SEE function "phong_steps3d"
float4 phong(
			float3 vnormal,
			float3 vpos,
			//varying mat4 mvMatrix,
			float3 LightDirection1,
			float3 specularColor,
			float3 diffuseColor,
			float3 ambientColor,
			float Smoothness
		) {

	float3 n = normalize(vnormal);
	float4 diffuse = (float4)(0.0f);
	float4 specular = (float4)(0.0f);
				
	// ambient
	float4 ambient = (float4)(ambientColor, 1.0f);
	// diffuse
	float4 kd = (float4)(diffuseColor, 1.0f);
	float3 lightDir = normalize(-LightDirection1);
	float NdotL = dot(n, lightDir);
				
	if (NdotL > 0.0f)
		diffuse = kd * NdotL;
				
	// specular
	float4 ks = (float4)(specularColor, 1.0f);
	float3 rVector = normalize(2.0f * n * NdotL - lightDir);
	float3 viewVector = normalize(-vpos);
	float RdotV = dot(rVector, viewVector);
				
	if (RdotV > 0.0f)
		specular = ks * pow(RdotV, Smoothness);
	float4 gl_FragColor = ambient + diffuse + specular;

	return gl_FragColor;
}

//float	fresnel ( float ca ) {
//	return (r0 + (1.0 - r0)*pow ( 1.0 - ca, 5.0)) / ca;
//}

// NOT IN USE - SEE function "CookTorrance_GGX"
float4 cook_torrance_steps3d (
			float3 l,
			float3 h,
			float3 v,
			float3 n,
			
			float r0, 
			float roughness,

			//uniform sampler2D	lookupMap

			float4 diffColor1, float4 specColor1
		)
		
{
	//const float4	diffColor = (float4) ( 0.5, 0.0, 0.0, 1.0 );
	const float4	specColor = (float4) ( 0.7, 0.7, 0.0, 1.0 );

	const float4	diffColor = (float4) ( 0.5, 0.5, 0.5, 1.0 );
	//const float4	specColor = (float4) ( 0.7, 0.7, 0.7, 1.0 );
	const float	e         = 2.7182818284;
	//const float pi        = 3.1415926;

	float3	n2   = normalize ( n );
	float3	l2   = normalize ( l );
	float3	v2   = normalize ( v );
	float3	h2   = normalize ( h );
	float	nh   = dot ( n2, h2 );
	float	nv   = dot ( n2, v2 );
	float	nl   = dot ( n2, l2 );
	//float	d    = texture2D ( lookupMap, vec2 ( roughness, nh ) ).x;
///*	
	float	r2   = roughness * roughness;
	float	nh2  = nh * nh;
	float   ex   = - (1.0f - nh2)/(nh2 * r2);
	float	d    = pow ( e, ex ) / (r2*nh2*nh2);
//*/	
	float	f    = mix ( pow ( 1.0f - nv, 5.0f ), 1.0f, r0 );		// Fresnel
	float	x    = 2.0f * nh / dot ( v2, h2 );
	float	g    = min ( 1.0f, min ( x * nl, x * nv ) );			// Geometry attenuation
	float	ct   = d*f*g / nv;
	
	float4	diff = diffColor1 * max ( 0.0f, nl );
	float4	spec = specColor * max ( 0.0f, ct );

	float4 gl_FragColor = diff + spec; //

	return gl_FragColor;
}

// NOT IN USE 
float4 phong_steps3d ( float3 l, float3 v, float3 n,
						float4 diffColor1, float4 specColor1)
{
	//const float4	diffColor = (float4) ( 0.5, 0.0, 0.0, 1.0 );
	const float4	specColor = (float4) ( 0.7, 0.7, 0.0, 1.0 );

	//const float4	diffColor = (float4) ( 0.5, 0.5, 0.5, 1.0 );
	//const float4	specColor = (float4) ( 0.7, 0.7, 0.7, 1.0 );

	//const float4	diffColor = (float4) ( 0.584314,  0.686275, 0.172549, 1.0 );
	//const float4	specColor = (float4) ( 0.584314,  0.686275, 0.172549, 1.0 );

	
	

	const float	specPower = 30.0f;

	float3	n2   = normalize ( n );
	float3	l2   = normalize ( l );
	float3	v2   = normalize ( v );
	float3	r    = reflect ( -v2, n2 );
	float4	diff = diffColor1 * max ( dot ( n2, l2 ), 0.0f );
	float4	spec = specColor * pow ( max ( dot ( l2, r ), 0.0f ), specPower );

	float4 gl_FragColor = diff + spec;

	return gl_FragColor;
}

float GGX_PartialGeometry(float cosThetaN, float alpha)
{
    float cosTheta_sqr = clamp( cosThetaN*cosThetaN, 0.0f, 1.0f ); // saturate(cosThetaN*cosThetaN);
    float tan2 = ( 1 - cosTheta_sqr ) / cosTheta_sqr;
    float GP = 2 / ( 1 + sqrt( 1 + alpha * alpha * tan2 ) );
    return GP;
}

float GGX_Distribution(float cosThetaNH, float alpha)
{
    float alpha2 = alpha * alpha;
    float NH_sqr = clamp( cosThetaNH * cosThetaNH, 0.0f, 1.0f ); // saturate(cosThetaNH * cosThetaNH);
    float den = NH_sqr * alpha2 + (1.0f - NH_sqr);
    return alpha2 / ( M_PI_F * den * den );
}

float3 FresnelSchlick(float3 F0, float cosTheta) {
    //return F0 + (1.0f - F0) * pow(1.0f - saturate(cosTheta), 5.0f);
	return F0 + (1.0f - F0) * pow(1.0f - clamp(cosTheta, 0.0f, 1.0f ), 5.0f);
}

// IN USE
float3 CookTorrance_GGX(float3 n, float3 l, float3 v, struct Material_pbr m) {
    n = normalize(n);
    v = normalize(v);
    l = normalize(l);
    float3 h = normalize(v+l);
    //precompute dots
    float NL = dot(n, l);
    if (NL <= 0.0f) return 0.0f;
    float NV = dot(n, v);
    if (NV <= 0.0f) return 0.0f;
    float NH = dot(n, h);
    float HV = dot(h, v);
    
    //precompute roughness square
    float roug_sqr = m.roughness*m.roughness;
    
    //calc coefficients
    float G = GGX_PartialGeometry(NV, roug_sqr) * GGX_PartialGeometry(NL, roug_sqr);
    float D = GGX_Distribution(NH, roug_sqr);
    float3 F = FresnelSchlick(m.f0, HV);

    //mix
    float3 specK = G*D*F*0.25f/(NV+0.001f);    
    float3 diffK = clamp(1.0f-F, 0.0f, 1.0f ); //saturate(1.0-F);
    return max(0.0, m.albedo*diffK*NL/M_PI_F + specK);
}