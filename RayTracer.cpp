#include "RayTracer.h"

#include <GL/glut.h>
#include <cmath>
#include <malloc.h>

// Includes
#include <iostream>
#include <string>

#include <fstream>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include "common.h"
//#include "vectors_math.h"

#include "BVH.h"
#include "BVH2.h"

//#include "Util.h"
//#include "Platform.h"
//#include "BVHNode.h"
//#include "SplitBVHBuilder.h"

#include "BVH_Cuda.h"

#include "vectors_math.cpp"

//#include "AABB.h"

#include "Camera.h"

#include "ColladaLoader.h"

#define REFRESH_DELAY	  10 //ms

#define WIDTH  1024 // 1280 // 640 // 
#define HEIGHT 768  // 1080 // 720  // 480 // 

Mesh mesh1;
BVH bvh;

FW::BVH2 bvh2;

BVH_Cuda bvh_cuda;

int winwidth = WIDTH;
int winheight = HEIGHT;

GLuint texnum;

int frames = 0;
int t0 = 0, te;

float light_x = -23;
float light_y = 200; // 25; // 100; // 
float light_z = 3;
float3 light_pos ( light_x, light_y, light_z ) ; 

float light_color[3] = {1,1,1};
float3 light_color1 ( light_color[0], light_color[1], light_color[2] ) ; 

float viewRotation[3];
float viewTranslation[3] = {0.0, 0.0, -4.0f};
float invViewMatrix[12];

float delta_t = 0;

// Camera parameters -----------------------------
float3 a; float3 b; float3 c; 
float3 campos; 
float camera_rotation = 0;
float camera_distance = 400; //75; // 25; //
float camera_height = 150; // 25; // 10; // 
bool animate = false; // true; // 

Camera cam;

const int w = WIDTH;  // 512;
const int h = HEIGHT; // 512;
const unsigned int image_width   = w;
const unsigned int image_height  = h;
int out_data[ w * h ] ;

int total_number_of_triangles = 0;

// Scene bounding box ----------------------------
float3 scene_aabbox_min;
float3 scene_aabbox_max;

        //OpenCL objects
        cl_context          context;
        cl_device_id        *devices;
        cl_command_queue    commandQueue;
        cl_program program;
        cl_kernel  kernel;

		CLCommandArgs   *sampleArgs;   /**< CLCommand argument class */
		SDKTimer    *sampleTimer;      /**< SDKTimer object */

		cl_bool reqdExtSupport;
        SDKDeviceInfo
        deviceInfo;            /**< Structure to store device information*/

		
cl_int ciErrNum;
cl_mem pbo_cl;

struct TRI {
	float4 v0, v1, v2;
} ; 

struct Params {

	float4 a; 
	float4 b; 
	float4 c; 
   float4 campos;
   float4 light_pos;
   float4 light_color;
   float4 scene_aabb_min; 
   float4 scene_aabb_max;

   Params ( ) {};

   void init ( float3 a, float3 b, float3 c, 
			 float3 campos,
		   float3 light_pos,
		   float3 light_color,
		   float3 scene_aabb_min, 
		   float3 scene_aabb_max ) {

	   /*
	   this-> a = a ; 
	   this-> b = b ; 
	   this-> c = c ; 

	   this-> campos = campos;

	   this-> light_pos = light_pos;
	   this-> light_color = light_color;

	   this-> scene_aabb_min = scene_aabb_min; 
	   this-> scene_aabb_max = scene_aabb_max;
	   */

	   this->a = float4(a.x, a.y, a.z, 1.0f);
	   this->b = float4(b.x, b.y, b.z, 1.0f);
	   this->c = float4(c.x, c.y, c.z, 1.0f);

	   this->campos = float4(campos.x, campos.y, campos.z, 1.0f);

	   this->light_pos = float4(light_pos.x, light_pos.y, light_pos.z, 1.0f);
	   this->light_color = float4(light_color.x, light_color.y, light_color.z, 1.0f);

	   this->scene_aabb_min = float4(scene_aabb_min.x, scene_aabb_min.y, scene_aabb_min.z, 1.0f);
	   this->scene_aabb_max = float4(scene_aabb_max.x, scene_aabb_max.y, scene_aabb_max.z, 1.0f);
   } ;
} ; 

Params params; 

struct RayTraceData {
	cl_mem source ; 
	//cl_mem dest ; 

	cl_mem params ; 

	cl_mem mesh_vertices ;
	cl_mem mesh_indices ;

	cl_mem bvh_nodes ;
	cl_mem bvh_tris_indices ; 

	cl_mem temp ; 

	// for smothing triangle reflections
	cl_mem mesh_normals;
	cl_mem mesh_normals_indices;

	cl_mem mesh_materials;
	cl_mem mesh_triangle_index_to_material_index;

	RayTraceData () {
		source = NULL ;
		//dest = NULL ;

		params = NULL ;

		mesh_vertices = NULL ;
		mesh_indices = NULL ;

		bvh_nodes = NULL ;
		bvh_tris_indices = NULL ;

		temp = NULL ;

		// for smothing triangle reflections
		mesh_normals = NULL;
		mesh_normals_indices = NULL;

		mesh_materials = NULL;
		mesh_triangle_index_to_material_index = NULL;
	}

	~RayTraceData () {
		// Free device memory 
		if ( NULL != source )	clReleaseMemObject( source ); 
		//if ( NULL != dest )		clReleaseMemObject( dest );
		if ( NULL != params )	clReleaseMemObject( params );

		if ( NULL != mesh_vertices )	clReleaseMemObject( mesh_vertices );
		if ( NULL != mesh_indices )	clReleaseMemObject( mesh_indices );

		if ( NULL != bvh_nodes )	clReleaseMemObject( bvh_nodes );
		if ( NULL != bvh_tris_indices )	clReleaseMemObject( bvh_tris_indices );

		if ( NULL != temp )	clReleaseMemObject( temp );

		// for smothing triangle reflections
		if (NULL != mesh_normals)	clReleaseMemObject(mesh_normals);
		if (NULL != mesh_normals_indices)	clReleaseMemObject(mesh_normals_indices);

		if (NULL != mesh_materials)	clReleaseMemObject(mesh_materials);
		if (NULL != mesh_triangle_index_to_material_index)	clReleaseMemObject(mesh_triangle_index_to_material_index);
	} ; 
} ; 

RayTraceData rtData ;


RayTracer::RayTracer(void)
{
}


RayTracer::~RayTracer(void)
{
}
///*
void
GLInit()
{
    glClearColor(0.0 ,0.0, 0.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);
    glClear(GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
}

void updateCamera();

void render();

void update_delta () {

	//update the delta time for animation
	static int lastFrameTime = 0;

	//static long t0 = GetTickCount () ;
	//static int frames = 0 ;

	if (lastFrameTime == 0)
	{
		lastFrameTime = glutGet(GLUT_ELAPSED_TIME);
	}

	int now = glutGet(GLUT_ELAPSED_TIME);
	int elapsedMilliseconds = now - lastFrameTime;

	//printf ( "\r now = %d, lastFrameTime = %d, elapsedMilliseconds = %d ", now, lastFrameTime, elapsedMilliseconds );
	//getchar();

	//if ( elapsedMilliseconds >= 1000 ) {
		delta_t = elapsedMilliseconds / 1000.0f ; // / 100; // delta_t = 0.01; //
		lastFrameTime = now;
	//}
} ; 

void draw_bvh ();

void displayGL () {

	//update_delta ();
	updateCamera();
	render();

	//draw_bvh ();

	glutSwapBuffers();
} ; 

void raytrace_cpu_pixel ( int x, int y ) ; 
void raytrace_cpu_pixel_bvh ( int x, int y );

void raytrace_cpu ( bool test_mode = false ) {

	//srand (time(NULL));

	if ( !test_mode ) {
		for ( int y =0; y<image_height; ++y ) {
			for ( int x =0; x<image_width; ++x ) {
				//raytrace_cpu_pixel ( x, y ) ; 
			//raytrace_cpu_pixel_bvh ( x, y ) ; 
				//out_data[y * image_width + x] = 0x005555;
				//out_data[y * image_width + x] = rand() % 256 << 16 | rand() % 256 << 8 | rand() % 256 << 0;
				//out_data[y * image_width + x] = frames << 16 | frames << 8 | frames << 0;
			}
		}

		raytrace_cpu_pixel_bvh ( 400, 400 ) ; 		

		return ; 
	}

	// OK
		for ( int y =517; y<517+300; ++y ) 
			for ( int x =563-150; x<563+100; ++x ) 
				raytrace_cpu_pixel ( x, y ) ;
} ; 

size_t shrRoundUp ( int loc, int w ) {
	return ( w / loc + ( w % loc > 0 ? 1 : 0 ) ) * loc;
} ; 

vector <int> temp;

int raytrace_gpgpu () {

	const size_t local_ws[] = {8, 8};
    const size_t global_ws[] = { shrRoundUp(local_ws[0], w), shrRoundUp(local_ws[1], h)};
	//const size_t global_ws[] = { 1024, 768 };

	ciErrNum = clEnqueueNDRangeKernel(commandQueue, kernel, 2, NULL, global_ws, local_ws, 0, 0, 0); // gridSize, localSize

	CHECK_OPENCL_ERROR(ciErrNum, "clEnqueueNDRangeKernel failed.)"); // status

    ciErrNum = clFinish(commandQueue); // status
    CHECK_OPENCL_ERROR(ciErrNum, "clFinish failed."); // status

	ciErrNum = clEnqueueReadBuffer(commandQueue, pbo_cl, CL_TRUE, 0, sizeof(unsigned int) * image_height * image_width, &out_data[0], 0, NULL, NULL);        
	CHECK_OPENCL_ERROR(ciErrNum, "clEnqueueReadBuffer." );


	//temp.resize( bvh_cuda.bvh_nodes.size() * 4 * sizeof ( int ) );

	//ciErrNum = clEnqueueReadBuffer(commandQueue, rtData.temp, CL_TRUE, 0, bvh_cuda.bvh_nodes.size() * sizeof ( BVH_Node_ ), &temp[0], 0, NULL, NULL);        
	//CHECK_OPENCL_ERROR(ciErrNum, "clEnqueueReadBuffer." );

	/*
	int flag1 = 0, flag2 = 0;
	for ( int i=0; i<bvh_cuda.bvh_nodes.size(); ++i ) {

		if ( temp[i*4+0] != bvh_cuda.bvh_nodes[i].offset_left || temp[i*4+1] != bvh_cuda.bvh_nodes[i].offset_right )
			flag1++;

		if ( temp[i*4+2] != bvh_cuda.bvh_nodes[i].offset_tris || temp[i*4+3] != bvh_cuda.bvh_nodes[i].num_tris )
			flag2++;
	}
	//printf ( "\n temp[0] = %d, temp[1] = %d", temp[0], temp[1] );
	printf ( "\n flag1 = %d, flag2 = %d", flag1, flag2 );
	getchar();
	*/

	/*
	printf ( "\n flag1 = %d", temp[0] );
	for ( int i=0; i<temp[0]; ++i ) {
		printf ( "\n stack_count = %d, node_offset = %d, offset_left = %d, offset_right = %d", temp[4+i*4], temp[4+i*4+1], temp[4+i*4+2], temp[4+i*4+3] );
		getchar();
	}
	*/

	/*
	printf( "\n ray->ori = (%.2f, %.2f, %.2f), tHit = %.2f, ray->dir = (%.2f, %.2f, %.2f)", 
		(float&)temp[4*1+0], (float&)temp[4*1+1], (float&)temp[4*1+2], (float&)temp[4*1+3],
		(float&)temp[4*2+0], (float&)temp[4*2+1], (float&)temp[4+2], (float&)temp[4*2+3]  );

	getchar();
	*/

	/*printf ( "\n flag1 = %d", temp[0] );
	for ( int i=0; i<temp[0]; ++i ) {
		printf ( "\n stack_count = %d, node_offset = %d, offset_left = %d, offset_right = %d, intersect0 = %d, intersect1 = %d, tspan0.x = %.2f, tspan1.x = %.2f ", 
			temp[4+i*4*2], temp[4+i*4*2+1], temp[4+i*4*2+2], temp[4+i*4*2+3], 
			temp[4+i*4*2+4*1 + 0], temp[4+i*4*2+4*1 + 1], (float&)temp[4+i*4*2+4*1 + 2], (float&)temp[4+i*4*2+4*1 + 3]	);
		getchar();
	}*/

	return SDK_SUCCESS;
} ; 

void draw_bvh();

void render()
{
	glClear(GL_COLOR_BUFFER_BIT);

    glBindTexture(GL_TEXTURE_2D, texnum);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);

    unsigned char bitmap[WIDTH * HEIGHT * 4]; // rgba unsigned bytes
	//unsigned char bitmap[1024 * 768 * 4]; // rgba unsigned bytes

	//bitmap[(1 + 1 * WIDTH) * 4 + 0] = 0;
	//bitmap[(4 + 3* WIDTH) * 4 + 0] = 0;

	//for(int y = 0; y < 1; y++) // 4	
	//	for(int x = 0; x < 10; x++) // 5
	//		bitmap[(x + y* 1024) * 4 + 0] = 0;

    //double m, r, g, b;

    /*for(int y = 0; y < HEIGHT; y++)
    {
        for(int x = 0; x < WIDTH; x++)
        {
            //r = g = b = 128;

			//bitmap[(x + y * WIDTH) * 4 + 0] = 0;

            //bitmap[(x + y * WIDTH) * 4 + 0] = (unsigned char)(r * 255);
            //bitmap[(x + y * WIDTH) * 4 + 1] = (unsigned char)(g * 255);
            //bitmap[(x + y * WIDTH) * 4 + 2] = (unsigned char)(b * 255);
            //bitmap[(x + y * WIDTH) * 4 + 3] = 255;
        }
    }*/

	long t0 = GetTickCount () ; 

	//raytrace_cpu(false);
	raytrace_gpgpu ();

	long t1 = GetTickCount () ; 

	//printf ( "\r render time = %d ms", (t1 - t0) );

    //glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, bitmap);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, WIDTH, HEIGHT, 0, GL_RGBA, GL_UNSIGNED_BYTE, (unsigned char *)&out_data[0]);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glEnable(GL_TEXTURE_2D);

    // setup 2d pixel plotting camera
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0f, (GLdouble) winwidth, 0.0f, (GLdouble) winheight, 0.0f, 1.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0, 0, winwidth, winheight);

    glBegin(GL_QUADS);

    glColor3f(1.0, 1.0, 1.0);

    glTexCoord2f(0.0, 0.0);
    glVertex2i(0, 0);

    glTexCoord2f(1.0, 0.0);
    glVertex2i(winwidth, 0);

    glTexCoord2f(1.0, 1.0);
    glVertex2i(winwidth, winheight);

    glTexCoord2f(0.0, 1.0);
    glVertex2i(0, winheight);

    glEnd();

    glDisable(GL_TEXTURE_2D);

    glFlush();
    //glutSwapBuffers();
}

void timerEvent(int value)
{
    glutPostRedisplay();
	glutTimerFunc(REFRESH_DELAY, timerEvent,0);
}

void update()
{
    //me->runCLKernels();

    // redraw
    glutPostRedisplay();

	if ( 0 == t0 ) {
		t0 = glutGet(GLUT_ELAPSED_TIME);

		delta_t = 0.01f;

        frames = 0;

		return ; 
	} ;

    frames++;

    te = glutGet(GLUT_ELAPSED_TIME);

	//printf ( "\n te = %d, t0 = %d, elapsedMilliseconds = ", te, t0 );
	//getchar();

    // every second approximately
    if (te - t0 >= 1000)
    {
        char title[250];
        sprintf(title, "Zabolotnov Vadim Alexandrovich ( itmanager85@ya.ru ),   Ray trace demo    %.1f fps", (1000.0*frames/(te-t0)));
        glutSetWindowTitle(title);

		delta_t = (te - t0) / 1000.0f / frames;

        frames = 0;
        t0 = te;
    }

	//printf ( "\r update, frames = %d", frames );

	//camera_rotation += 1.25f * 0.25f/ 100.0f;
}

int setupCL();
int cleanup();

void init_arg1 ();
void initCLVolume2();

void initRayTrace ();

int button_curr = 0, state_curr = 0;
int x_curr = 0, y_curr = 0;

void mouse( int button, int state, int x, int y ) {
	button_curr = button;
	state_curr = state;
	x_curr = x; 
	y_curr = y;

	/*if ((button == 3) || (button == 4)) // It's a wheel event
    {
		printf ("\n test !! wheel ");
       // Each wheel event reports like a button click, GLUT_DOWN then GLUT_UP
       if (GLUT_DOWN == state) {
		   camera_distance += (3 == button) ? -5 : 5 ;
	   }
	}*/
} ; 

//void wheel(int wheel, int direction, int x, int y) {
//	printf ("\n test !! wheel ");
//} ; 

void motion( int x, int y ) {
	if ( GLUT_RIGHT_BUTTON == button_curr && GLUT_DOWN == state_curr ) {
		camera_rotation += (x - x_curr) * 0.25f/ 100.0f;

		cam.add_rotate((x - x_curr) * 0.25f / 100.0f, (y - y_curr) * 0.25f / 100.0f);

	} else if ( GLUT_LEFT_BUTTON == button_curr && GLUT_DOWN == state_curr ) {
		camera_distance += (y - y_curr) ;

		cam.add_radius(y - y_curr);
	}

	x_curr = x; 
	y_curr = y;
} ; 

int
main(int argc, char * argv[])
{
	if ( SDK_SUCCESS != setupCL() ) {
		getchar();
		return -1;	
	}
	initRayTrace ();

	init_arg1 ();
	initCLVolume2();

	// Run in  graphical window if requested
	glutInit(&argc, argv);
	glutInitWindowPosition(100,10);
	glutInitWindowSize ( WIDTH, HEIGHT ) ; // (1024,768);
	glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE );
	glutCreateWindow("Ray Tracer ");
	GLInit();
	glutDisplayFunc(displayGL ); //glutDisplayFunc(render);
	glutIdleFunc(update);
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	//glutMouseWheelFunc(wheel);
	//glutReshapeFunc(reshape);
	//glutKeyboardFunc(keyboard);

	//glutTimerFunc(REFRESH_DELAY, timerEvent,0);

	glutMainLoop();

	cleanup();
}
//*/

void updateCamera()
{
	
	/*
	campos = make_float3( 
		cos(camera_rotation) * camera_distance,
		camera_height,
		-sin(camera_rotation) * camera_distance
		);

	//camera.pos = campos ; 

	//printf ( "\n campos = ( %.2f, %.2f, %.2f ) ", campos.x, campos.y,campos.z );
	//printf ( "\n cam.pos = ( %.2f, %.2f, %.2f ) ", camera.pos.x, camera.pos.y, camera.pos.z );

	float3 cam_dir = -campos;
	cam_dir = normalize(cam_dir);
	float3 cam_up  = make_float3(0,1,0);
	float3 cam_right = cross(cam_dir,cam_up);
	cam_right = normalize(cam_right);

	cam_up = -cross(cam_dir,cam_right);
	cam_up = normalize(cam_up);
	*/

	campos = cam.eye;
	float3 cam_dir = cam.camera_direction;
	float3 cam_right = cam.camera_right;
	float3 cam_up = cam.camera_up;
	
	float FOV = 60.0f;
	float theta = (FOV*3.1415*0.5) / 180.0f;
	float half_width = tanf(theta);
	float aspect = (float)image_width / (float)image_height;

	float u0 = -half_width * aspect;
	float v0 = -half_width;
	float u1 =  half_width * aspect;
	float v1 =  half_width;
	float dist_to_image = 1;

	a = (u1-u0)*cam_right;
	b = (v1-v0)*cam_up;
	c = campos + u0*cam_right + v0*cam_up + dist_to_image*cam_dir;
	
	if(animate) {
		camera_rotation += 0.25f * delta_t;

		cam.add_rotate( -0.25f*1.5f * delta_t, 0.0f );
		//cam.add_rotate(-0.25f*1.5f / 100.0f, 0.0f);
	}

	//printf ( "\n camera_rotation = %.2f ", camera_rotation );

	//a = cam.camera_right;
	//b = cam.camera_up;
	//c = cam.camera_direction;

	//campos = cam.eye;

	params.init ( a, b, c, campos, light_pos, light_color1, scene_aabbox_min, scene_aabbox_max ); 

	clEnqueueWriteBuffer( commandQueue, rtData.params, CL_FALSE, 0, sizeof( Params ), &params, 0, NULL, NULL);
}


void raytrace_cpu_pixel ( int x, int y ) {

	const float xf = (x-0.5)/((float)image_width); // w
	const float yf = (y-0.5)/((float)image_height); // h

	int ray_depth = 0;
	bool continue_path = true;

	float3 t1 = c+(a*xf);
	float3 t2 = b*yf;
	float3 image_pos = t1 + t2;

	spec :: Ray r( image_pos,image_pos - campos );
	//r.print () ; 

	spec:: HitRecord hit_r;

	float t_min,t_max;
	continue_path = spec:: RayBoxIntersection(scene_aabbox_min, scene_aabbox_max, r.o, r.inv_dir,t_min, t_max);
	hit_r.color = make_float3(0,0,0);

	//float3 light_pos =  make_float3(light_x,light_y,light_z) ; 

	// hack to display the light source we simple make a ray sphere intersection and 
	// compare the depth with the found t value from the triangles
	float sphere_t;
	bool sphere_hit = spec:: RaySphereIntersection(r, light_pos,2.0,sphere_t); 

	//float3 light_color1 = make_float3( light_color[0], light_color[1], light_color[2] ) ; 
	
	if(sphere_hit && sphere_t > 0.001)
	{
		if(!continue_path)
		{
			hit_r.color =  light_color1;
		}
		sphere_hit = true;
	}

	int count = 0 ; 

	while(continue_path && ray_depth < 1) // 4 
	{
			for(int i = 0; i < total_number_of_triangles; i++)
			{
				//atomicAdd ( &rt_p.pixels_recomputed[2], 1 ) ; 
				//if ( !b_is_recompute && 1 == count && index_tri == i ) continue ; 

				float4 v0 = triangles [i*3]; // make_float4 ( 1,1,1,1 ) ; // 
				float4 e1 = triangles [i*3+1]; // make_float4 ( 1,1,1,1 ) ; // 
				float4 e2 = triangles [i*3+2]; // make_float4 ( 1,1,1,1 ) ; // 

				float t = spec:: RayTriangleIntersection(r, make_float3(v0.x,v0.y,v0.z),make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));

				if(t < hit_r.t && t > 0.001)
				{
					hit_r.t = t; 
					hit_r.hit_index = i;

					//if ( 1 == count )
					//	pvdata [ y * w + x ] = i ; 
					//pm_triangle_index
					//pm_lvl_triangle_index [ count ] = i ; 

					//find = true ; 

					//break ; 
				}
			}

		if(sphere_hit && sphere_t < hit_r.t) // if(sphere_hit && sphere_t < hit_r_t) // 
		{
			hit_r.color += light_color1;
			continue_path = false;
			break;
		}

		if(hit_r.hit_index >= 0)
		{

			ray_depth++;
			
			// create the normal
			float4 e1 = triangles [ hit_r.hit_index*3+1 ];
			float4 e2 = triangles [ hit_r.hit_index*3+2 ];

			hit_r.normal = cross(make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));
			//hit_r.normal = cross( triangles [ hit_r.hit_index*3+1 ], triangles [ hit_r.hit_index*3+2 ] ) ; 

			hit_r.normal = normalize(hit_r.normal);

			/*if ( width/2 == x && height/2 == y ) {
				printf ( "\n float4 = (%.2f, %.2f, %.2f )", e1.x, e1.y, e1.z );
				printf ( "\n float4 = (%.2f, %.2f, %.2f )", e2.x, e2.y, e2.z );
				printf ( "\n float4 = (%.2f, %.2f, %.2f )", hit_r.normal.x, hit_r.normal.y, hit_r.normal.z );
			}*/

			// calculate simple diffuse light
			float3 hitpoint = r.o + r.dir *hit_r.t;
			float3 L = light_pos - hitpoint; // rayTraceExtrime
			float dist_to_light = length(L);
			
			L = normalize(L);
			float diffuse_light = fmaxf1( dot(L,hit_r.normal), 0.0);
			diffuse_light = fminf1( (diffuse_light),1.0);

			/*if ( width/2 == x && height/2 == y ) {
				printf ( "\n float = (%.2f)", dot(L, hit_r.normal) );
				printf ( "\n float2 = (%.2f, %.2f)", dist_to_light, diffuse_light );
				getchar();
			}*/

			//calculate simple specular light
			float3 H = L + (-r.dir);
			H = normalize(H);
			float specular_light = powf(fmaxf1(dot(H,hit_r.normal),0.0),25.0f);

			diffuse_light  *=  16.0/dist_to_light;
			specular_light *=  16.0/dist_to_light;

			clamp(diffuse_light, 0.0f, 1.0f);
			clamp(specular_light, 0.0f, 1.0f);

			// rayTraceExtrime
			hit_r.color += light_color1 * diffuse_light + make_float3(1.0,1.0,1.0)*specular_light*0.2 + make_float3(0.2,0.2,0.2); 
		}
		else
		{
			continue_path = false;
			hit_r.color += make_float3(0.5,0.5,0.95*yf+0.3);
		}
	}

	if(ray_depth >= 1 || sphere_hit)
	{
		ray_depth = imax1(ray_depth,1);
		hit_r.color /= ray_depth; // normalize the colors
	}
	else
	{
		hit_r.color = make_float3(0.5,0.5,yf+0.3);
	}

	int val = spec:: rgbToInt(hit_r.color.x*255,hit_r.color.y*255,hit_r.color.z*255);
	out_data[y * image_width + x] = val;
} ; 


void loadObj(const std::string filename, TriangleMesh &mesh);

void init_bvh_test0 () {

	mesh1.add_tri ( make_float4( -10, 0, 0, 1.0f ), make_float4( -10, 10, 0, 1.0f ), make_float4( -1, 0, 0, 1.0f ) );
	mesh1.add_tri ( make_float4( 10, 0, 0, 1.0f ), make_float4( 10, 10, 0, 1.0f ), make_float4( 1, 0, 0, 1.0f ) );

	mesh1.add_tri ( make_float4( 0, 0, -10, 1.0f ), make_float4( 0, 10, -10, 1.0f ), make_float4( 0, 0, -1, 1.0f ) );
	mesh1.add_tri ( make_float4( 0, 0, 10, 1.0f ), make_float4( 0, 10, 10, 1.0f ), make_float4( 0, 0, 1, 1.0f ) );

	//mesh1.add_tri ( make_float4( 10, -5, 0, 1.0f ), make_float4( 10, -5, 10, 1.0f ), make_float4( 0, -5, 0, 1.0f ) );
	//mesh1.add_tri ( make_float4( 10, -5, 10, 1.0f ), make_float4( 0, -5, 10, 1.0f ), make_float4( 0, -5, 0, 1.0f ) );
	mesh1.add_tri ( make_float4( -10, -5, -10, 1.0f ), make_float4( -10, -5, 10, 1.0f ), make_float4( 10, -5, -10, 1.0f ) );
	mesh1.add_tri ( make_float4( -10, -5, 10, 1.0f ), make_float4( 10, -5, 10, 1.0f ), make_float4( 10, -5, -10, 1.0f ) );

	//bvh.build( mesh1 ); 
} ; 

void init_bvh_test1 () {

	
	mesh1.add_tri ( make_float4( -10, 0, 0, 1.0f ), make_float4( -10, 10, 0, 1.0f ), make_float4( 0, 0, 0, 1.0f ) );
	mesh1.add_tri ( make_float4( 10, 0, 0, 1.0f ), make_float4( 10, 10, 0, 1.0f ), make_float4( 0, 0, 0, 1.0f ) );

	mesh1.add_tri ( make_float4( 0, 0, -10, 1.0f ), make_float4( 0, 10, -10, 1.0f ), make_float4( 0, 0, 0, 1.0f ) );
	 //mesh1.add_tri ( make_float4( 0, 0, 10, 1.0f ), make_float4( 0, 10, 10, 1.0f ), make_float4( 0, 0, 0, 1.0f ) );
	
	///*
	mesh1.add_tri ( make_float4( -10, 0, -10, 1.0f ), make_float4( -10, 0, 10, 1.0f ), make_float4( 10, 0, -10, 1.0f ) );
	mesh1.add_tri ( make_float4( -10, 0, 10, 1.0f ), make_float4( 10, 0, 10, 1.0f ), make_float4( 10, 0, -10, 1.0f ) );
	//*/

	//bvh.build( mesh1 ); 
} ; 

void initRayTrace () {

	ColladaLoader loader;

	loader.load("data/collada/cubes2.DAE");
	//loader.load("G:/Dev/collada/test/cubes2_update1.DAE");

	 //int size = window_width * window_height ; 
	//int size = image_width * image_height ; 

	// next we load a simple obj file and upload the triangles to an 1D texture.
	//loadObj("data/cube.obj",mesh);
	//loadObj("data/models/cubes_v02.obj",mesh);
	//loadObj("data/models/cubes2.obj", mesh); // with vertex normals
	//loadObj("data/models/snake_head_v02.obj",mesh);
	//loadObj("data/conference.obj",mesh); // conference // cube // cubes // 
	//loadObj("data/models/chess_v1.obj",mesh);
//	loadObj("data/models/chess2.obj",mesh);
//	loadObj("data/sphere.obj",sphere);

	//vector<float4> triangles;

//	printf ( "\n mesh.faces.size() = %d", mesh.faces.size() );
//	printf ( "\n sphere.faces.size() = %d \n", sphere.faces.size() );

//	mesh1.init ( mesh );
	mesh1.init( loader );
	printf("\n mesh1.faces.size() = %d \n", mesh1.indices.size()/3);

	//init_bvh_test1 ();
//	bvh.build( mesh1 ); 
	//bvh_cuda.build ( bvh.root );	

	bvh2.setMesh( &mesh1 );
	bvh_cuda.build_from_bvh2 ( bvh2 );

	//AABB aabb;
	//printf ( "\n sizeof (AABB aabb) = %d", sizeof (aabb) );

	//BVH_Node_ test;
	//printf ( "\n sizeof (BVH_Node_ test) = %d", sizeof (test) );

	//printf ( "\n bvh.build( mesh1 ) " );
	//printf ( "\n depth = %d, root, tris.size() = %d", 0, bvh.root->tris.size() );
	//printf ( "\n depth = %d, left, tris.size() = %d", 1, bvh.root->pLeft->tris.size() );
	//printf ( "\n depth = %d, right, tris.size() = %d", 1, bvh.root->pRight->tris.size() );

	//if ( false )
	for(unsigned int i = 0; i < mesh.faces.size(); i++)
	{
		float3 v0 = mesh.verts[mesh.faces[i].v[0]-1];
		float3 v1 = mesh.verts[mesh.faces[i].v[1]-1];
		float3 v2 = mesh.verts[mesh.faces[i].v[2]-1];
		triangles.push_back(make_float4(v0.x,v0.y,v0.z,0));
		triangles.push_back(make_float4(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z,0)); // notice we store the edges instead of vertex points, to save some calculations in the 
		triangles.push_back(make_float4(v2.x-v0.x, v2.y-v0.y, v2.z-v0.z,0)); // ray triangle intersection test.
	}

	//if ( false )
	for(unsigned int i = 0; i < sphere.faces.size(); i++)
	{
		float3 v0 = sphere.verts[sphere.faces[i].v[0]-1];
		float3 v1 = sphere.verts[sphere.faces[i].v[1]-1];
		float3 v2 = sphere.verts[sphere.faces[i].v[2]-1];
		triangles.push_back(make_float4(v0.x,v0.y,v0.z,0));
		triangles.push_back(make_float4(v1.x-v0.x, v1.y-v0.y, v1.z-v0.z,1)); // notice we store the edges instead of vertex points, to save some calculations in the 
		triangles.push_back(make_float4(v2.x-v0.x, v2.y-v0.y, v2.z-v0.z,0)); // ray triangle intersection test.
	}

//	cout << "\ntotal number of triangles check:" << mesh.faces.size() + sphere.faces.size() << " == " << triangles.size()/3 << endl;

	size_t triangle_size = triangles.size() * sizeof(float4);
	total_number_of_triangles = triangles.size()/3;

	//printf ( "\n triangles.size() = %d \n ", triangles.size() );
	//printf ( "\n triangles[0] = %.2f \n ", triangles[0] );
	//printf ( "\n triangles[total_number_of_triangles-1] = %.2f \n ", triangles[total_number_of_triangles-1] );

	if (total_number_of_triangles > 0) //if(triangle_size > 0)
	{
		//cutilSafeCall( cudaMalloc((void **)&dev_triangle_p, triangle_size));
		//cudaMemcpy(dev_triangle_p,&triangles[0],triangle_size,cudaMemcpyHostToDevice);
		//bindTriangles(dev_triangle_p, total_number_of_triangles);

		rtData.source = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, triangle_size, &triangles[0], 0);

		/*vector <float3> mesh_verts;
		mesh_verts.reserve( mesh.verts.size() + sphere.verts.size() );
		mesh_verts.insert( mesh_verts.end(), mesh.verts.begin(), mesh.verts.end() );
		mesh_verts.insert( mesh_verts.end(), sphere.verts.begin(), sphere.verts.end() );

		rtData.mesh_vertices = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh_verts.size() * sizeof(float3), &mesh_verts[0], 0);

		vector <float3> mesh_indices;
		mesh_indices.reserve( mesh.faces.size() + sphere.faces.size() );
		mesh_indices.insert( mesh_verts.end(), mesh.faces.begin(), mesh.faces.end() );
		mesh_indices.insert( mesh_verts.end(), sphere.faces.begin(), sphere.faces.end() );

		rtData.mesh_indices = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh_verts.size() * sizeof(int), &mesh_indices[0], 0);*/
	}

	if ( mesh1.indices.size() >= 3 ) { 

		rtData.mesh_vertices = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh1.vertices.size() * sizeof(float4), &mesh1.vertices[0], 0);
		rtData.mesh_indices = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh1.indices.size() * sizeof(int), &mesh1.indices[0], 0);

		rtData.bvh_nodes = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, bvh_cuda.bvh_nodes.size() * sizeof(BVH_Node_), &bvh_cuda.bvh_nodes[0], 0);
		rtData.bvh_tris_indices = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, bvh_cuda.tri_indices.size() * sizeof(int), &bvh_cuda.tri_indices[0], 0);

		if (mesh1.normals_indices.size() > 0) {
			//printf("\n Good normals!!");
			//getchar();

			rtData.mesh_normals = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh1.normals.size() * sizeof(float4), &mesh1.normals[0], 0);
			rtData.mesh_normals_indices = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh1.normals_indices.size() * sizeof(int), &mesh1.normals_indices[0], 0);
		}

		if (mesh1.triangle_index_to_material_index.size() > 0) {
			rtData.mesh_materials = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh1.materials.size() * sizeof(Material), &mesh1.materials[0], 0);
			rtData.mesh_triangle_index_to_material_index = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, mesh1.triangle_index_to_material_index.size() * sizeof(int), &mesh1.triangle_index_to_material_index[0], 0);
		}
	}

	//rtData.temp = clCreateBuffer( context, CL_MEM_WRITE_ONLY, bvh_cuda.bvh_nodes.size() * sizeof(BVH_Node_), 0, 0);
	rtData.temp = clCreateBuffer(context, CL_MEM_WRITE_ONLY, (64+1)*sizeof(float4), 0, 0);

	rtData.params = clCreateBuffer( context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, sizeof(Params), &params, 0);

	//float3 mid3 = ( mesh.bounding_box[0] + mesh.bounding_box[1] ) / 2;
//	scene_aabbox_min = mesh.bounding_box[0]; //+ mid3;
//	scene_aabbox_max = mesh.bounding_box[1]; //+ mid3;

	scene_aabbox_min = mesh1.scene_aabbox_min;
	scene_aabbox_max = mesh1.scene_aabbox_max;

	//light_pos = ( mesh.bounding_box[0] + mesh.bounding_box[1] ) / 1.5f;

	/*
	scene_aabbox_min.x = min(scene_aabbox_min.x,sphere.bounding_box[0].x);
	scene_aabbox_min.y = min(scene_aabbox_min.y,sphere.bounding_box[0].y);
	scene_aabbox_min.z = min(scene_aabbox_min.z,sphere.bounding_box[0].z);

	scene_aabbox_max.x = max(scene_aabbox_max.x,sphere.bounding_box[1].x);
	scene_aabbox_max.y = max(scene_aabbox_max.y,sphere.bounding_box[1].y);
	scene_aabbox_max.z = max(scene_aabbox_max.z,sphere.bounding_box[1].z);
	*/
	//rtData.dest = clCreateBuffer( cxGPUContext, CL_MEM_WRITE_ONLY, size, 0, 0);
} ; 

// load a simple obj file without normals or tex-coords
void loadObj( const std::string filename, TriangleMesh &mesh )
{

	std::ifstream in(filename.c_str());

	if(!in.good())
	{
		cout  << "ERROR: loading obj:(" << filename << ") file is not good" << "\n";
		getchar();
		exit(0);
	}

	char buffer[256], str[255];
	float f1,f2,f3;

	while(!in.getline(buffer,255).eof())
	{
		buffer[255]='\0';

		sscanf_s(buffer,"%s",str,255);

		if (buffer[0] == 'v' && buffer[1] == 'n' && (buffer[2] == ' ' || buffer[2] == 32)) {

			if (sscanf(buffer, "vn %f %f %f", &f1, &f2, &f3) == 3) {
				mesh.norms.push_back(make_float3(f1, f2, f3));
				//mesh.norms.push_back(make_float3(1, 1, 1));

			}
			else {
				cout << "ERROR: normals not in wanted format in OBJLoader" << "\n";
				getchar();
				exit(-1);
			}

		// reading a vertex
		} else if (buffer[0]=='v' && (buffer[1]==' '  || buffer[1]==32) ) {
		
			if ( sscanf(buffer,"v %f %f %f",&f1,&f2,&f3)==3)
			{
				mesh.verts.push_back(make_float3(f1,f2,f3));
			}
			else
			{
				cout << "ERROR: vertex not in wanted format in OBJLoader" << "\n";
				getchar();
				exit(-1);
			}

		// reading FaceMtls 
		} else if (buffer[0]=='f' && (buffer[1]==' ' || buffer[1]==32) )
		{
			TriangleFace f;
			int nt = sscanf(buffer,"f %d %d %d",&f.v[0],&f.v[1],&f.v[2]);

			//int3 vn = (int3)(0,0,0);

			if (nt != 3) {

				// for smothing triangle reflections
				nt = sscanf(buffer, "f %d//%d %d//%d %d//%d", &f.v[0], &f.vn[0], &f.v[1], &f.vn[1], &f.v[2], &f.vn[2]);

				if (nt != 6) {
					cout << "ERROR: I don't know the format of that FaceMtl" << "\n";
					getchar();
					exit(-1);
				}
			} 

			mesh.faces.push_back(f);
		}
	}

	// calculate the bounding box
	mesh.bounding_box[0] = make_float3(1000000,1000000,1000000);
	mesh.bounding_box[1] = make_float3(-1000000,-1000000,-1000000);
	for(unsigned int i = 0; i < mesh.verts.size(); i++)
	{
		//update min value
		mesh.bounding_box[0].x = imin1(mesh.verts[i].x,mesh.bounding_box[0].x);
		mesh.bounding_box[0].y = imin1(mesh.verts[i].y,mesh.bounding_box[0].y);
		mesh.bounding_box[0].z = imin1(mesh.verts[i].z,mesh.bounding_box[0].z);

		//update max value
		mesh.bounding_box[1].x = imax1(mesh.verts[i].x,mesh.bounding_box[1].x);
		mesh.bounding_box[1].y = imax1(mesh.verts[i].y,mesh.bounding_box[1].y);
		mesh.bounding_box[1].z = imax1(mesh.verts[i].z,mesh.bounding_box[1].z);

	}

	cout << "obj file loaded: number of faces:" << mesh.faces.size() << " number of vertices:" << mesh.verts.size() << " number of normals:" << mesh.norms.size() << endl;
	cout << "obj bounding box: min:(" << mesh.bounding_box[0].x << "," << mesh.bounding_box[0].y << "," << mesh.bounding_box[0].z <<") max:" 
		<< mesh.bounding_box[1].x << "," << mesh.bounding_box[1].y << "," << mesh.bounding_box[1].z <<")" << endl;
}

int setupCL2()
{
	sampleArgs = new CLCommandArgs() ;
    sampleTimer = new SDKTimer();
    sampleArgs->sampleVerStr = SAMPLE_VERSION;

    cl_int status = CL_SUCCESS;
    cl_device_type dType;

    if(sampleArgs->deviceType.compare("cpu") == 0)
    {
        dType = CL_DEVICE_TYPE_CPU;
    }
    else //sampleArgs->deviceType = "gpu"
    {
        dType = CL_DEVICE_TYPE_GPU;
        if(sampleArgs->isThereGPU() == false)
        {
            std::cout << "GPU not found. Falling back to CPU device" << std::endl;
            dType = CL_DEVICE_TYPE_CPU;
        }
    }

    //
    // Have a look at the available platforms and pick either
    // the AMD one if available or a reasonable default.
    //
    cl_platform_id platform = NULL;
    int retValue = getPlatform(platform, sampleArgs->platformId,
                               sampleArgs->isPlatformEnabled());
    CHECK_ERROR(retValue, SDK_SUCCESS, "getPlatform() failed");

    // Display available devices.
    retValue = displayDevices(platform, dType);
    CHECK_ERROR(retValue, SDK_SUCCESS, "displayDevices() failed");


    // If we could find our platform, use it. Otherwise use just available platform.

    cl_context_properties cps[3] =
    {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)platform,
        0
    };

    context = clCreateContextFromType(
                  cps,
                  dType,
                  NULL,
                  NULL,
                  &status);
    CHECK_OPENCL_ERROR( status, "clCreateContextFromType failed.");

    // getting device on which to run the sample
    status = getDevices(context, &devices, sampleArgs->deviceId,
                        sampleArgs->isDeviceIdEnabled());
    CHECK_ERROR(status, SDK_SUCCESS, "getDevices() failed");

    {
        // The block is to move the declaration of prop closer to its use
        cl_command_queue_properties prop = 0;
        commandQueue = clCreateCommandQueue(
                           context,
                           devices[sampleArgs->deviceId],
                           prop,
                           &status);
        CHECK_OPENCL_ERROR( status, "clCreateCommandQueue failed.");
    }

    //Set device info of given cl_device_id
    retValue = deviceInfo.setDeviceInfo(devices[sampleArgs->deviceId]);
    CHECK_ERROR(retValue, 0, "SDKDeviceInfo::setDeviceInfo() failed");


    std::string buildOptions = std::string("");
    // Check if cl_khr_fp64 extension is supported
    if(strstr(deviceInfo.extensions, "cl_khr_fp64"))
    {
        buildOptions.append("-D KHR_DP_EXTENSION");
    }
    else
    {
        // Check if cl_amd_fp64 extension is supported
        if(!strstr(deviceInfo.extensions, "cl_amd_fp64"))
        {
            reqdExtSupport = false;
            OPENCL_EXPECTED_ERROR("Device does not support cl_amd_fp64 extension!");
        }
    }

	// create a CL program using the kernel source
    buildProgramData buildData;
    buildData.kernelName = std::string("volumeRender.cl"); // buildData.kernelName = std::string("FluidSimulation2D_Kernels.cl");
    buildData.devices = devices;
    buildData.deviceId = sampleArgs->deviceId;
    buildData.flagsStr = std::string("");
    if(sampleArgs->isLoadBinaryEnabled())
    {
        buildData.binaryName = std::string(sampleArgs->loadBinary.c_str());
    }

    if(sampleArgs->isComplierFlagsSpecified())
    {
        buildData.flagsFileName = std::string(sampleArgs->flags.c_str());
    }

    retValue = buildOpenCLProgram(program, context, buildData);
    CHECK_ERROR(retValue, 0, "buildOpenCLProgram() failed");

    // get a kernel object handle for a kernel with the given name
    kernel = clCreateKernel(
                 program,
                 "raytracer_bvh", //"raytracer", // "lbm"
                 &status);
    CHECK_OPENCL_ERROR(status, "clCreateKernel failed.");

	return SDK_SUCCESS;
}

void init_arg1 () {

	ciErrNum = CL_SUCCESS;

	pbo_cl = clCreateBuffer( context, CL_MEM_WRITE_ONLY, image_width * image_height * sizeof(GLubyte) * 4, NULL, &ciErrNum);

    ciErrNum |= clSetKernelArg ( kernel, 0, sizeof(cl_mem), (void *) &pbo_cl ) ;
    ciErrNum |= clSetKernelArg ( kernel, 1, sizeof(unsigned int), &image_width );
    ciErrNum |= clSetKernelArg ( kernel, 2, sizeof(unsigned int), &image_height );
} ; 


void initCLVolume2() {

	cl_int error;

		error = clSetKernelArg( kernel, 3, sizeof(cl_mem), &rtData.source ); assert(!error);
		error = clSetKernelArg( kernel, 4, sizeof(cl_int), &total_number_of_triangles ); assert(!error);
		error = clSetKernelArg( kernel, 5, sizeof(cl_mem), &rtData.params ); assert(!error);
		
		error = clSetKernelArg( kernel, 6, sizeof(cl_mem), &rtData.mesh_vertices ); assert(!error);
		error = clSetKernelArg( kernel, 7, sizeof(cl_mem), &rtData.mesh_indices ); assert(!error);

		error = clSetKernelArg( kernel, 8, sizeof(cl_mem), &rtData.bvh_nodes ); assert(!error);
		error = clSetKernelArg( kernel, 9, sizeof(cl_mem), &rtData.bvh_tris_indices ); assert(!error);

		int num_bvh_tris = bvh_cuda.tri_indices.size();
		error = clSetKernelArg( kernel, 10, sizeof(cl_int), &num_bvh_tris ); assert(!error);

		int num_bvh_nodes = bvh_cuda.bvh_nodes.size();
		error = clSetKernelArg( kernel, 11, sizeof(cl_int), &num_bvh_nodes ); assert(!error);

		error = clSetKernelArg( kernel, 12, sizeof(cl_mem), &rtData.temp ); assert(!error);

		error = clSetKernelArg(kernel, 13, sizeof(cl_mem), &rtData.mesh_normals); assert(!error);
		error = clSetKernelArg(kernel, 14, sizeof(cl_mem), &rtData.mesh_normals_indices); assert(!error);

		error = clSetKernelArg(kernel, 15, sizeof(cl_mem), &rtData.mesh_materials); assert(!error);
		error = clSetKernelArg(kernel, 16, sizeof(cl_mem), &rtData.mesh_triangle_index_to_material_index); assert(!error);
}

int cleanup()
{
	delete sampleArgs;   /**< CLCommand argument class */
	delete sampleTimer;      /**< SDKTimer object */

    // Releases OpenCL resources (Context, Memory etc.)
    cl_int status;

    status = clReleaseKernel(kernel);
    CHECK_OPENCL_ERROR(status, "clReleaseKernel failed.(kernel)");

    status = clReleaseProgram(program);
    CHECK_OPENCL_ERROR(status, "clReleaseProgram failed.(program)");

    //status = clReleaseMemObject(d_if0);
    //CHECK_OPENCL_ERROR(status, "clReleaseMemObject failed.(d_if0)");

    status = clReleaseCommandQueue(commandQueue);
    CHECK_OPENCL_ERROR(status, "clReleaseCommandQueue failed.(commandQueue)");

    status = clReleaseContext(context);
    CHECK_OPENCL_ERROR(status, "clReleaseContext failed. (context)");

    return SDK_SUCCESS;
}

float2 ray_box ( spec::Ray &r, AABB & box ) {

	float3 t0 = ( make_float3 ( box.min.x, box.min.y, box.min.z ) - r.o) / r.dir;
    float3 t1 = ( make_float3 ( box.max.x, box.max.y, box.max.z ) - r.o) / r.dir;

	float3 tmin = fminf1(t0,t1);
    float3 tmax = fmaxf1(t0,t1);

    float tmin1 = fmaxf1( fmaxf1( tmin.x, tmin.y ), tmin.z );
    float tmax1 = fminf1( fminf1( tmax.x, tmax.y ), tmax.z );

	return make_float2( tmin1, tmax1 );
} ; 

//template <class T> void swap1(T& a, T& b) { T t = a; a = b; b = t; }

const float TMIN = 0.001;

bool traceRecursive ( Mesh & mesh, BVH_Node * pNode, spec :: Ray & ray, float &tHit, TRI & tri ) { // find_tri

	//printf ( "\n pNode->tris.size() = %d ", pNode->tris.size() );
	//printf ( "\n bvh.build( mesh1 ) " );
	//printf ( "\n depth = %d, root, tris.size() = %d", 0, bvh.root->tris.size() );
	//printf ( "\n depth = %d, left, tris.size() = %d", 1, bvh.root->pLeft->tris.size() );
	//printf ( "\n depth = %d, right, tris.size() = %d", 1, bvh.root->pRight->tris.size() );
	//getchar();

	if ( pNode->isLeaf ) {

		//return false;

		//printf ( "\n pNode->isLeaf, pNode->tris.size() = %d ", pNode->tris.size() );
		//getchar();

		bool flag = false ; 
		for ( int i=0; i<pNode->tris.size(); ++i ) {

			int & tri1 = pNode->tris[i];

			//printf ( "\n i = %d, tri1 = %d ", i, tri1 );
			//getchar();

			float4 v0 = mesh.vertices [ mesh.indices[ tri1 + 0 ] ];
			float4 e1 = mesh.vertices [ mesh.indices[ tri1 + 1 ] ] - v0;
			float4 e2 = mesh.vertices [ mesh.indices[ tri1 + 2 ] ] - v0;

			float t = spec::RayTriangleIntersection(ray, make_float3(v0.x,v0.y,v0.z), make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));

			if(t < tHit && t > TMIN) {
			
				tHit = t; 

				tri.v0 = v0;
				tri.v1 = e1;
				tri.v2 = e2;

				//if(!needClosestHit)
                //    return;

				//return true; 

				flag = true ; 

				//printf ( "\n flag = true " );
				//getchar();
			}
		}

		return flag ; 

	} else {

		//printf ( "\n bvh.build( mesh1 ) " );
		//printf ( "\n depth = %d, root, tris.size() = %d", 0, pNode->tris.size() );
		//printf ( "\n depth = %d, left, tris.size() = %d", 1, pNode->pLeft->tris.size() );
		//printf ( "\n depth = %d, right, tris.size() = %d", 1, pNode->pRight->tris.size() );
		//getchar();

		BVH_Node* child0 = pNode->pLeft;
        BVH_Node* child1 = pNode->pRight;

		float2 tspan0 = ray_box ( ray, child0->bounds );
		float2 tspan1 = ray_box ( ray, child1->bounds );

		bool intersect0 = (tspan0.x <= tspan0.y ) && ( tspan0.y >= TMIN) && (tspan0.x <= tHit);
		bool intersect1 = (tspan1.x <= tspan1.y ) && ( tspan1.y >= TMIN) && (tspan1.x <= tHit);

		if(intersect0 && intersect1)
			if( tspan0.x > tspan1.x )
			{
				swap1(tspan0,tspan1);
				swap1(child0,child1);
			}

		bool flag = false;

		if(intersect0)
			flag = traceRecursive(mesh, child0,ray, tHit, tri ); // ,needClosestHit);

		//if(result.hit() && !needClosestHit)
		//	return;

//      if(tspan1[TMIN] <= ray.tmax)    // this test helps only about 1-2%
		if(intersect1)
			flag |= traceRecursive(mesh, child1,ray, tHit, tri ); // ,needClosestHit);

		return flag;
	}

	return false ; 
} ; 

bool traceRecursive_cuda_stack ( Mesh & mesh, int node_offset, spec :: Ray & ray, float &tHit, TRI & tri ) { // find_tri

	bool print_flag = true ; // false; // 
	bool print_flag2 = false; // true ; // 

	if ( print_flag2 ) {
		printf( "\n ray->ori = (%.2f, %.2f, %.2f), tHit = %.2f, ray->dir = (%.2f, %.2f, %.2f)", 
			ray.o.x, ray.o.y, ray.o.z, tHit,
			ray.dir.x, ray.dir.y, ray.dir.z  );

		getchar();
	}

	const int stack_size = 64; // 1024; // 

	int stack[stack_size];
	int stack_count = 1;
	stack[0] = 0;

	//if ( bvh_cuda.bvh_nodes.size() < node_offset + 1 ) return false; 

	bool flag = false ; 

	/*
	for(int i = 0; i < mesh1.indices.size()/3; i++) //for(int i = 0; i < total_number_of_triangles; i++)
			{
				//atomicAdd ( &rt_p.pixels_recomputed[2], 1 ) ; 
				//if ( !b_is_recompute && 1 == count && index_tri == i ) continue ; 

				float4 v0 = mesh1.vertices[ mesh1.indices [i*3] ];
				float4 e1 = mesh1.vertices[ mesh1.indices [i*3+1] ] - v0;
				float4 e2 = mesh1.vertices[ mesh1.indices [i*3+2] ] - v0;

				//float4 v0 = triangles [i*3]; // make_float4 ( 1,1,1,1 ) ; // 
				//float4 e1 = triangles [i*3+1]; // make_float4 ( 1,1,1,1 ) ; // 
				//float4 e2 = triangles [i*3+2]; // make_float4 ( 1,1,1,1 ) ; // 

				float t = spec:: RayTriangleIntersection(ray, make_float3(v0.x,v0.y,v0.z),make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));

				if(t < tHit && t > 0.001)
				{
					tHit = t; 
					//hit_r.hit_index = i;

					tri.v0 = v0;
					tri.v1 = e1;
					tri.v2 = e2;

					flag = true ; 

					//break ; 
				}
			}

	//return flag ; 
	*/

	if ( print_flag ) {
		printf( "\n traceRecursive_cuda_stack" );
	}
	
	int flag2 = 0;
	while ( stack_count > 0 ) {

		if ( print_flag ) {
			printf( "\n stack_count = %d , node_index = %d , ", stack_count, stack[stack_count - 1] );
		}

		if ( stack_count > stack_size ) return false;
		if ( flag2++ > 1000 ) return false;

		BVH_Node_ & node = bvh_cuda.bvh_nodes[ stack[stack_count - 1] ]; // node_offset

		if ( node.offset_left < 0 ) { //if ( pNode->isLeaf ) {

			if ( print_flag ) {
				printf( "\n node.offset_left < 0" );
				getchar();
			}
						
			//return false;

			//printf ( "\n pNode->isLeaf, pNode->tris.size() = %d ", pNode->tris.size() );
			//getchar();

			bool flag4 = false ; 
			int tri_index4 = -1;
			for ( int i=0; i<node.num_tris; ++i ) {

				int & tri1 = bvh_cuda.tri_indices[ node.offset_tris + i ]; // pNode->tris[i];

				//printf ( "\n i = %d, tri1 = %d ", i, tri1 );
				//getchar();

				float4 v0 = mesh.vertices [ mesh.indices[ tri1 + 0 ] ];
				float4 e1 = mesh.vertices [ mesh.indices[ tri1 + 1 ] ] - v0;
				float4 e2 = mesh.vertices [ mesh.indices[ tri1 + 2 ] ] - v0;

				float t = spec::RayTriangleIntersection(ray, make_float3(v0.x,v0.y,v0.z), make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));

				if(t < tHit && t > TMIN) {
			
					tHit = t; 

					tri.v0 = v0;
					tri.v1 = e1;
					tri.v2 = e2;

					//if(!needClosestHit)
					//    return;

					//return true; 

					flag4 = true ; 
					tri_index4 = tri1;

					//printf ( "\n flag = true " );
					//getchar();
				}
			}

			if ( print_flag ) {
				if ( flag4 ) {
					printf( "\n tri_index4 = %d ", tri_index4 );
					getchar();
				}
			}

			//return flag ; 
			--stack_count;

		} else {

			//printf ( "\n bvh.build( mesh1 ) " );
			//printf ( "\n depth = %d, root, tris.size() = %d", 0, pNode->tris.size() );
			//printf ( "\n depth = %d, left, tris.size() = %d", 1, pNode->pLeft->tris.size() );
			//printf ( "\n depth = %d, right, tris.size() = %d", 1, pNode->pRight->tris.size() );
			//getchar();
			
			int offset_left = node.offset_left;
			int offset_right = node.offset_right;

			if ( print_flag ) {
				printf( "\n offset_left = %d, offset_right = %d ", offset_left, offset_right );
				getchar();
			}

			BVH_Node_ child0 = bvh_cuda.bvh_nodes [ node.offset_left ];
			BVH_Node_ child1 = bvh_cuda.bvh_nodes [ node.offset_right ];

			float2 tspan0 = ray_box ( ray, child0.aabb );
			float2 tspan1 = ray_box ( ray, child1.aabb );

			bool intersect0 = (tspan0.x <= tspan0.y ) && ( tspan0.y >= TMIN) && (tspan0.x <= tHit);
			bool intersect1 = (tspan1.x <= tspan1.y ) && ( tspan1.y >= TMIN) && (tspan1.x <= tHit);
			
			if(intersect0 && intersect1) {

				if( tspan0.x > tspan1.x ) {
					//swap1(tspan0,tspan1);
					swap1(offset_left,offset_right); //swap1(child0,child1);
				}

				stack[stack_count - 1] = offset_right;

				if ( stack_count >= stack_size ) return false; 
				stack[stack_count] = offset_left;
				++stack_count;

			} else if ( intersect0 ) {
				stack[stack_count - 1] = offset_left;
			} else if ( intersect1 ) {
				stack[stack_count - 1] = offset_right;
			} else {
				--stack_count;
			}
			//--stack_count;
			
			//bool flag = false;

			//if(intersect0)
			//	flag = traceRecursive_cuda(mesh, offset_left,ray, tHit, tri ); // ,needClosestHit);

			//if(result.hit() && !needClosestHit)
			//	return;

	//      if(tspan1[TMIN] <= ray.tmax)    // this test helps only about 1-2%
			//if(intersect1)
			//	flag |= traceRecursive_cuda(mesh, offset_right,ray, tHit, tri ); // ,needClosestHit);

			//return flag;
		}
	}

	//return false ; 
	return flag;
}; 

bool traceRecursive_cuda ( Mesh & mesh, int node_offset, spec :: Ray & ray, float &tHit, TRI & tri ) { // find_tri

	//printf ( "\n pNode->tris.size() = %d ", pNode->tris.size() );
	//printf ( "\n bvh.build( mesh1 ) " );
	//printf ( "\n depth = %d, root, tris.size() = %d", 0, bvh.root->tris.size() );
	//printf ( "\n depth = %d, left, tris.size() = %d", 1, bvh.root->pLeft->tris.size() );
	//printf ( "\n depth = %d, right, tris.size() = %d", 1, bvh.root->pRight->tris.size() );
	//getchar();

	/*
		bool flag = false ; 
		for ( int i=0; i<bvh_cuda.tri_indices.size(); ++i ) {

			int & tri1 = bvh_cuda.tri_indices[ i ]; // pNode->tris[i];

			//printf ( "\n i = %d, tri1 = %d ", i, tri1 );
			//getchar();

			float4 v0 = mesh.vertices [ mesh.indices[ tri1 + 0 ] ];
			float4 e1 = mesh.vertices [ mesh.indices[ tri1 + 1 ] ] - v0;
			float4 e2 = mesh.vertices [ mesh.indices[ tri1 + 2 ] ] - v0;

			float t = spec::RayTriangleIntersection(ray, make_float3(v0.x,v0.y,v0.z), make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));

			if(t < tHit && t > TMIN) {
			
				tHit = t; 

				tri.v0 = v0;
				tri.v1 = e1;
				tri.v2 = e2;

				//if(!needClosestHit)
                //    return;

				//return true; 

				flag = true ; 

				//printf ( "\n flag = true " );
				//getchar();
			}
		}

		return flag ; 
	*/

	if ( bvh_cuda.bvh_nodes.size() < node_offset + 1 ) return false; 

	BVH_Node_ & node = bvh_cuda.bvh_nodes[node_offset];

	if ( node.offset_left < 0 ) { //if ( pNode->isLeaf ) {

		//return false;

		//printf ( "\n pNode->isLeaf, pNode->tris.size() = %d ", pNode->tris.size() );
		//getchar();

		bool flag = false ; 
		for ( int i=0; i<node.num_tris; ++i ) {

			int & tri1 = bvh_cuda.tri_indices[ node.offset_tris + i ]; // pNode->tris[i];

			//printf ( "\n i = %d, tri1 = %d ", i, tri1 );
			//getchar();

			float4 v0 = mesh.vertices [ mesh.indices[ tri1 + 0 ] ];
			float4 e1 = mesh.vertices [ mesh.indices[ tri1 + 1 ] ] - v0;
			float4 e2 = mesh.vertices [ mesh.indices[ tri1 + 2 ] ] - v0;

			float t = spec::RayTriangleIntersection(ray, make_float3(v0.x,v0.y,v0.z), make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));

			if(t < tHit && t > TMIN) {
			
				tHit = t; 

				tri.v0 = v0;
				tri.v1 = e1;
				tri.v2 = e2;

				//if(!needClosestHit)
                //    return;

				//return true; 

				flag = true ; 

				//printf ( "\n flag = true " );
				//getchar();
			}
		}

		return flag ; 

	} else {

		//printf ( "\n bvh.build( mesh1 ) " );
		//printf ( "\n depth = %d, root, tris.size() = %d", 0, pNode->tris.size() );
		//printf ( "\n depth = %d, left, tris.size() = %d", 1, pNode->pLeft->tris.size() );
		//printf ( "\n depth = %d, right, tris.size() = %d", 1, pNode->pRight->tris.size() );
		//getchar();

		int offset_left = node.offset_left;
		int offset_right = node.offset_right;

		BVH_Node_ child0 = bvh_cuda.bvh_nodes [ node.offset_left ];
        BVH_Node_ child1 = bvh_cuda.bvh_nodes [ node.offset_right ];

		float2 tspan0 = ray_box ( ray, child0.aabb );
		float2 tspan1 = ray_box ( ray, child1.aabb );

		bool intersect0 = (tspan0.x <= tspan0.y ) && ( tspan0.y >= TMIN) && (tspan0.x <= tHit);
		bool intersect1 = (tspan1.x <= tspan1.y ) && ( tspan1.y >= TMIN) && (tspan1.x <= tHit);

		if(intersect0 && intersect1)
			if( tspan0.x > tspan1.x )
			{
				swap1(tspan0,tspan1);
				swap1(offset_left,offset_right); //swap1(child0,child1);
			}

		bool flag = false;

		if(intersect0)
			flag = traceRecursive_cuda(mesh, offset_left,ray, tHit, tri ); // ,needClosestHit);

		//if(result.hit() && !needClosestHit)
		//	return;

//      if(tspan1[TMIN] <= ray.tmax)    // this test helps only about 1-2%
		if(intersect1)
			flag |= traceRecursive_cuda(mesh, offset_right,ray, tHit, tri ); // ,needClosestHit);

		return flag;
	}

	return false ; 
} ; 

void raytrace_cpu_pixel_bvh ( int x, int y ) {

	const float xf = (x-0.5)/((float)image_width); // w
	const float yf = (y-0.5)/((float)image_height); // h

	int ray_depth = 0;
	bool continue_path = true;

	float3 t1 = c+(a*xf);
	float3 t2 = b*yf;
	float3 image_pos = t1 + t2;

	spec :: Ray r( image_pos,image_pos - campos );
	//r.print () ; 

	spec:: HitRecord hit_r;

	float t_min,t_max;
	continue_path = spec:: RayBoxIntersection(scene_aabbox_min, scene_aabbox_max, r.o, r.inv_dir,t_min, t_max);
	hit_r.color = make_float3(0,0,0);

	//float3 light_pos =  make_float3(light_x,light_y,light_z) ; 

	// hack to display the light source we simple make a ray sphere intersection and 
	// compare the depth with the found t value from the triangles
	float sphere_t;
	bool sphere_hit = spec:: RaySphereIntersection(r, light_pos,2.0,sphere_t); 

	//float3 light_color1 = make_float3( light_color[0], light_color[1], light_color[2] ) ; 
	
	if(sphere_hit && sphere_t > 0.001)
	{
		if(!continue_path)
		{
			hit_r.color =  light_color1;
		}
		sphere_hit = true;
	}

	int count = 0 ; 

	TRI tri;

	while(continue_path && ray_depth < 1) // 4 
	{
			//bool flag = traceRecursive ( mesh1, bvh.root, r, hit_r.t, tri );
			
			bool flag = traceRecursive_cuda_stack ( mesh1, 0, r, hit_r.t, tri );

			/*for(int i = 0; i < mesh1.indices.size()/3; i++) //for(int i = 0; i < total_number_of_triangles; i++)
			{
				//atomicAdd ( &rt_p.pixels_recomputed[2], 1 ) ; 
				//if ( !b_is_recompute && 1 == count && index_tri == i ) continue ; 

				float4 v0 = mesh1.vertices[ mesh1.indices [i*3] ];
				float4 e1 = mesh1.vertices[ mesh1.indices [i*3+1] ] - v0;
				float4 e2 = mesh1.vertices[ mesh1.indices [i*3+2] ] - v0;

				//float4 v0 = triangles [i*3]; // make_float4 ( 1,1,1,1 ) ; // 
				//float4 e1 = triangles [i*3+1]; // make_float4 ( 1,1,1,1 ) ; // 
				//float4 e2 = triangles [i*3+2]; // make_float4 ( 1,1,1,1 ) ; // 

				float t = spec:: RayTriangleIntersection(r, make_float3(v0.x,v0.y,v0.z),make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));

				if(t < hit_r.t && t > 0.001)
				{
					hit_r.t = t; 
					hit_r.hit_index = i;

					//if ( 1 == count )
					//	pvdata [ y * w + x ] = i ; 
					//pm_triangle_index
					//pm_lvl_triangle_index [ count ] = i ; 

					//find = true ; 

					//break ; 
				}
			}*/

		if(sphere_hit && sphere_t < hit_r.t) // if(sphere_hit && sphere_t < hit_r_t) // 
		{
			hit_r.color += light_color1;
			continue_path = false;
			break;
		}

		if ( flag )
		//if(hit_r.hit_index >= 0)
		{

			ray_depth++;
			
			// create the normal
			//float4 e1 = triangles [ hit_r.hit_index*3+1 ];
			//float4 e2 = triangles [ hit_r.hit_index*3+2 ];

			//float4 v0 = mesh1.vertices[ mesh1.indices [hit_r.hit_index*3] ];
			//float4 e1 = mesh1.vertices[ mesh1.indices [hit_r.hit_index*3+1] ] - v0;
			//float4 e2 = mesh1.vertices[ mesh1.indices [hit_r.hit_index*3+2] ] - v0;

			float4 e1 = tri.v1 ; // - tri.v0;
			float4 e2 = tri.v2 ; // - tri.v0;

			hit_r.normal = cross(make_float3(e1.x,e1.y,e1.z), make_float3(e2.x,e2.y,e2.z));
			//hit_r.normal = cross( triangles [ hit_r.hit_index*3+1 ], triangles [ hit_r.hit_index*3+2 ] ) ; 

			hit_r.normal = normalize(hit_r.normal);

			/*if ( width/2 == x && height/2 == y ) {
				printf ( "\n float4 = (%.2f, %.2f, %.2f )", e1.x, e1.y, e1.z );
				printf ( "\n float4 = (%.2f, %.2f, %.2f )", e2.x, e2.y, e2.z );
				printf ( "\n float4 = (%.2f, %.2f, %.2f )", hit_r.normal.x, hit_r.normal.y, hit_r.normal.z );
			}*/

			// calculate simple diffuse light
			float3 hitpoint = r.o + r.dir *hit_r.t;
			float3 L = light_pos - hitpoint; // rayTraceExtrime
			float dist_to_light = length(L);
			
			L = normalize(L);
			float diffuse_light = fmaxf1( dot(L,hit_r.normal), 0.0);
			diffuse_light = fminf1( (diffuse_light),1.0);

			/*if ( width/2 == x && height/2 == y ) {
				printf ( "\n float = (%.2f)", dot(L, hit_r.normal) );
				printf ( "\n float2 = (%.2f, %.2f)", dist_to_light, diffuse_light );
				getchar();
			}*/

			//calculate simple specular light
			float3 H = L + (-r.dir);
			H = normalize(H);
			float specular_light = powf(fmaxf1(dot(H,hit_r.normal),0.0),25.0f);

			diffuse_light  *=  16.0/dist_to_light;
			specular_light *=  16.0/dist_to_light;

			clamp(diffuse_light, 0.0f, 1.0f);
			clamp(specular_light, 0.0f, 1.0f);

			// rayTraceExtrime
			hit_r.color += light_color1 * diffuse_light + make_float3(1.0,1.0,1.0)*specular_light*0.2 + make_float3(0.2,0.2,0.2); 
		}
		else
		{
			continue_path = false;
			hit_r.color += make_float3(0.5,0.5,0.95*yf+0.3);
		}
	}

	if(ray_depth >= 1 || sphere_hit)
	{
		ray_depth = imax1(ray_depth,1);
		hit_r.color /= ray_depth; // normalize the colors
	}
	else
	{
		hit_r.color = make_float3(0.5,0.5,yf+0.3);
	}

	int val = spec:: rgbToInt(hit_r.color.x*255,hit_r.color.y*255,hit_r.color.z*255);
	out_data[y * image_width + x] = val;
} ; 

void draw_bvh_box( float3 & v0, float3 & v1, float3 & color ) {

	glLineWidth(1);       // ширину линии

	glColor3f (color.x, color.y, color.z);		// текущий цвет примитивов

	glBegin (GL_LINE_LOOP);
		glVertex3f( v0.x, v0.y, v0.z );
		glVertex3f( v0.x, v1.y, v0.z );
		glVertex3f( v1.x, v1.y, v0.z );
		glVertex3f( v1.x, v0.y, v0.z );
	glEnd();

	glBegin (GL_LINE_LOOP);
		glVertex3f( v0.x, v0.y, v1.z );
		glVertex3f( v0.x, v1.y, v1.z );
		glVertex3f( v1.x, v1.y, v1.z );
		glVertex3f( v1.x, v0.y, v1.z );
	glEnd();

	glBegin (GL_LINE_LOOP);
		glVertex3f( v0.x, v0.y, v0.z );
		glVertex3f( v0.x, v1.y, v0.z );
		glVertex3f( v0.x, v1.y, v1.z );
		glVertex3f( v0.x, v0.y, v1.z );
	glEnd();

	glBegin (GL_LINE_LOOP);
		glVertex3f( v1.x, v0.y, v0.z );
		glVertex3f( v1.x, v1.y, v0.z );
		glVertex3f( v1.x, v1.y, v1.z );
		glVertex3f( v1.x, v0.y, v1.z );
	glEnd();

	glBegin (GL_LINE_LOOP);
		glVertex3f( v0.x, v0.y, v0.z );
		glVertex3f( v0.x, v0.y, v1.z );
		glVertex3f( v1.x, v0.y, v1.z );
		glVertex3f( v1.x, v0.y, v0.z );
	glEnd();

	glBegin (GL_LINE_LOOP);
		glVertex3f( v0.x, v1.y, v0.z );
		glVertex3f( v0.x, v1.y, v1.z );
		glVertex3f( v1.x, v1.y, v1.z );
		glVertex3f( v1.x, v1.y, v0.z );
	glEnd();
} ; 

void draw_test0 ();
void set_matrix ();

float3 get_bvh_node_color( int depth ) {
	switch ( depth % 4 ) {
		case 0 : return make_float3(1, 0, 0);
		case 1 : return make_float3(0, 1, 0);
		case 2 : return make_float3(0, 0, 1);
		case 3 : return make_float3(0.5f, 0.5f, 0.5f);
	}
} ; 

void draw_bvh_node ( int depth, BVH_Node * pNode ) {

	if ( NULL == pNode ) return;
	
	//glLineWidth( (5 - depth*2) >= 1 ? (5 - depth*2) : 1 );       // ширину линии
	draw_bvh_box( get_float3(pNode->bounds.min), get_float3(pNode->bounds.max), get_bvh_node_color(depth) );

	if ( pNode->isLeaf ) return;

	draw_bvh_node( ++depth, pNode->pLeft );
	draw_bvh_node( depth, pNode->pRight );
} ; 

void draw_bvh () {

	//draw_test0 ();

	set_matrix ();

	//return;
	//glClear(GL_COLOR_BUFFER_BIT);

	//draw_bvh_box( make_float3(-10, -10, -10), make_float3(10, 10, 10), make_float3(1,0,0) );
	draw_bvh_node( 0, bvh.root );

	glFlush();
    //glutSwapBuffers();
} ; 

void draw_test0 () {

	set_matrix ();

	glClear(GL_COLOR_BUFFER_BIT);

	glLineWidth(1);       // ширину линии
                      // устанавливаем 1
	 glBegin(GL_LINES);
	  glColor3d(1,0,0);     // красный цвет
	  glVertex3d(-4.5,3,0); // первая линия
	  glVertex3d(-3,3,30);
	  glColor3d(0,1,0);     // зеленый
	  glVertex3d(-3,3.3,0); // вторая линия
	  glVertex3d(-40,3.4,0);
	 glEnd();

	glFlush();
    glutSwapBuffers();
} ; 

void draw_line ( float &v0, float &v1 ) {

} ;

void set_matrix () {

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	//gluPerspective(50.0, 1.0, 3.0, 7.0);

	float aspect = (float)image_width / (float)image_height;
	gluPerspective( 60.0, aspect, 1.0, exp(80.0f) );

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt( campos.x, campos.y, campos.z, // 0.0, 0.0, 5.0,
			  0.0, 0.0, 0.0,
			  0.0, 1.0, 0.0);

	glViewport(0, 0, image_width, image_height ) ; // winwidth, winheight);
} ; 

void reshape(int w, int h)
{
	//image_width = ::w = w;
	//image_height = ::h = h;
	
    //glViewport(0, 0, w, h);
       
    //set_matrix ();
}

///
//  Create an OpenCL context on the first available platform using
//  either a GPU or CPU depending on what is available.
//
cl_context CreateContext()
{
    cl_int errNum;
    cl_uint numPlatforms;
    cl_platform_id firstPlatformId;
    cl_context context = NULL;

    // First, select an OpenCL platform to run on.  For this example, we
    // simply choose the first available platform.  Normally, you would
    // query for all available platforms and select the most appropriate one.
    errNum = clGetPlatformIDs(1, &firstPlatformId, &numPlatforms);
    if (errNum != CL_SUCCESS || numPlatforms <= 0)
    {
        std::cerr << "Failed to find any OpenCL platforms." << std::endl;
        return NULL;
    }

    // Next, create an OpenCL context on the platform.  Attempt to
    // create a GPU-based context, and if that fails, try to create
    // a CPU-based context.
    cl_context_properties contextProperties[] =
    {
        CL_CONTEXT_PLATFORM,
        (cl_context_properties)firstPlatformId,
        0
    };
    context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_GPU,
                                      NULL, NULL, &errNum);
    if (errNum != CL_SUCCESS)
    {
        std::cout << "Could not create GPU context, trying CPU..." << std::endl;
        context = clCreateContextFromType(contextProperties, CL_DEVICE_TYPE_CPU,
                                          NULL, NULL, &errNum);
        if (errNum != CL_SUCCESS)
        {
            std::cerr << "Failed to create an OpenCL GPU or CPU context." << std::endl;
            return NULL;
        }
    }

    return context;
}

///
//  Create a command queue on the first device available on the
//  context
//
cl_command_queue CreateCommandQueue(cl_context context, cl_device_id *device)
{
    cl_int errNum;
    cl_device_id *devices;
    cl_command_queue commandQueue = NULL;
    size_t deviceBufferSize = -1;

    // First get the size of the devices buffer
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, 0, NULL, &deviceBufferSize);
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Failed call to clGetContextInfo(...,GL_CONTEXT_DEVICES,...)";
        return NULL;
    }

    if (deviceBufferSize <= 0)
    {
        std::cerr << "No devices available.";
        return NULL;
    }

    // Allocate memory for the devices buffer
    devices = new cl_device_id[deviceBufferSize / sizeof(cl_device_id)];
    errNum = clGetContextInfo(context, CL_CONTEXT_DEVICES, deviceBufferSize, devices, NULL);
    if (errNum != CL_SUCCESS)
    {
        delete [] devices;
        std::cerr << "Failed to get device IDs";
        return NULL;
    }

    // In this example, we just choose the first available device.  In a
    // real program, you would likely use all available devices or choose
    // the highest performance device based on OpenCL device queries
    commandQueue = clCreateCommandQueue(context, devices[0], 0, NULL);
    if (commandQueue == NULL)
    {
        delete [] devices;
        std::cerr << "Failed to create commandQueue for device 0";
        return NULL;
    }

    *device = devices[0];
    delete [] devices;
    return commandQueue;
}

///
//  Create an OpenCL program from the kernel source file
//
cl_program CreateProgram(cl_context context, cl_device_id device, const char* fileName)
{
    cl_int errNum;
    cl_program program;

    std::ifstream kernelFile(fileName, std::ios::in);
    if (!kernelFile.is_open())
    {
        std::cerr << "Failed to open file for reading: " << fileName << std::endl;
        return NULL;
    }

    std::ostringstream oss;
    oss << kernelFile.rdbuf();

    std::string srcStdStr = oss.str();
    const char *srcStr = srcStdStr.c_str();
    program = clCreateProgramWithSource(context, 1,
                                        (const char**)&srcStr,
                                        NULL, NULL);
    if (program == NULL)
    {
        std::cerr << "Failed to create CL program from source." << std::endl;
        return NULL;
    }

    errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (errNum != CL_SUCCESS)
    {
        // Determine the reason for the error
        char buildLog[16384];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                              sizeof(buildLog), buildLog, NULL);

        std::cerr << "Error in kernel: " << std::endl;
        std::cerr << buildLog;
        clReleaseProgram(program);
        return NULL;
    }

    return program;
}

///
//  Attempt to create the program object from a cached binary.  Note that
//  on first run this will fail because the binary has not yet been created.
//
cl_program CreateProgramFromBinary(cl_context context, cl_device_id device, const char* fileName)
{
    FILE *fp = fopen(fileName, "rb");
    if (fp == NULL)
    {
        return NULL;
    }

    // Determine the size of the binary
    size_t binarySize;
    fseek(fp, 0, SEEK_END);
    binarySize = ftell(fp);
    rewind(fp);

    unsigned char *programBinary = new unsigned char[binarySize];
    fread(programBinary, 1, binarySize, fp);
    fclose(fp);

    cl_int errNum = 0;
    cl_program program;
    cl_int binaryStatus;

    program = clCreateProgramWithBinary(context,
                                        1,
                                        &device,
                                        &binarySize,
                                        (const unsigned char**)&programBinary,
                                        &binaryStatus,
                                        &errNum);
    delete [] programBinary;
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Error loading program binary." << std::endl;
        return NULL;
    }

    if (binaryStatus != CL_SUCCESS)
    {
        std::cerr << "Invalid binary for device" << std::endl;
        return NULL;
    }

    errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (errNum != CL_SUCCESS)
    {
        // Determine the reason for the error
        char buildLog[16384];
        clGetProgramBuildInfo(program, device, CL_PROGRAM_BUILD_LOG,
                              sizeof(buildLog), buildLog, NULL);

        std::cerr << "Error in program: " << std::endl;
        std::cerr << buildLog << std::endl;
        clReleaseProgram(program);
        return NULL;
    }

    return program;
}

//
///
//  Retreive program binary for all of the devices attached to the
//  program an and store the one for the device passed in
//
bool SaveProgramBinary(cl_program program, cl_device_id device, const char* fileName)
{
    cl_uint numDevices = 0;
    cl_int errNum;

    // 1 - Query for number of devices attached to program
    errNum = clGetProgramInfo(program, CL_PROGRAM_NUM_DEVICES, sizeof(cl_uint),
                              &numDevices, NULL);
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Error querying for number of devices." << std::endl;
        return false;
    }

    // 2 - Get all of the Device IDs
    cl_device_id *devices = new cl_device_id[numDevices];
    errNum = clGetProgramInfo(program, CL_PROGRAM_DEVICES,
                              sizeof(cl_device_id) * numDevices,
                              devices, NULL);
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Error querying for devices." << std::endl;
        delete [] devices;
        return false;
    }

    // 3 - Determine the size of each program binary
    size_t *programBinarySizes = new size_t [numDevices];
    errNum = clGetProgramInfo(program, CL_PROGRAM_BINARY_SIZES,
                              sizeof(size_t) * numDevices,
                              programBinarySizes, NULL);
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Error querying for program binary sizes." << std::endl;
        delete [] devices;
        delete [] programBinarySizes;
        return false;
    }

    unsigned char **programBinaries = new unsigned char*[numDevices];
    for (cl_uint i = 0; i < numDevices; i++)
    {
        programBinaries[i] = new unsigned char[programBinarySizes[i]];
    }

    // 4 - Get all of the program binaries
    errNum = clGetProgramInfo(program, CL_PROGRAM_BINARIES, sizeof(unsigned char*) * numDevices,
                              programBinaries, NULL);
    if (errNum != CL_SUCCESS)
    {
        std::cerr << "Error querying for program binaries" << std::endl;

        delete [] devices;
        delete [] programBinarySizes;
        for (cl_uint i = 0; i < numDevices; i++)
        {
            delete [] programBinaries[i];
        }
        delete [] programBinaries;
        return false;
    }

    // 5 - Finally store the binaries for the device requested out to disk for future reading.
    for (cl_uint i = 0; i < numDevices; i++)
    {
        // Store the binary just for the device requested.  In a scenario where
        // multiple devices were being used you would save all of the binaries out here.
        if (devices[i] == device)
        {
            FILE *fp = fopen(fileName, "wb");
            fwrite(programBinaries[i], 1, programBinarySizes[i], fp);
            fclose(fp);
            break;
        }
    }

    // Cleanup
    delete [] devices;
    delete [] programBinarySizes;
    for (cl_uint i = 0; i < numDevices; i++)
    {
        delete [] programBinaries[i];
    }
    delete [] programBinaries;
    return true;
}

///
//  Cleanup any created OpenCL resources
//
void Cleanup(cl_context context, cl_command_queue commandQueue,
             cl_program program, cl_kernel kernel ) //, cl_mem memObjects[3])
{
    //for (int i = 0; i < 3; i++)
    //{
    //    if (memObjects[i] != 0)
    //        clReleaseMemObject(memObjects[i]);
    //}
    if (commandQueue != 0)
        clReleaseCommandQueue(commandQueue);

    if (kernel != 0)
        clReleaseKernel(kernel);

    if (program != 0)
        clReleaseProgram(program);

    if (context != 0)
        clReleaseContext(context);

}

int setupCL()
{
	cl_device_id device = 0;

	std::cerr << "\n setupCL()" << std::endl;

	 // Create an OpenCL context on first available platform
    context = CreateContext();
    if (context == NULL)
    {
        std::cerr << "Failed to create OpenCL context." << std::endl;
        return 1;
    }

    // Create a command-queue on the first device available
    // on the created context
    commandQueue = CreateCommandQueue(context, &device);
    if (commandQueue == NULL)
    {
        Cleanup(context, commandQueue, program, kernel ) ; // , memObjects);
        return 1;
    }

    // Create OpenCL program - first attempt to load cached binary.
    //  If that is not available, then create the program from source
    //  and store the binary for future use.
    std::cout << "Attempting to create program from binary..." << std::endl;
    program = CreateProgramFromBinary(context, device, "volumeRender.cl.bin"); // "HelloWorld.cl.bin"
    if (true || program == NULL) // if (program == NULL) // 
    {
        std::cout << "Binary not loaded, create from source..." << std::endl;
        program = CreateProgram(context, device, "volumeRender.cl"); // "HelloWorld.cl"
        if (program == NULL)
        {
            Cleanup(context, commandQueue, program, kernel ) ; // , memObjects);
            return 1;
        }

        std::cout << "Save program binary for future run..." << std::endl;
        if (SaveProgramBinary(program, device, "volumeRender.cl.bin") == false) // "HelloWorld.cl.bin"
        {
            std::cerr << "Failed to write program binary" << std::endl;
            Cleanup(context, commandQueue, program, kernel ) ; // , memObjects);
            return 1;
        }
    }
    else
    {
        std::cout << "Read program from binary." << std::endl;
    }

    // Create OpenCL kernel
    kernel = clCreateKernel(program, "raytracer_bvh", NULL); // "hello_kernel"
    if (kernel == NULL)
    {
        std::cerr << "Failed to create kernel" << std::endl;
        Cleanup(context, commandQueue, program, kernel ) ; // , memObjects);
        return 1;
    }

	std::cerr << "\n setupCL() good" << std::endl;

	return SDK_SUCCESS;
}