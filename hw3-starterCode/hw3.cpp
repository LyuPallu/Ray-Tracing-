/* **************************
 * CSCI 420
 * Assignment 3 Raytracer
 * Name: lyup
 * *************************
*/

#ifdef WIN32
  #include <windows.h>
#endif

#if defined(WIN32) || defined(linux)
  #include <GL/gl.h>
  #include <GL/glut.h>
#elif defined(__APPLE__)
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef WIN32
  #define strcasecmp _stricmp
#endif

#include <imageIO.h>
#include <cmath>
#include <string>
#include <ctime>
#include <iostream>
using namespace std;

//for limitations
#define MAX_TRIANGLES 20000
#define MAX_SPHERES 100
#define MAX_LIGHTS 100

//store the file name
char * filename = NULL;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2

int mode = MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

//for Sphere triangle
#define PI 3.14159265
#define INFINITY 1e8

//
bool antialiasing_mode = false;
int reflection_mode = 0;

///////////////define structures//////
//did nothing
struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};
struct Triangle
{
  Vertex v[3];
};
struct Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
};
struct Light
{
  double position[3];
  double color[3];
};

//// done
//to store rgb
unsigned char pixel_buffer[HEIGHT][WIDTH][3];

struct Color
{
  double r;
  double g;
  double b;

  Color() : 
    r(0), g(0), b(0) {}
  
  Color(double r, double g, double b) : 
    r(r), g(g), b(b) {}

  Color& operator+= (const Color& other) {
		r += other.r;
		if (r > 1.0f) r = 1.0f;
		else if (r < 0.0f) r = 0.0f;
		g += other.g;
		if (g > 1.0f) g = 1.0f;
		else if (g < 0.0f) g = 0.0f;
		b += other.b;
		if (b > 1.0f) b = 1.0f;
		else if (b < 0.0f) b = 0.0f;
		return *this;
	}
  Color operator* (double scalar) const { 
    return Color(scalar * r, scalar * g, scalar * b); }

};

struct Vector
{
  double x;
  double y;
  double z;

  Vector():
    x(0), y(0),z(0) {}
  Vector(double x, double y, double z):
    x(x), y(y), z(z) {} 

  //operator reload
  Vector operator+ (const Vector& vec) const { 
    return Vector(x + vec.x, y + vec.y, z + vec.z); }

  Vector operator- (const Vector& vec) const {
    return Vector(x - vec.x, y - vec.y, z - vec.z);}

  Vector operator* (double scalar) const { 
    return Vector(scalar * x, scalar * y, scalar * z); }

  double dot(const Vector& vec) { 
    return x * vec.x + y * vec.y + z * vec.z; }

	double magnitude() { 
    return sqrt(x * x + y * y + z * z); }

	Vector& cross(const Vector& vec) {
		// using cross product to calculate normal: a x b= (a2b3-a3b2)i+(a3b1-a1b3)j+(a1b2-a2b1)k
		double vx = y * vec.z - z * vec.y;
		double vy = z * vec.x - x * vec.z;
		double vz = x * vec.y - y * vec.x;
		return Vector(vx, vy, vz);
	}

	Vector& normalize() {
		double distance = sqrt(x * x + y * y + z * z);
		if (distance != 0) {
			x = x / distance;
			y = y / distance;
			z = z / distance;
		}
		return *this;
	}

};

////////////////////////////
//ray
struct Ray
{
  Vector origin;
	Vector direction;

  Ray() {}
	Ray(const Vector& origin, const Vector& direction) : origin(origin), direction(direction) {}

  bool intersectTriangle(const Triangle& triangle, Vector& intersection, double& distance);
	bool intersectSphere(const Sphere& sphere, Vector& intersection, double& distance);
	
};

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles = 0;
int num_spheres = 0;
int num_lights = 0;

//for triangle intersection
Vector triangle_intersection(0, 0, INFINITY);
int index_triangle = -1; // record the index of triangle, which is intersected with ray

//for sphere intersection
Vector sphere_intersection(0, 0, INFINITY);
int index_sphere = -1; // record the index of sphere, which is intersected with ray

///////////// call function /////////////
//load scene
int loadScene(char* argv);
void parse_doubles(FILE* file, const char* check, double p[3]);
void parse_rad(FILE* file, double* r);
void parse_shi(FILE* file, double* shi);
void parse_check(const char* expected, char* found);

//for pixels
void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//for triangle
Vector normalAtTriangle(Triangle triangle, Vector intersection);
//Triangle Phong Shading
Color createTriangleShading(Triangle& triangle, Light& light, Vector& intersection);
Color detectTriangleIntersection(Ray& ray, Color& color, double& closest_dist);

//for sphere
Vector normalAtSphere(Sphere sphere, Vector intersection);
//Sphere Phong Shading
Color createSphereShading(Sphere& sphere, Light& light, Vector& intersection);																			  
Color detectSphereIntersection(Ray& ray, Color& color, double& closest_dist);


//reflection ray
Ray createReflectRay(Vector normal, Vector intersection, Vector light_dir) {
	Vector direction = Vector(-light_dir.x, -light_dir.y, -light_dir.z);
	Vector reflect_dir = normal * (direction.dot(normal));
	reflect_dir = (reflect_dir * 2 - direction).normalize();
	Vector reflect_origin = intersection + reflect_dir * 0.001f ; // avoid the same intersection
	return Ray(reflect_origin, reflect_dir);
}

//////////////////////////////////////////////
///////////// triangle part///////////
//normal
Vector normalAtTriangle(Triangle triangle, Vector intersection) {
	//using Barycentric coordinates to calculate normal
	Vector v_1 = Vector(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
	Vector v_2 = Vector(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
	Vector v_3 = Vector(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
	Vector v_inter = intersection;

	Vector v_123_normal = (v_2 - v_1).cross(v_3 - v_1);
	double v_123_area = 0.5f * v_123_normal.magnitude();

	Vector v_23inter_normal = (v_2 - v_inter).cross(v_3 - v_inter);
	double v_23inter_area = 0.5f * v_23inter_normal.magnitude();

	Vector v_13inter_normal = (v_inter - v_1).cross(v_3 - v_1);
	double v_13inter_area = 0.5f * v_13inter_normal.magnitude();

	//calculte normal
	double alpha = v_23inter_area / v_123_area;
	double beta = v_13inter_area / v_123_area;
	double gamma = 1.0 - alpha - beta;
	
	Vector normal = Vector(
		alpha*triangle.v[0].normal[0] + beta*triangle.v[1].normal[0] + gamma*triangle.v[2].normal[0],
		alpha*triangle.v[0].normal[1] + beta*triangle.v[1].normal[1] + gamma*triangle.v[2].normal[1],
		alpha*triangle.v[0].normal[2] + beta*triangle.v[1].normal[2] + gamma*triangle.v[2].normal[2]
	);
	normal.normalize();
	return normal;
}

//triangle intersection
bool Ray::intersectTriangle(const Triangle& triangle, Vector& intersection, double& distance) {
	Vector v_1(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
	Vector v_2(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
	Vector v_3(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);

	Vector normal = (v_2 - v_1).cross(v_3 - v_1);
	normal.normalize();

	double normal_dir = normal.dot(direction); // check whether ray is parallel to the plane
	if (std::abs(normal_dir) <= 1e-15) {	   // ray parallel to plane, no intersection
		return false;
	}
	
	double d = -(v_1.dot(normal));// v_1(x,y,z),normal(a,b,c), ax + by + cz + d = 0;
	double t = -(normal.dot(origin) + d) / normal_dir;
	if (t < 0) {
		return false;
	}
	intersection = origin + direction * t;
	distance = t;
	// in-out test
	// new triangle is composed with intersection, and other 2 nodes from the old triangle
	// if the new triangle's normal and the old triangle's normal don't share the same direction,
	// a.k.a, intersection is not in the old triangle, then it won't be drawn.
	if (normal.dot((v_2 - v_1).cross(intersection - v_1)) < 0 ||
		normal.dot((v_3 - v_2).cross(intersection - v_2)) < 0 ||
		normal.dot((v_1 - v_3).cross(intersection - v_3)) < 0)
	{
		return false;
	}
	return true;
}

//triangle phong shading
Color createTriangleShading(Triangle& triangle, Light& light, Vector& intersection) {
	Color color;
	//Barycentric coordinates for triangles
	Vector v_1 = Vector(triangle.v[0].position[0], triangle.v[0].position[1], triangle.v[0].position[2]);
	Vector v_2 = Vector(triangle.v[1].position[0], triangle.v[1].position[1], triangle.v[1].position[2]);
	Vector v_3 = Vector(triangle.v[2].position[0], triangle.v[2].position[1], triangle.v[2].position[2]);
	Vector v_4 = intersection;

	Vector v_123_normal = (v_2 - v_1).cross(v_3 - v_1);
	double v_123_area = 0.5f * v_123_normal.magnitude();

	Vector v_234_normal = (v_2 - v_4).cross(v_3 - v_4);
	double v_234_area = 0.5f * v_234_normal.magnitude();

	Vector v_134_normal = (v_4 - v_1).cross(v_3 - v_1);
	double v_134_area = 0.5f * v_134_normal.magnitude();

	//calcuate normal 
	double alpha = v_234_area / v_123_area;
	double beta = v_134_area / v_123_area;
	double gamma = 1.0 - alpha - beta;

	Vector normal = Vector(
		alpha*triangle.v[0].normal[0] + beta*triangle.v[1].normal[0] + gamma*triangle.v[2].normal[0],
		alpha*triangle.v[0].normal[1] + beta*triangle.v[1].normal[1] + gamma*triangle.v[2].normal[1],
		alpha*triangle.v[0].normal[2] + beta*triangle.v[1].normal[2] + gamma*triangle.v[2].normal[2]
	);
	normal.normalize();

	//phong shading part
	Vector light_pos = Vector(light.position[0], light.position[1], light.position[2]);
	Vector light_dir = (light_pos - intersection).normalize();

	Vector reflect_dir = normal * (light_dir.dot(normal));
	reflect_dir = (reflect_dir * 2 - light_dir).normalize();

	Vector eye_dir = Vector(-intersection.x, -intersection.y, -intersection.z);
	eye_dir.normalize();

	//final result to 0 and 1.
	double diffuse_factor = light_dir.dot(normal);
	if (diffuse_factor > 1.0f) 
		diffuse_factor = 1.0f;
	else if (diffuse_factor < 0.0f) 
		diffuse_factor = 0.0f;

	double specular_factor = reflect_dir.dot(eye_dir);
	if (specular_factor > 1.0f) 
		specular_factor = 1.0f;
	else if (specular_factor < 0.0f) 
		specular_factor = 0.0f;

	//calcuate kd ks
	Color kd = Color(
		alpha*triangle.v[0].color_diffuse[0] + beta*triangle.v[1].color_diffuse[0] + gamma*triangle.v[2].color_diffuse[0],
		alpha*triangle.v[0].color_diffuse[1] + beta*triangle.v[1].color_diffuse[1] + gamma*triangle.v[2].color_diffuse[1],
		alpha*triangle.v[0].color_diffuse[2] + beta*triangle.v[1].color_diffuse[2] + gamma*triangle.v[2].color_diffuse[2]
	);
	Color ks = Color(
		alpha*triangle.v[0].color_specular[0] + beta*triangle.v[1].color_specular[0] + gamma*triangle.v[2].color_specular[0],
		alpha*triangle.v[0].color_specular[1] + beta*triangle.v[1].color_specular[1] + gamma*triangle.v[2].color_specular[1],
		alpha*triangle.v[0].color_specular[2] + beta*triangle.v[1].color_specular[2] + gamma*triangle.v[2].color_specular[2]
	);
	double sh = alpha*triangle.v[0].shininess + beta*triangle.v[1].shininess + gamma*triangle.v[2].shininess;
	Color lightColor(light.color[0], light.color[1], light.color[2]);

	color.r = lightColor.r * (kd.r * diffuse_factor + ks.r * pow(specular_factor, sh));
	color.g = lightColor.g * (kd.g * diffuse_factor + ks.g * pow(specular_factor, sh));
	color.b = lightColor.b * (kd.b * diffuse_factor + ks.b * pow(specular_factor, sh));

	return color;
}

Color detectTriangleIntersection(Ray& ray, Color& color, double& closest_dist) {
	Color pixel_color = color;
	triangle_intersection = Vector(0, 0, INFINITY);
	// indicate the index of triangle, which is intersected with ray and has minimum distance
	index_triangle = -1;
	double distance = INFINITY;
	for (int i = 0; i < num_triangles; ++i) {
		//when it cannot satisfy the if condition, the intersection has to calculate, and change the right answer.
		Vector tmp_intersect;
		if (ray.intersectTriangle(triangles[i], tmp_intersect, distance)) {
			if (distance < closest_dist && distance > 0){
				closest_dist = distance;
				triangle_intersection = tmp_intersect;
				index_triangle = i;
			}
		}
	}
	// if there is no triangle intersected with the ray, return the original color
	if (index_triangle == -1) {
		return pixel_color;
	}

	pixel_color = Color(0, 0, 0);  // set by default (0,0,0) -- black
	for (int j = 0; j < num_lights; ++j) {
		Vector light_pos(lights[j].position[0], lights[j].position[1], lights[j].position[2]);
		Vector light_dir = (light_pos - triangle_intersection).normalize();
		Ray shadow(triangle_intersection, light_dir); // shadow ray

		 // If the object has light on it, add color to it
		bool hasLight = true;

		//check shadow ray's collisions against every other object
		//triangles
		for (int k = 0; k < num_triangles; ++k) {
			Vector shadow_intersection;
			if (shadow.intersectTriangle(triangles[k], shadow_intersection, distance) && k != index_triangle) {
				Vector dist_shadow = shadow_intersection - triangle_intersection;
				Vector dist_light = light_pos - triangle_intersection;
				// Make sure that the shadow intersection point does not pass the light
				if (dist_shadow.magnitude() < dist_light.magnitude()) {
					hasLight = false;
					break;
				}
			}
		}

		//spheres
		for (int k = 0; k < num_spheres; ++k) {
			Vector shadow_intersection;
			if (shadow.intersectSphere(spheres[k], shadow_intersection, distance)) {
				Vector dist_shadow = shadow_intersection - triangle_intersection;
				Vector dist_light = light_pos - triangle_intersection;
				// Make sure that the shadow intersection point does not pass the light
				if (dist_shadow.magnitude() < dist_light.magnitude()) {
					hasLight = false;
					break;
				}
			}
		}
		
		if (hasLight) {
			pixel_color += createTriangleShading(triangles[index_triangle], lights[j], triangle_intersection);
		}
	}
	return pixel_color;
}


/////////////sphere part///////////////
//normal
//vertex's normal at sphere is equal to its position.
Vector normalAtSphere(Sphere sphere, Vector intersection) {
	Vector normal;
	Vector center(sphere.position[0], sphere.position[1], sphere.position[2]);
	normal = (intersection - center).normalize();
	return normal;
}

// sphere intersection
bool Ray::intersectSphere(const Sphere& sphere, Vector& intersection, double& distance) {
	Vector center(sphere.position[0], sphere.position[1], sphere.position[2]);
	Vector dist = origin - center;
	double r = sphere.radius;

	double dir_ = direction.dot(direction); 
	double dir_dist = 2 * direction.dot(dist);
	double dir_dist_r = dist.dot(dist) - r * r;

	double delta = dir_dist * dir_dist - 4 * dir_ * dir_dist_r;

	double t0 = 0.0f, t1 = 0.0f, t = 0.0f;

	if (delta < 0) {
		return false;
	}

	if (std::abs(delta) <= 1e-15) {
		t0 = -dir_dist / (2 * dir_);
		t1 = t0;
	}
	else {
		//This formula insures that the quantities added for q have the same sign, 
		//avoiding catastropic cancellation.
		//more stable when implemented on computers.
		double q = (dir_dist > 0) ? -0.5f*(dir_dist + sqrt(delta)) : -0.5f*(dir_dist - sqrt(delta));
		t0 = q / dir_;
		t1 = dir_dist_r / q;
	}
	
	if (t0 > t1) { 
		double tmp = t0;
		t0 = t1;
		t1 = tmp;
	}

	if (t0 < 0) {
		t0 = t1;
		if (t0 < 0) {
			return false;
		}
	}
	t = t0;
	intersection = origin + direction * t;
	distance = t;
	return true;
}

//sphere phong shading
Color createSphereShading(Sphere& sphere, Light& light, Vector& intersection) {
	Color color;
	Vector center = Vector(sphere.position[0], sphere.position[1], sphere.position[2]);
	Vector light_pos = Vector(light.position[0], light.position[1], light.position[2]);

	
	Vector normal = (intersection - center).normalize();
	Vector light_dir = (light_pos - intersection).normalize();
	Vector reflect_dir = normal * (light_dir.dot(normal));
	reflect_dir = (reflect_dir * 2 - light_dir).normalize();

	Vector eye_dir = Vector(-intersection.x, -intersection.y, -intersection.z);
	eye_dir.normalize();

	double diffuse_factor = light_dir.dot(normal);
	if (diffuse_factor > 1.0f) 
		diffuse_factor = 1.0f;
	else if (diffuse_factor < 0.0f) 
		diffuse_factor = 0.0f;

	double specular_factor = reflect_dir.dot(eye_dir);
	if (specular_factor > 1.0f) 
		specular_factor = 1.0f;
	else if (specular_factor < 0.0f) 
		specular_factor = 0.0f;

	//calcuate kd ks shininess
	Color lightColor(light.color[0], light.color[1], light.color[2]);
	Color kd(sphere.color_diffuse[0], sphere.color_diffuse[1], sphere.color_diffuse[2]);
	Color ks(sphere.color_specular[0], sphere.color_specular[1], sphere.color_specular[2]);
	double sh = sphere.shininess;

	color.r = lightColor.r * (kd.r * diffuse_factor + ks.r * pow(specular_factor, sh));
	color.g = lightColor.g * (kd.g * diffuse_factor + ks.g * pow(specular_factor, sh));
	color.b = lightColor.b * (kd.b * diffuse_factor + ks.b * pow(specular_factor, sh));

	return color;
}

Color detectSphereIntersection(Ray& ray, Color& color, double& closest_dist) {
	Color pixel_color = color;
	sphere_intersection = Vector(0, 0, INFINITY);
	
	index_sphere = -1;
	double distance = INFINITY;
	for (int i = 0; i < num_spheres; ++i) {
		Vector tmp_intersect;
		if (ray.intersectSphere(spheres[i], tmp_intersect, distance)) {
			if (distance < closest_dist && distance > 0) {
				closest_dist = distance;
				sphere_intersection = tmp_intersect;
				index_sphere = i;
			}
		}
	}
	if (index_sphere == -1) {
		return pixel_color;
	}

	//set black
	pixel_color = Color(0, 0, 0);  

	for (int j = 0; j < num_lights; ++j) {
		Vector light_pos(lights[j].position[0], lights[j].position[1], lights[j].position[2]);
		Vector light_dir = (light_pos - sphere_intersection).normalize();
		Ray shadow(sphere_intersection, light_dir); // shadow ray

		bool hasLight = true;

		//check shadow ray's collisions against every other object
		//spheres
		for (int k = 0; k < num_spheres; ++k) {
			Vector shadow_intersection;
			if (shadow.intersectSphere(spheres[k], shadow_intersection, distance) && k != index_sphere) {
				Vector dist_shadow = shadow_intersection - sphere_intersection;
				Vector dist_light = light_pos - sphere_intersection;

				// justify whether the shadow ray could reach the light
				if (dist_shadow.magnitude() < dist_light.magnitude()) { // if cannot, set black
					hasLight = false;
					break;
				}
			}
		}
		//triangles
		for (int k = 0; k < num_triangles; ++k) {
			Vector shadow_intersection;
			if (shadow.intersectTriangle(triangles[k], shadow_intersection, distance)) {
				Vector dist_shadow = shadow_intersection - sphere_intersection;
				Vector dist_light = light_pos - sphere_intersection;

				// Make sure that the shadow intersection point does not pass the light
				if (dist_shadow.magnitude() < dist_light.magnitude()) {
					hasLight = false;
					break;
				}
			}
		}
		if (hasLight) {
			pixel_color += createSphereShading(spheres[index_sphere], lights[j], sphere_intersection);
		}
	}
	return pixel_color;
}
////////////////
Color createColor(Ray& ray, int count_reflect) {
	if (count_reflect > reflection_mode) {
		return Color(0, 0, 0);
	}
	Color color;
	Color lit_color, reflect_color;
	double closest_dist = INFINITY;

	lit_color = detectSphereIntersection(ray, lit_color, closest_dist);
	lit_color = detectTriangleIntersection(ray, lit_color, closest_dist);

	//reflection part
	//sphere
	if (index_sphere != -1) {
		Vector normal_reflect = normalAtSphere(spheres[index_sphere], sphere_intersection);
		Ray reflectRay = createReflectRay(normal_reflect, sphere_intersection, ray.direction);
		double reflect_dist;
		reflect_color = createColor(reflectRay, count_reflect + 1);
	}
	//triangle
	if (index_triangle != -1) {
		Vector normal_reflect = normalAtTriangle(triangles[index_triangle], triangle_intersection);
		Ray reflectRay = createReflectRay(normal_reflect, triangle_intersection, ray.direction);
		double reflect_dist;
		reflect_color = createColor(reflectRay, count_reflect + 1);
	}
	
	color += lit_color;
	color += reflect_color;
	return color;
}
Ray createRay(double x, double y) {
	double ratio = (double)WIDTH / double(HEIGHT);
	double fov_tan = tan((fov / 2) * (PI / 180.0));

	//screen coordinates, -1 <= x, y <= 1
	double screen_x = 2 * (x / (double)WIDTH) - 1;
	double screen_y = 2 * (y / (double)HEIGHT) - 1;

	//position on frontal view
	double actual_x = screen_x*ratio*fov_tan;
	double actual_y = screen_y*fov_tan;
	double actual_z = -1;

	Vector origin; // camera is at (0,0,0)
	Vector direction(actual_x, actual_y, actual_z);
	direction.normalize();

	return Ray(origin, direction);
}
////////////////////////////////////////

//MODIFY THIS FUNCTION
void draw_scene()
{
  //a simple test output
  for(unsigned int x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(unsigned int y=0; y<HEIGHT; y++)
    {
	  //antialiasing
      if (antialiasing_mode) { 
				//average colors of 4 adjacent rays to calculate the final color
				Ray rays[4];
				Color colors[4];
				double r = 0.0, g = 0.0, b = 0.0;
				rays[0] = createRay(x - 0.25, y - 0.25);
				colors[0] = createColor(rays[0], 0);
				rays[1] = createRay(x - 0.25, y + 0.25);
				colors[1] = createColor(rays[1], 0);
				rays[2] = createRay(x + 0.25, y + 0.25);
				colors[2] = createColor(rays[2], 0);
				rays[3] = createRay(x + 0.25, y - 0.25);
				colors[3] = createColor(rays[3], 0);
				for (int i = 0; i < 4; ++i) {
					r += colors[i].r;
					g += colors[i].g;
					b += colors[i].b;
				}
				r = r / 4;
				g = g / 4;
				b = b / 4;
				plot_pixel(x, y, r * 255, g * 255, b * 255);
			}
			// without antialiasing
			else {
				Ray ray = createRay(x, y);
				Color color = createColor(ray, 0);
				plot_pixel(x, y, color.r * 255, color.g * 255, color.b * 255);
			}
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); 
  fflush(stdout);
}
////////////////////////////
//did nothing
void plot_pixel_display(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  glColor3f(((float)r) / 255.0f, ((float)g) / 255.0f, ((float)b) / 255.0f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  pixel_buffer[y][x][0] = r;
  pixel_buffer[y][x][1] = g;
  pixel_buffer[y][x][2] = b;
}

void plot_pixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
    plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  printf("Saving JPEG file: %s\n", filename);

  ImageIO img(WIDTH, HEIGHT, 3, &pixel_buffer[0][0][0]);
  if (img.save(filename, ImageIO::FORMAT_JPEG) != ImageIO::OK)
    printf("Error in Saving\n");
  else 
    printf("File saved Successfully\n");
}

void parse_check(const char *expected, char *found)
{
  if(strcasecmp(expected,found))
  {
    printf("Expected '%s ' found '%s '\n", expected, found);
    printf("Parse error, abnormal abortion\n");
    exit(0);
  }
}

void parse_doubles(FILE* file, const char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE *file, double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE *file, double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE * file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i", &number_of_objects);

  printf("number of objects: %i\n",number_of_objects);

  parse_doubles(file,"amb:",ambient_light);

  for(int i=0; i<number_of_objects; i++)
  {
    fscanf(file,"%s\n",type);
    printf("%s\n",type);
    if(strcasecmp(type,"triangle")==0)
    {
      printf("found triangle\n");
      for(int j=0;j < 3;j++)
      {
        parse_doubles(file,"pos:",t.v[j].position);
        parse_doubles(file,"nor:",t.v[j].normal);
        parse_doubles(file,"dif:",t.v[j].color_diffuse);
        parse_doubles(file,"spe:",t.v[j].color_specular);
        parse_shi(file,&t.v[j].shininess);
      }

      if(num_triangles == MAX_TRIANGLES)
      {
        printf("too many triangles, you should increase MAX_TRIANGLES!\n");
        exit(0);
      }
      triangles[num_triangles++] = t;
    }
    else if(strcasecmp(type,"sphere")==0)
    {
      printf("found sphere\n");

      parse_doubles(file,"pos:",s.position);
      parse_rad(file,&s.radius);
      parse_doubles(file,"dif:",s.color_diffuse);
      parse_doubles(file,"spe:",s.color_specular);
      parse_shi(file,&s.shininess);

      if(num_spheres == MAX_SPHERES)
      {
        printf("too many spheres, you should increase MAX_SPHERES!\n");
        exit(0);
      }
      spheres[num_spheres++] = s;
    }
    else if(strcasecmp(type,"light")==0)
    {
      printf("found light\n");
      parse_doubles(file,"pos:",l.position);
      parse_doubles(file,"col:",l.color);

      if(num_lights == MAX_LIGHTS)
      {
        printf("too many lights, you should increase MAX_LIGHTS!\n");
        exit(0);
      }
      lights[num_lights++] = l;
    }
    else
    {
      printf("unknown type in scene description:\n%s\n",type);
      exit(0);
    }
  }
  return 0;
}

void display()
{
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  clock_t start;
	double duration;
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
    start = clock();
    draw_scene();
    if(mode == MODE_JPEG)
      save_jpg();
    duration = (std::clock() - start) / (double)CLOCKS_PER_SEC;
		std::cout << "Time: " << duration << '\n';
  }
  once=1;
}

int main(int argc, char ** argv)
{
  if ((argc < 2) || (argc > 5))
  {  
    printf ("Usage: %s <input scenefile> [output jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
  {
    mode = MODE_JPEG;
    filename = argv[2];
  }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  ///////////////////
  if (argc == 5) {
		std::string str1(argv[3]);
		if (str1 == "antialiasing") {
			printf("antialiasing");
			antialiasing_mode = true;
		}
		else {
			printf("without antialiasing\n");
			antialiasing_mode = false;
		}
		std::string str2(argv[4]);
		if (str2 == "reflection") {
			reflection_mode = 3;
		}
		else {
			printf("without reflection\n");
			reflection_mode = 0;
		}
		mode = MODE_JPEG;
		filename = argv[2];
	}
	if (argc == 4) {
		std::string str(argv[3]);
		if (str == "antialiasing") {
			antialiasing_mode = true;
		}
		else if (str == "reflection") {
			reflection_mode = 3;
		}
		else {
			printf("without antialiasing\n");
			antialiasing_mode = false;
			reflection_mode = 0;
		}
		mode = MODE_JPEG;
		filename = argv[2];
	}

  /////////
  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  #ifdef __APPLE__
    // This is needed on recent Mac OS X versions to correctly display the window.
    glutReshapeWindow(WIDTH - 1, HEIGHT - 1);
  #endif
  glutDisplayFunc(display);
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}

