#include "./code/PathTracer.h"

#include <time.h>
#include <iostream>
#include <cstdlib>

#include "./code/sphere.h"
#include "./code/hitable_list.h"
#include "./code/camera.h"
#include "./code/material.h"
#include "./code/texture.h"
#include "./code/triangle.h"
#include "./code/model.h"
#include "./code/aabb.h"
#include "./code/bvh.h"
#include "./code/skybox.h"

#include "../include/tbb/parallel_for.h"

//#include "model.h"
#define STB_IMAGE_IMPLEMENTATION
#include "./code/stb_image.h"

using namespace tbb;
using namespace std;

hitable_list *amodel() {
	hitable **list = new hitable*[10];
	int i = 0;

	//int nx, ny, nn;
	//unsigned char *tex_data = stbi_load("./resource/picture/robot.jpg", &nx, &ny, &nn, 0);
	model *board = new model("./resource/model/board.obj", vec3(3, 3, 0), vec3(1, 1, 1), new lambertian(new constant_texture(vec3(1, 1, 1))));
	board->rotate(vec3(0, 1, 0), 90.0f);
	board->scale(vec3(3, 3, 3));

	list[i++] = board;
	list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new light(new constant_texture(vec3(1, 1, 1))));
	list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));

	return new hitable_list(list, i);
}

hitable_list *atriangle() {
	hitable **list = new hitable*[1];

	//list[0] = new sphere(vec3(0, 0, 0), 2, mat);
	int nx, ny, nn;
	//unsigned char *tex_data = stbi_load("tiled.jpg", &nx, &ny, &nn, 0);
	unsigned char *tex_data = stbi_load("./resource/picture/earthmap.jpg", &nx, &ny, &nn, 0);
	material *mat = new lambertian(new image_texture(tex_data, nx, ny));

	list[0] = new triangle(vec3(1, 1, 1), vec3(0, 1, 0), vec3(2, 2, 2), new lambertian(new constant_texture(vec3(0.4, 0.2, 0.1))));
	return new hitable_list(list, 1);
}

hitable_list *earth() {
	hitable **list = new hitable*[1];

	int nx, ny, nn;
	//unsigned char *tex_data = stbi_load("tiled.jpg", &nx, &ny, &nn, 0);
	unsigned char *tex_data = stbi_load("./resource/picture/earthmap.jpg", &nx, &ny, &nn, 0);
	material *mat = new lambertian(new image_texture(tex_data, nx, ny));

	list[0] = new sphere(vec3(0, 0, 0), 2, mat);
	return new hitable_list(list, 1); 
}
hitable_list *two_perlin_spheres() {
	texture *pertext = new noise_texture(4);
	hitable **list = new hitable*[2];
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(pertext));
	list[1] = new triangle(vec3(1, 1, 1), vec3(0, 1, 0), vec3(2, 2, 2), new light(new constant_texture(vec3(0.4, 0.2, 0.1))));
	return new hitable_list(list, 2);
}

hitable_list *random_scene() {
	int n = 500;
	hitable **list = new hitable*[n + 1];
	
	list[0] = new sphere(vec3(0, -1000, 0), 1000, new lambertian(new constant_texture(vec3(1, 1, 1))));
	int i = 1;
	for (int a = -11; a < 11; a++) {
		for (int b = -11; b < 11; b++) {
			float choose_mat = drand48();
			vec3 center(a + 0.9*drand48(), 0.2, b + 0.9*drand48());
			if ((center - vec3(4, 0.2, 0)).length() > 0.9) {
				if (choose_mat < 0.8) {
					list[i++] = new sphere(center, 0.2, new light(new constant_texture(vec3(drand48()*drand48(), drand48()*drand48(), drand48()*drand48()))));
				}
				else if (choose_mat < 0.95) {
					list[i++] = new sphere(center, 0.2, new metal(vec3(0.5*(1 + drand48()), 0.5*(1 + drand48()), 0.5*(1+ drand48())), 0.5*drand48()));
				}
				else {
					list[i++] = new sphere(center, 0.2, new dielectric(1.5));
				}
			}
		}
	}

	list[i++] = new sphere(vec3(0, 1, 0), 1.0, new dielectric(1.5));
	list[i++] = new sphere(vec3(-4, 1, 0), 1.0, new light(new constant_texture(vec3(1, 1, 1))));
	list[i++] = new sphere(vec3(4, 1, 0), 1.0, new dielectric(1.5));

	model *board = new model("./resource/model/board.obj", vec3(0, 1, 0), vec3(1, 1, 1), new lambertian(new constant_texture(vec3(0.4, 0.2, 0.1))));
	board->scale(vec3(0.5, 0.5, 0.5));
	board->translate(vec3(0, 2, 0));
	list[i++] = board;
	
	
	return new hitable_list(list, i);
}

hitable_list *final_scene() {
	hitable **list = new hitable*[20];
	int i = 0;

	list[i++] = new sphere(vec3(0, -1000, 0), 1000, new metal(vec3(0.7, 0.6, 0.5), 0.0));

	int nx, ny, nn;
	unsigned char *tex_data1 = stbi_load("./resource/picture/earthmap.jpg", &nx, &ny, &nn, 0);
	list[i++] = new sphere(vec3(1.3, 1, -6), 1.0, new lambertian(new image_texture(tex_data1, nx, ny)));

	texture *perlintex = new noise_texture(4);
	list[i++] = new sphere(vec3(4, 2, 0), 2, new lambertian(perlintex));

	list[i++] = new sphere(vec3(3.5, 0.8, -7), 0.8, new light(new constant_texture(vec3(0.8, 0.7, 0.2))));
	
	unsigned char *tex_data2 = stbi_load("./resource/picture/robot.jpg", &nx, &ny, &nn, 0);
	model *robot = new model("./resource/model/board.obj", vec3(-3, 2.3, -7), vec3(1, 1, 1), new lambertian(new constant_texture(vec3(0.9, 0.9, 0.9))));
	robot->rotate(vec3(0, 1, 0), 90.0f);
	robot->scale(vec3(5, 5, 5));
	list[i++] = robot;
	
	return new hitable_list(list, i);

}

vec3 PathTracer::color(const ray& r, hitable *world, int depth) {
	hit_record rec;
	if (world->hit(r, 0.001, FLT_MAX, rec)) {
		ray scattered;
		vec3 attenuation;
		vec3 emitted = rec.mat_ptr->emitted(rec.texcoord.x, rec.texcoord.y, rec.p);
		if (depth < 50 && rec.mat_ptr->scatter(r, rec, attenuation, scattered)) {
			return emitted + attenuation * color(scattered, world, depth + 1);
		}
		else {
			return emitted;
		}
	}
	else {
		//vec3 unit_direction = unit_vector(r.direction());
		//float t = 0.5*(unit_direction.y + 1.0);
		//return (1.0 - t) * vec3(1.0, 1.0, 1.0) + t * vec3(0.5, 0.7, 1.0);
		vec3 tr = skybox_str->sampleBackground(r);
		//vec3 tr(0, 0, 0);
		return tr;
	}
}

PathTracer::PathTracer()
	: m_channel(4), m_width(800), m_height(600), m_image(nullptr) {}

PathTracer::~PathTracer()
{
	if (m_image != nullptr)
		m_image;
	m_image = nullptr;
}

void PathTracer::initialize(int width, int height)
{
	m_width = width;
	m_height = height;
	if (m_image != nullptr)
		delete m_image;

	// allocate pixel buffer, RGBA format.
	m_image = new unsigned char[width * height * m_channel];
}

unsigned char * PathTracer::render(double & timeConsuming)
{
	if (m_image == nullptr)
	{
		std::cout << "Must call initialize() before rendering.\n";
		return nullptr;
	}

	// record start time.
	double startFrame = clock();

	const string path = "./resource/skybox";

	int nx, ny, nn;
	unsigned char *tex_data = stbi_load("./resource/skybox/right.jpg", &nx, &ny, &nn, 0);
	texture* right = new image_texture(tex_data, nx, ny);
	tex_data = stbi_load("./resource/skybox/left.jpg", &nx, &ny, &nn, 0);
	texture* left = new image_texture(tex_data, nx, ny);
	tex_data = stbi_load("./resource/skybox/top.jpg", &nx, &ny, &nn, 0);
	texture* top = new image_texture(tex_data, nx, ny);
	tex_data = stbi_load("./resource/skybox/bottom.jpg", &nx, &ny, &nn, 0);
	texture* bottom = new image_texture(tex_data, nx, ny);
	tex_data = stbi_load("./resource/skybox/back.jpg", &nx, &ny, &nn, 0);
	texture* back = new image_texture(tex_data, nx, ny);
	tex_data = stbi_load("./resource/skybox/front.jpg", &nx, &ny, &nn, 0);
	texture* front = new image_texture(tex_data, nx, ny);
	
	skybox_str = new skybox({ front, back, left, right, top, bottom });
	
	int ns = 10;
	//hitable_list *objects = amodel();
	//hitable_list *objects = atriangle();
	//hitable_list *objects = earth();
	//hitable_list *objects = random_scene();
	//hitable_list *objects = two_perlin_spheres();
	
	hitable_list *objects = final_scene();
	
	bvh_node *root = nullptr;

	for (int x = 0; x < objects->list_size; ++x)
	{
		objects->list[x]->preRendering();
	}
	if (root) delete root;
	root = new bvh_node(objects->list, objects->list_size, 0.0f, FLT_MAX);

	hitable* hitableNode = reinterpret_cast<hitable*>(root);
	

	vec3 lookfrom(0, 4, -30);
	vec3 lookat(0, 3, 0);
	vec3 vup(0, 1, 0);
	float dist_to_focus = 10.0;
	float aperture = 0.0;
	camera cam(lookfrom, lookat, vup, 20, static_cast<float>(m_width) / static_cast<float>(m_height), aperture, dist_to_focus);

	// render the image pixel by pixel.
	parallel_for(blocked_range<size_t>(0, m_height * m_width, 10000), [&](blocked_range<size_t>& r) {
		for (size_t i = r.begin(); i != r.end(); ++i) {
			vec3 ray_color(0, 0, 0);
			size_t col = i % m_width;
			size_t row = i / m_width;
			for (int s = 0; s < ns; s++) {
				float u = static_cast<float>(col + drand48()) / static_cast<float>(m_width);
				float v = static_cast<float>(row + drand48()) / static_cast<float>(m_height);
				//float u = static_cast<float>(col) / static_cast<float>(m_width);
				//float v = static_cast<float>(row) / static_cast<float>(m_height);
				ray r = cam.get_ray(u, v);
				ray_color += color(r, hitableNode, 0);
			}
			ray_color /= static_cast<float>(ns);
			ray_color = vec3(sqrt(ray_color[0]), sqrt(ray_color[1]), sqrt(ray_color[2]));
			drawPixel(col, row, ray_color);
		}
	}, auto_partitioner());

	// record end time.
	double endFrame = clock();

	// calculate time consuming.
	timeConsuming = static_cast<double>(endFrame - startFrame) / CLOCKS_PER_SEC;

	return m_image;
}

void PathTracer::drawPixel(unsigned int x, unsigned int y, const vec3 & color)
{
	// Check out 
	if (x < 0 || x >= m_width || y < 0 || y >= m_height)
		return;
	// x is column's index, y is row's index.
	unsigned int index = (y * m_width + x) * m_channel;
	// store the pixel.
	// red component.
	m_image[index + 0] = static_cast<unsigned char>(255 * color.x);
	// green component.
	m_image[index + 1] = static_cast<unsigned char>(255 * color.y);
	// blue component.
	m_image[index + 2] = static_cast<unsigned char>(255 * color.z);
	// alpha component.
	m_image[index + 3] = static_cast<unsigned char>(255);
}
