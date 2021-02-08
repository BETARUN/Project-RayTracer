#pragma once
#ifndef HITABLEH
#define HITABLEH

#include "ray.h"
#include "aabb.h"

class material;

void get_sphere_uv(const vec3& p, vec2 &tex) {
	float phi = atan2(p.z, p.x);
	float theta = asin(p.y);
	tex = vec2(1 - (phi + M_PI) / (2 * M_PI), (theta + M_PI / 2) / M_PI);
}

struct hit_record {
	float t;
	vec3 p;
	vec3 normal;
	material* mat_ptr;
	vec2 texcoord;
};

struct vertex {
	vec3 p;
	vec3 normal;
	vec2 texcoord;
};

class hitable {
public:
	virtual bool hit(const ray& ray,const float t_min,const float t_max, hit_record& rec) const = 0;
	virtual bool bounding_box(const float &t0, const float &t1, aabb &box) const = 0;
	virtual void preRendering() {}
};

#endif
