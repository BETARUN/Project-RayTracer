#pragma once
#ifndef TRIANGLEH
#define TRIANGLEH

#include "hitable.h"

class triangle : public hitable {
public:
	triangle() {}
	triangle(vec3 p0, vec3 p1, vec3 p2, material* m) : point0(p0), point1(p1), point2(p2), mat_ptr(m) {
		normal = cross((p1 - p0), (p2 - p0));
		normal.make_unit_vector();
	};
	virtual bool hit(const ray& ray, float t_min, float t_max, hit_record& rec) const;
	virtual bool bounding_box(const float &t0, const float &t1, aabb &box) const;
	
	vec3 normal;
	vec3 point0, point1, point2;
	material* mat_ptr;
};

bool triangle::hit(const ray& ray, float t_min, float t_max, hit_record& ret) const {
	float n_dot_dir = dot(normal, ray.direction());
	if (equal(n_dot_dir, 0.0))
		return false;
	float d = -dot(normal, point0);
	float t = -(dot(normal, ray.origin()) + d) / n_dot_dir;
	if (t < t_min || t > t_max)
		return false;
	ret.t = t;
	ret.p = ray.point_at_parameter(t);
	ret.mat_ptr = mat_ptr;
	vec3 r = ret.p - point0;
	vec3 q1 = point1 - point0;
	vec3 q2 = point2 - point0;
	float q1_squaredLen = q1.length() * q1.length();
	float q2_squaredLen = q2.length() * q2.length();
	float q1_dot_q2 = dot(q1, q2);
	float r_dot_q1 = dot(r, q1);
	float r_dot_q2 = dot(r, q2);
	float determinant = 1.0f / (q1_squaredLen * q2_squaredLen - q1_dot_q2 * q1_dot_q2);

	float omega1 = determinant * (q2_squaredLen * r_dot_q1 - q1_dot_q2 * r_dot_q2);
	float omega2 = determinant * (-q1_dot_q2 * r_dot_q1 + q1_squaredLen * r_dot_q2);
	if (omega1 + omega2 > 1.0f || omega1 < 0.0f || omega2 < 0.0f)
		return false;
	return true;
}

bool triangle::bounding_box(const float &t0, const float &t1, aabb &box) const {
	vec3 minp, maxp;
	minp.x = min(point0.x, min(point1.x, point2.x));
	minp.y = min(point0.y, min(point1.y, point2.y));
	minp.z = min(point0.z, min(point1.z, point2.z));
	maxp.x = max(point0.x, max(point1.x, point2.x));
	maxp.y = max(point0.y, max(point1.y, point2.y));
	maxp.z = max(point0.z, max(point1.z, point2.z));
	box = aabb(minp, maxp);
	return true;
}

#endif
