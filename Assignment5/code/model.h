#pragma once

#pragma once
#ifndef MODELH
#define MODELH

//#include "mathutils.h"
#include "hitable.h"
#include "transform.h"
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

class model : public hitable {
public:

	model(const std::string &path, vec3 pos, vec3 len, material* mat);

	bool triangleHit(const ray &r, const float &t_min, const float &t_max, hit_record &ret, const vertex &p0, const vertex &p1, const vertex &p2, const vec3 &normal) const;
	virtual bool hit(const ray &ray, float t_min, float t_max, hit_record &rec) const;
	virtual bool bounding_box(const float &t0, const float &t1, aabb &box) const;
	virtual void preRendering();

	void loadObjFile(const std::string &path);

	void scale(const vec3 &ds) { transformation.scale(ds); }
	void translate(const vec3 &dt) { transformation.translate(dt); }
	void rotate(const vec3 &axis, float angle) { transformation.rotate(axis, angle); }

	m_transform transformation;
	material* mat_ptr;
	std::vector<vertex> vertices;
	std::vector<unsigned int> indices;
	std::vector<vec3> face_normal;
	aabb m_box;

	vec3 model_center;
	vec3 model_scale;
};

model::model(const std::string &path, vec3 pos, vec3 len, material* mat)
{
	mat_ptr = mat;
	loadObjFile(path);
	translate(pos - model_center);
	float length = max(max(model_scale.x, model_scale.y), model_scale.z);
	scale(/*len*/ vec3(1.0f / length, 1.0f / length, 1.0f / length));
}

bool model::hit(const ray &ray, float t_min, float t_max, hit_record &rec)const {

	hit_record tmpRec;
	bool hitSome = false;
	float closestSoFar = t_max;
	for (int x = 0; x < indices.size(); x += 3)
	{
		int index1 = indices[x + 0];
		int index2 = indices[x + 1];
		int index3 = indices[x + 2];
		if (triangleHit(ray, t_min, closestSoFar, tmpRec,
			vertices[index1],
			vertices[index2],
			vertices[index3],
			face_normal[x / 3]))
		{
			hitSome = true;
			closestSoFar = tmpRec.t;
			rec = tmpRec;
		}
	}
	return hitSome;
}

bool model::triangleHit(const ray &ray, const float &t_min, const float &t_max, hit_record &ret, const vertex &p0, const vertex &p1, const vertex &p2, const vec3 &normal) const {

	float n_dot_dir = dot(normal, ray.direction());
	if (equal(n_dot_dir, 0.0))
		return false;
	float d = -dot(normal, p0.p);
	float t = -(dot(normal, ray.origin()) + d) / n_dot_dir;
	if (t < t_min || t > t_max)
		return false;
	ret.t = t;
	ret.p = ray.point_at_parameter(t);
	ret.mat_ptr = mat_ptr;
	vec3 r = ret.p - p0.p;
	vec3 q1 = p1.p - p0.p;
	vec3 q2 = p2.p - p0.p;
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
	ret.normal = p0.normal * (1.0f - omega1 - omega2) + p1.normal * omega1 + p2.normal * omega2;
	ret.texcoord = p0.texcoord * (1.0f - omega1 - omega2) + p1.texcoord * omega1 + p2.texcoord * omega2;
	if (dot(ret.normal, ray.direction()) > 0.0f)
		ret.normal = -ret.normal;
	return true;
}

bool model::bounding_box(const float &t0, const float &t1, aabb &box) const {
	(void)t0;
	(void)t1;
	box = m_box;
	return true;
}

void model::loadObjFile(const std::string &path)
{

	// obj loader.
	ifstream in;
	in.open(path, ifstream::in);
	if (in.fail())
	{
		std::cout << "Fail to load obj->" << path << endl;
	}
	string line;
	vector<vec3> model_vertices;
	vector<vec3> model_normals;
	vector<vec2> model_texcoords;
	vec3 minPoint(+FLT_MAX, +FLT_MAX, +FLT_MAX);
	vec3 maxPoint(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	while (!in.eof())
	{
		getline(in, line);
		istringstream iss(line.c_str());
		char trash;
		//vertex
		if (!line.compare(0, 2, "v "))
		{
			iss >> trash;
			vec3 model_vertex;
			iss >> model_vertex.x;
			iss >> model_vertex.y;
			iss >> model_vertex.z;
			model_vertices.push_back(model_vertex);

			if (minPoint.x > model_vertex.x)minPoint.x = model_vertex.x;
			if (minPoint.y > model_vertex.y)minPoint.y = model_vertex.y;
			if (minPoint.z > model_vertex.z)minPoint.z = model_vertex.z;
			if (maxPoint.x < model_vertex.x)maxPoint.x = model_vertex.x;
			if (maxPoint.y < model_vertex.y)maxPoint.y = model_vertex.y;
			if (maxPoint.z < model_vertex.z)maxPoint.z = model_vertex.z;
		}
		// normal
		else if (!line.compare(0, 3, "vn "))
		{
			iss >> trash >> trash;
			vec3 model_normal;
			iss >> model_normal.x;
			iss >> model_normal.y;
			iss >> model_normal.z;
			model_normal.make_unit_vector();
			model_normals.push_back(model_normal);
		}
		// texcoord
		else if (!line.compare(0, 3, "vt "))
		{
			iss >> trash >> trash;
			vec2 model_texcoord;
			iss >> model_texcoord.x;
			iss >> model_texcoord.y;
			model_texcoords.push_back(model_texcoord);
		}
		// face
		else if (!line.compare(0, 2, "f "))
		{
			iss >> trash;
			int index[3];
			while (iss >> index[0] >> trash >> index[1] >> trash >> index[2])
			{
				vertex data;
				data.p = model_vertices[index[0] - 1];
				data.texcoord = model_texcoords[index[1] - 1];
				data.normal = model_normals[index[2] - 1];
				indices.push_back(vertices.size());
				vertices.push_back(data);
			}
		}
	}
	in.close();
	model_center = (maxPoint + minPoint) / 2;
	model_scale.x = maxPoint.x - minPoint.x;
	model_scale.y = maxPoint.y - minPoint.y;
	model_scale.z = maxPoint.z - minPoint.z;

}

inline vec3& vec3::operator=(const vec4 &v) {
	e[0] = v.e[0];
	e[1] = v.e[1];
	e[2] = v.e[2];
	return *this;
}

void model::preRendering() {
	// transform and calculate aabb box.
	if (!transformation.getDirtry() || indices.empty())
		return;
	vec3 minPoint(+FLT_MAX, +FLT_MAX, +FLT_MAX);
	vec3 maxPoint(-FLT_MAX, -FLT_MAX, -FLT_MAX);
	matrix4 modelMatrix = transformation.toMatrix();
	matrix4 invModelMatrix = transformation.toInvMatrix();
	vec4 pos, nor;
	unsigned int p0, p1, p2;
	for (int x = 0; x < vertices.size(); ++x)
	{
		int index = x;
		pos = vertices[index].p;
		nor = vertices[index].normal;
		pos.w = 1.0f;
		nor.w = 0.0f;
		pos = modelMatrix * pos;
		nor = invModelMatrix * nor;
		vertices[index].p = pos;
		vertices[index].normal = nor;
		vertices[index].normal.make_unit_vector();
		minPoint.x = fmin(minPoint.x, pos.x);
		minPoint.y = fmin(minPoint.y, pos.y);
		minPoint.z = fmin(minPoint.z, pos.z);
		maxPoint.x = fmax(maxPoint.x, pos.x);
		maxPoint.y = fmax(maxPoint.y, pos.y);
		maxPoint.z = fmax(maxPoint.z, pos.z);
	}
	if (equal(minPoint.x, maxPoint.x))
	{
		minPoint.x -= 0.0001f;
		maxPoint.x += 0.0001f;
	}
	if (equal(minPoint.y, maxPoint.y))
	{
		minPoint.y -= 0.0001f;
		maxPoint.y += 0.0001f;
	}
	if (equal(minPoint.z, maxPoint.z))
	{
		minPoint.z -= 0.0001f;
		maxPoint.z += 0.0001f;
	}

	// face normal
	for (int x = 0; x < indices.size(); x += 3)
	{
		int index0 = indices[x + 0];
		int index1 = indices[x + 1];
		int index2 = indices[x + 2];
		vec3 normal = cross((vertices[index1].p - vertices[index0].p),
			vertices[index2].p - vertices[index0].p);
		normal.make_unit_vector();
		face_normal.push_back(normal);
	}

	m_box = aabb(minPoint, maxPoint);

}

#endif

