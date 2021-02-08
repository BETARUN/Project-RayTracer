#pragma once

#ifndef SKYBOXH
#define SKYBOXH

#include "triangle.h"
#include "model.h"
#include "texture.h"

class skybox
{
private:
	// front->back->left->right->top->bottom.
	std::vector<vertex> m_vertices;
	std::vector<unsigned int> m_indices;
	std::vector<texture*> m_cubemap;

public:
	skybox(const std::vector<texture*> &skymap);
	~skybox();

	vec3 sampleBackground(const ray &ray);

private:
	bool triangleHit(const ray &ray, const float &t_min, const float &t_max,
		hit_record &ret, const vertex &p0, const vertex &p1,
		const vertex &p2, const vec3 &normal) const;

};

skybox::skybox(const std::vector<texture*> &skymap)
{
	m_cubemap.push_back(skymap[0]);
	m_cubemap.push_back(skymap[1]);
	m_cubemap.push_back(skymap[2]);
	m_cubemap.push_back(skymap[3]);
	m_cubemap.push_back(skymap[4]);
	m_cubemap.push_back(skymap[5]);

	m_vertices.resize(24);
	m_indices.resize(36);

	const float size = 1.0f;
	// front
	m_vertices[0].p = vec3(+size, +size, +size);
	m_vertices[0].normal = vec3(+0, +0, -1);
	m_vertices[0].texcoord = vec2(0.0, 1.0);
	m_vertices[1].p = vec3(-size, +size, +size);
	m_vertices[1].normal = vec3(+0, +0, -1);
	m_vertices[1].texcoord = vec2(1.0, 1.0);
	m_vertices[2].p = vec3(-size, -size, +size);
	m_vertices[2].normal = vec3(+0, +0, -1);
	m_vertices[2].texcoord = vec2(1.0, 0.0);
	m_vertices[3].p = vec3(+size, -size, +size);
	m_vertices[3].normal = vec3(+0, +0, -1);
	m_vertices[3].texcoord = vec2(0.0, 0.0);
	m_indices[0] = 0; m_indices[1] = 1; m_indices[2] = 2;
	m_indices[3] = 0; m_indices[4] = 2; m_indices[5] = 3;

	// back
	m_vertices[4].p = vec3(+size, +size, -size);
	m_vertices[4].normal = vec3(+0, +0, +1);
	m_vertices[4].texcoord = vec2(1.0, 1.0);
	m_vertices[5].p = vec3(+size, -size, -size);
	m_vertices[5].normal = vec3(+0, +0, +1);
	m_vertices[5].texcoord = vec2(1.0, 0.0);
	m_vertices[6].p = vec3(-size, -size, -size);
	m_vertices[6].normal = vec3(+0, +0, +1);
	m_vertices[6].texcoord = vec2(0.0, 0.0);
	m_vertices[7].p = vec3(-size, +size, -size);
	m_vertices[7].normal = vec3(+0, +0, +1);
	m_vertices[7].texcoord = vec2(0.0, 1.0);
	m_indices[6] = 4; m_indices[7] = 5; m_indices[8] = 6;
	m_indices[9] = 4; m_indices[10] = 6; m_indices[11] = 7;

	// left
	m_vertices[8].p = vec3(-size, +size, +size);
	m_vertices[8].normal = vec3(+1, +0, +0);
	m_vertices[8].texcoord = vec2(0.0, 1.0);
	m_vertices[9].p = vec3(-size, +size, -size);
	m_vertices[9].normal = vec3(+1, +0, +0);
	m_vertices[9].texcoord = vec2(1.0, 1.0);
	m_vertices[10].p = vec3(-size, -size, -size);
	m_vertices[10].normal = vec3(+1, +0, +0);
	m_vertices[10].texcoord = vec2(1.0, 0.0);
	m_vertices[11].p = vec3(-size, -size, +size);
	m_vertices[11].normal = vec3(+1, +0, +0);
	m_vertices[11].texcoord = vec2(0.0, 0.0);
	m_indices[12] = 8; m_indices[13] = 9; m_indices[14] = 10;
	m_indices[15] = 8; m_indices[16] = 10; m_indices[17] = 11;

	// right
	m_vertices[12].p = vec3(+size, +size, -size);
	m_vertices[12].normal = vec3(-1, +0, +0);
	m_vertices[12].texcoord = vec2(0.0, 1.0);
	m_vertices[13].p = vec3(+size, +size, +size);
	m_vertices[13].normal = vec3(-1, +0, +0);
	m_vertices[13].texcoord = vec2(1.0, 1.0);
	m_vertices[14].p = vec3(+size, -size, +size);
	m_vertices[14].normal = vec3(-1, +0, +0);
	m_vertices[14].texcoord = vec2(1.0, 0.0);
	m_vertices[15].p = vec3(+size, -size, -size);
	m_vertices[15].normal = vec3(-1, +0, +0);
	m_vertices[15].texcoord = vec2(0.0, 0.0);
	m_indices[18] = 12; m_indices[19] = 13; m_indices[20] = 14;
	m_indices[21] = 12; m_indices[22] = 14; m_indices[23] = 15;

	// top
	m_vertices[16].p = vec3(+size, +size, -size);
	m_vertices[16].normal = vec3(+0, -1, +0);
	m_vertices[16].texcoord = vec2(1.0, 0.0);
	m_vertices[17].p = vec3(-size, +size, -size);
	m_vertices[17].normal = vec3(+0, -1, +0);
	m_vertices[17].texcoord = vec2(0.0, 0.0);
	m_vertices[18].p = vec3(-size, +size, +size);
	m_vertices[18].normal = vec3(+0, -1, +0);
	m_vertices[18].texcoord = vec2(0.0, 1.0);
	m_vertices[19].p = vec3(+size, +size, +size);
	m_vertices[19].normal = vec3(+0, -1, +0);
	m_vertices[19].texcoord = vec2(1.0, 1.0);
	m_indices[24] = 16; m_indices[25] = 17; m_indices[26] = 18;
	m_indices[27] = 16; m_indices[28] = 18; m_indices[29] = 19;

	// bottom
	m_vertices[20].p = vec3(+size, -size, -size);
	m_vertices[20].normal = vec3(+0, +1, +0);
	m_vertices[20].texcoord = vec2(1.0, 1.0);
	m_vertices[21].p = vec3(+size, -size, +size);
	m_vertices[21].normal = vec3(+0, +1, +0);
	m_vertices[21].texcoord = vec2(1.0, 0.0);
	m_vertices[22].p = vec3(-size, -size, +size);
	m_vertices[22].normal = vec3(+0, +1, +0);
	m_vertices[22].texcoord = vec2(0.0, 0.0);
	m_vertices[23].p = vec3(-size, -size, -size);
	m_vertices[23].normal = vec3(+0, +1, +0);
	m_vertices[23].texcoord = vec2(0.0, 1.0);
	m_indices[30] = 20; m_indices[31] = 21; m_indices[32] = 22;
	m_indices[33] = 20; m_indices[34] = 22; m_indices[35] = 23;
}

skybox::~skybox()
{

}

vec3 skybox::sampleBackground(const ray &skybox_ray)
{
	hit_record rec;
	ray r(vec3(0, 0, 0), skybox_ray.direction());
	int index = -1;
	for (int x = 0; x < m_indices.size(); x += 3)
	{
		int index1 = m_indices[x + 0];
		int index2 = m_indices[x + 1];
		int index3 = m_indices[x + 2];
		if (triangleHit(r, 0.001f, FLT_MAX, rec,
			m_vertices[index1], m_vertices[index2], m_vertices[index3],
			m_vertices[index1].normal))
		{
			index = x;
			break;
		}
	}

	if (index != -1)
	{
		int map = index / 6;
		return m_cubemap[map]->value(rec.texcoord.x, rec.texcoord.y, rec.p);
	}
	else
		return vec3(0.0, 0.0, 0.0);
}

bool skybox::triangleHit(const ray &ray, const float &t_min, const float &t_max,
	hit_record &ret, const vertex &p0, const vertex &p1,
	const vertex &p2, const vec3 &normal) const
{
	float n_dot_dir = dot(normal, ray.direction());
	// no intersection.
	if (equal(n_dot_dir, 0.0))
		return false;
	float d = -dot(normal, p0.p);
	float t = -(dot(normal, ray.origin()) + d) / n_dot_dir;
	if (t < t_min || t > t_max)
		return false;
	ret.t = t;
	ret.p = ray.point_at_parameter(t);
	ret.mat_ptr = nullptr;
	// judge inside or not.
	vec3 r = ret.p - p0.p;
	vec3 q1 = p1.p - p0.p;
	vec3 q2 = p2.p - p0.p;
	float q1_squaredLen = q1.length() * q1.length();
	float q2_squaredLen = q2.length() * q2.length();
	float q1_dot_q2 = dot(q1, q2);
	float r_dot_q1 = dot(r, q1);
	float r_dot_q2 = dot(r, q2);
	float determinant = 1.0f / (q1_squaredLen * q2_squaredLen - q1_dot_q2 * q1_dot_q2);

	//
	float omega1 = determinant * (q2_squaredLen * r_dot_q1 - q1_dot_q2 * r_dot_q2);
	float omega2 = determinant * (-q1_dot_q2 * r_dot_q1 + q1_squaredLen * r_dot_q2);
	if (omega1 + omega2 > 1.0f || omega1 < 0.0f || omega2 < 0.0f)
		return false;
	ret.texcoord = p0.texcoord * (1.0f - omega1 - omega2) + p1.texcoord * omega1 + p2.texcoord * omega2;
	return true;
}

#endif