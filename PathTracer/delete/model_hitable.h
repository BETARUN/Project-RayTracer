#pragma once
#ifndef MODELHITABLEH
#define MODELHITABLEH

#include "mesh_hitable.h"
#include <fstream>
#include <sstream>

using namespace std;

class model_hitable : public mesh_hitable
{
private:
	vec3 model_center;
	vec3 model_scale;

public:
	model_hitable(const std::string &path, vec3 pos, vec3 len, material* mat);
	virtual ~model_hitable() = default;

private:
	// Obj file loader.
	void loadObjFile(const std::string &path);
};

model_hitable::model_hitable(const std::string &path, vec3 pos, vec3 len, material* mat)
{
	mat_ptr = mat;
	loadObjFile(path);
	translate(pos - model_center);
	scale(len /** Vector3D(1.0f/m_scale.x, 1.0f/m_scale.y, 1.0f/m_scale.z)*/);
}


void model_hitable::loadObjFile(const std::string &path)
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

#endif