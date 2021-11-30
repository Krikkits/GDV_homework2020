/**
 * Computer Graphics Lecture WS 2009/2010 Ulm University
 * Creator: Manuel Finckh
 * Email:   manuel.finckh@uni-ulm.de
 */

#include "scene.h"
#include "utils/vec.h"
#include "utils/fileio.h"
#include <fstream>
#include <string>
#include <cstdio>
#include <iostream>
#include <vector>

using namespace std;

Scene::Scene(const char* sceneFile, const char* envFile)
{
	LoadScn(sceneFile);
	LoadEnv(envFile);
}

Scene::~Scene()
{
	delete[] triangles;
	delete[] normals;
	delete[] material;
	delete[] mat_index;
	if (environment != 0)
		delete[] environment;
	delete[] cam;
}

bool Scene::LoadScn(const char * file)
{
	string sceneFilename = string(file) + string(".scn");
	string meshFilename = string(file) + string(".ra2");
	string normalFilename = string(file) + string(".n");

	ifstream sceneFile(sceneFilename.c_str());
	ifstream meshFile(sceneFilename.c_str());
	ifstream normalFile(sceneFilename.c_str());

	if (!sceneFile || !meshFile || !normalFile)
		return false;

	meshFile.close();
	normalFile.close();

	// Load mesh and normals
	float* tris;
	float* norms;
	num_tris = load_float_data(meshFilename.c_str(), tris) / (3 * 3);
	triangles = (Triangle*) tris;
	load_float_data(normalFilename.c_str(), norms);
	normals = (Vec3*) norms;

	mat_index = new int[num_tris];

	vector<Material> tmpMaterials;
	int ResX = 300, ResY = 300;

	string line;
	while (getline(sceneFile, line))
	{
		if (line.find("#") == 0 || line.length() == 0)
			continue;
		if (line.find("Width") == 0)
		{
			sscanf(line.c_str(), "%*s %d", &ResX);
			continue;
		}
		if (line.find("Height") == 0)
		{
			sscanf(line.c_str(), "%*s %d", &ResY);
			continue;
		}
		if (line.find("Material") == 0)
		{
			int start;
			int end;
			sscanf(line.c_str(), "%*s %d %d", &start, &end);
			for (int t = start; t < end; t++)
				mat_index[t] = tmpMaterials.size();

			Material mat;
			getline(sceneFile, line);
			sscanf(line.c_str(), "%*s %f %f %f", &mat.color_d[0],
					&mat.color_d[1], &mat.color_d[2]);
			getline(sceneFile, line);
			sscanf(line.c_str(), "%*s %f %f %f", &mat.color_e[0],
					&mat.color_e[1], &mat.color_e[2]);
			tmpMaterials.push_back(mat);

			continue;
		}
	}

	material = new Material[tmpMaterials.size()];
	for (unsigned int m = 0; m < tmpMaterials.size(); m++)
		material[m] = tmpMaterials[m];

	sceneFile.close();

	cam = new Cam(Triangle::getAABB(triangles, num_tris), ResX, ResY);
	return true;
}

bool Scene::LoadEnv(const char * file)
{
	environment = 0;

	if (file == 0)
		return false;

	float *data;

	bool success = load_image_hdr(file, data, env_x, env_y);
	environment = (Vec3*) data;

	return success;
}
