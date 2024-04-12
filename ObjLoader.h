#pragma once
#include <string>
#include <vector>
#include <QVector3D>
#include <QVector2D>
#include <QString>
#include "ObjectMesh.h"

class Face
{
public:
	std::vector<int> vertices;
	std::vector<int> uvs;
	std::vector<int> normals;
	
};


class ObjLoader
{
	std::string file_path;

	std::vector<QVector3D> normal_indices;
	std::vector<QVector3D> vertex_indices;
	std::vector<QVector2D> uv_indices;
	std::vector<Face> faces;

public:
	ObjectMesh loadFile(QString path);
	void loadFileContents();
	void parseLine(QString line);
	ObjectMesh generateMesh();
private:
	void flushArrays();
};

