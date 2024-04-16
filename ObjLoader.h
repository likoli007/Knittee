#pragma once
#include <string>
#include <vector>
#include <QVector3D>
#include <QVector2D>
#include <QString>
#include "ObjectMesh.h"
#include "Face.h"



class ObjLoader
{
	std::string file_path;

	std::vector<QVector3D> normal_indices;
	std::vector<QVector3D> vertex_list;
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

