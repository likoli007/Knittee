#pragma once
#include <string>
#include <vector>
#include <QVector3D>
#include <QVector2D>
#include <QString>
#include "ObjectMesh.h"
#include "Face.h"
#include <QFile>


class ObjHandler
{
	QString path;
	QFile file;

	std::vector<QVector3D> normal_indices;
	std::vector<QVector3D> vertex_list;
	std::vector<QVector2D> uv_indices;
	std::vector<Face> faces;

public:
	ObjectMesh loadFile();
	void loadFileContents();
	void parseLine(QString line);
	void setFilePath(QString path);
	void copyObjFileToProject(QString projectName);
	ObjectMesh generateMesh();
private:
	void flushArrays();
};
