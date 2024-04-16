#include "ObjLoader.h"

#include<string>
#include<fstream>
#include <QVector3D>
#include <QVector2D>
#include <QFile>
#include <QString>
#include <iostream>
#include <QMessageBox>
#include "ObjectMesh.h"



/* 
*	Function: open an.obj file, read each line and process it using parseLine(),
*		and generate an ObjectMesh object after parsing using generateMesh()
*	Return: an ObjectMesh object that includes the vertex, normals etc. information of the loaded object
*/
ObjectMesh ObjLoader::loadFile(QString path)
{
	QFile file(path);
	QString line;
	ObjectMesh result;

	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		QMessageBox::information(0, "error", file.errorString());
		return result;
	}

	QTextStream in(&file);

	while (!in.atEnd()) {
		//load lines and parse them
		line = in.readLine();
		parseLine(line);
	}
	file.close();

	result = generateMesh();
	flushArrays();
	return result;
}

/*
*	Function: parse each line of an .obj file, generate new vertex/normals/face information based on type of line and append it to relevant vectors  
*	Return: no return value, each call of this function either appends one new entry to the loader or skips the line (in case of comments, empty lines etc.)
*	TODO: currently only triangle faces are supported, add support for quads, normals for shading purposes?
*		maybe texture information is not needed?
*/
void ObjLoader::parseLine(QString line)
{
	QStringList list;
	std::string type;

	// check for empty string, if empty skip parsing (empty line in .obj)
	if (!line.isEmpty()) {
		// vertex line, has 3 float values - coordinates of the vertex in 3D space - extract them and add them to vertex_indices
		if (line.startsWith("v ")) {
			list = line.simplified().split(' ');

			QVector3D vertex3(list[1].toFloat(), list[2].toFloat(), list[3].toFloat());
			qDebug() << list;
			//vertex3.normalize();
			vertex_list.push_back(vertex3);
		}

		// vertex texture, currently being loaded but may not be necessary, as the final texture of the object will just be yarn-like
		else if (line.startsWith("vt ")) {
			list = line.split(" ");
	
			QVector2D vertex2(list[1].toFloat(), list[2].toFloat());
			uv_indices.push_back(vertex2);
		}
		// vertex normals, currently unused but loaded anyways
		else if (line.startsWith("vn ")) {
			list = line.split(" ");

			QVector3D vertex3(list[1].toFloat(), list[2].toFloat(), list[3].toFloat());
			normal_indices.push_back(vertex3);
		}
		// face information, store in Face object format
		else if (line.startsWith("f ")) {
			Face face;
			int vertex_index;

			// most basic format, 'f a b c' where a,b,c is an index of the vertex in the file starting from 1
			if (!line.contains("/")) {
				list = line.split(" ");
				face.indices.push_back(list[1].toInt());
				face.indices.push_back(list[2].toInt());
				face.indices.push_back(list[3].toInt());
				faces.push_back(face);
			}
		}
	}
}

/*
	Function: generate an ObjectMesh object from the information already loaded into the loader
		that is, add the vertex coordinates information to the ObjectMesh in the order given by the face information
	Return: ObjectMesh object that holds all the information of the loaded .obj file
	TODO: quad support
*/
ObjectMesh ObjLoader::generateMesh()
{
	//iterate through the faces and construct a final representation
	ObjectMesh result;
	unsigned int vertex_index;

	//result.indices = faces;
	result.vertices = vertex_list;
	
	for (Face f : faces) {
		for (unsigned int i = 0; i < f.indices.size(); i++)
		{
			vertex_index = f.indices[i];
			GLuint index = f.indices[i]-1;

			//qDebug() << index;

			result.indices.push_back(index);

			//QVector3D vertex = vertex_list[vertex_index - 1];
			//result.vertices.push_back(vertex);
		}
	}

	//qDebug() << result.vertices.size();
	//qDebug() << result.indices.size();
	return result;
}


/*
*	Function: after the information is stored into the resulting ObjectMesh, flush the arrays that were used in the creation,
*		done so that loading a new object does not result in a mishmash of the previous and new object
*   Return: no return value
*/
void ObjLoader::flushArrays() {
	normal_indices.clear();
	vertex_list.clear();
	uv_indices.clear();
	faces.clear();
}