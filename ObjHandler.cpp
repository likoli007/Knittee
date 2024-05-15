#include "ObjHandler.h"

#include<string>
#include<fstream>
#include <QVector3D>
#include <QVector2D>
#include <QFile>
#include <QString>
#include <iostream>
#include <QMessageBox>
#include "ObjectMesh.h"

#include <QCoreApplication>




void ObjHandler::normalizeMesh(ObjectMesh& mesh) {
	float minX = std::numeric_limits<float>::max();
	float minY = std::numeric_limits<float>::max();
	float minZ = std::numeric_limits<float>::max();
	float maxX = std::numeric_limits<float>::lowest();
	float maxY = std::numeric_limits<float>::lowest();
	float maxZ = std::numeric_limits<float>::lowest();

	for (const auto& vertex : mesh.vertices) {
		if (vertex.x() < minX) minX = vertex.x();
		if (vertex.y() < minY) minY = vertex.y();
		if (vertex.z() < minZ) minZ = vertex.z();
		if (vertex.x() > maxX) maxX = vertex.x();
		if (vertex.y() > maxY) maxY = vertex.y();
		if (vertex.z() > maxZ) maxZ = vertex.z();
	}


	float currentSizeX = maxX - minX;
	float currentSizeY = maxY - minY;
	float currentSizeZ = maxZ - minZ;

	float scaleFactor = desiredSize / std::max({ currentSizeX, currentSizeY, currentSizeZ });

	qDebug() << "current dimensions: " << currentSizeX << " " << currentSizeY << " " << currentSizeZ;
	qDebug() << "scale factor: " << scaleFactor;

	// Scale the mesh vertices
	for (auto& vertex : mesh.vertices) {
		vertex.setX(vertex.x() * scaleFactor);
		vertex.setY(vertex.y() * scaleFactor);
		vertex.setZ(vertex.z() * scaleFactor);
	}

}



/* 
*	Function: open an.obj file, read each line and process it using parseLine(),
*		and generate an ObjectMesh object after parsing using generateMesh()
*	Return: an ObjectMesh object that includes the vertex, normals etc. information of the loaded object
*/


void ObjHandler::copyObjFileToProject(QString projectName) {

	
	QString destinationFolderPath = QCoreApplication::applicationDirPath() + "/Projects/" + projectName;
	qDebug() << destinationFolderPath;
	//QDir().mkpath(destinationFolderPath);
	QFile file(path);
	QFile newFile(destinationFolderPath+ "/mesh.obj");
	
	
	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		QMessageBox::information(0, "error", file.errorString());
		return;
	}

	if (!newFile.open(QIODevice::WriteOnly | QIODevice::Text)) {
		QMessageBox::information(0, "error", newFile.errorString());
		return;
	}

	QTextStream in(&file);
	QTextStream out(&newFile);

	while (!in.atEnd()) {
		QString line = in.readLine();
		out << line << "\n";
	}

	file.close();
	newFile.close();

	
	

	//newFile.close();
}
void ObjHandler::setFilePath(QString path) {
	this->path = path;
	file.setFileName(path);
}

ObjectMesh ObjHandler::loadFile()
{
	//QFile file(path);
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
void ObjHandler::parseLine(QString line)
{
	QStringList list;
	std::string type;

	// check for empty string, if empty skip parsing (empty line in .obj)
	if (!line.isEmpty()) {
		// vertex line, has 3 float values - coordinates of the vertex in 3D space - extract them and add them to vertex_indices
		if (line.startsWith("v ")) {
			list = line.simplified().split(' ');

			QVector3D vertex3(list[1].toFloat(), list[2].toFloat(), list[3].toFloat());
			//qDebug() << list;
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
			//Face face;
			//int vertex_index;

			int vertices = line.split(" ").size()-1;

			// most basic format, 'f a b c' where a,b,c is an index of the vertex in the file starting from 1
			if (vertices == 3) {
				//qDebug() << "triangle";
				addTriangle(line);
				/*list = line.split(" ");
				face.indices.push_back(list[1].toInt());
				face.indices.push_back(list[2].toInt());
				face.indices.push_back(list[3].toInt());
				faces.push_back(face);*/
			}
			else if (vertices == 4) {
				addQuad(line);

				/*
				list = line.split(" ");
				face.indices.push_back(list[1].toInt());
				face.indices.push_back(list[2].toInt());
				face.indices.push_back(list[3].toInt());
				face.indices.push_back(list[4].toInt());
				faces.push_back(face);*/
			}
			else {
				qDebug() << "non triangle or quad!";
			}
		}
	}
}

void ObjHandler::addTriangle(QString line) {
	Face face;
	QStringList list = line.split(" ");
	if (list.size() == 4) {
		QString a = list[1];
		QString b = list[2];
		QString c = list[3];

		QStringList a_list = a.split("/");
		QStringList b_list = b.split("/");
		QStringList c_list = c.split("/");

		face.indices.push_back(a_list[0].toInt());
		face.indices.push_back(b_list[0].toInt());
		face.indices.push_back(c_list[0].toInt());
	}
	faces.push_back(face);
}

void ObjHandler::addQuad(QString line) {
	Face face1, face2;
	QStringList list = line.split(" ");
	if (list.size() == 5) {
		QString a = list[1];
		QString b = list[2];
		QString c = list[3];
		QString d = list[4];

		QStringList a_list = a.split("/");
		QStringList b_list = b.split("/");
		QStringList c_list = c.split("/");
		QStringList d_list = d.split("/");

		face1.indices.push_back(a_list[0].toInt());
		face1.indices.push_back(b_list[0].toInt());
		face1.indices.push_back(c_list[0].toInt());

		face2.indices.push_back(a_list[0].toInt());
		face2.indices.push_back(c_list[0].toInt());
		face2.indices.push_back(d_list[0].toInt());

		//qDebug() << "face1" << face1.indices[0] << face1.indices[1] << face1.indices[2];
		//qDebug() << "face2" << face2.indices[0] << face2.indices[1] << face2.indices[2];

		faces.push_back(face1);
		faces.push_back(face2);
	}
}



/*
	Function: generate an ObjectMesh object from the information already loaded into the loader
		that is, add the vertex coordinates information to the ObjectMesh in the order given by the face information
	Return: ObjectMesh object that holds all the information of the loaded .obj file
	TODO: quad support
*/
ObjectMesh ObjHandler::generateMesh()
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

	normalizeMesh(result);

	//qDebug() << result.vertices.size();
	//qDebug() << result.indices.size();
	return result;
}


/*
*	Function: after the information is stored into the resulting ObjectMesh, flush the arrays that were used in the creation,
*		done so that loading a new object does not result in a mishmash of the previous and new object
*   Return: no return value
*/
void ObjHandler::flushArrays() {
	normal_indices.clear();
	vertex_list.clear();
	uv_indices.clear();
	faces.clear();
}