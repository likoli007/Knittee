#pragma once
#include <vector>
#include <qvector2d.h>
#include <qvector3d.h>
#include <Qt3DRender>
#include "Face.h"
#include <glm/glm.hpp>





class ObjectMesh
{
public:
    std::vector<QVector3D> vertices;
    std::vector<GLuint> indices;

    //std::vector<QVector3D> normal_indices;
    //std::vector<QVector2D> uv_indices;
public:
    void clear();

};


//helper struct for autoknit algorithms, first iteration had QVector3D as vertices, however some distance calculations seemed
//to have given wrong results, so switched to glm::vec3 for the duration of this algorithm, will switch back ObjectMesh when finished
struct AutoKnitMesh {
	std::vector<glm::vec3> vertices;
	std::vector<glm::uvec3> triangles;

	void clear() {
		vertices.clear();
		triangles.clear();
	}
	AutoKnitMesh() {
		vertices.clear();
		triangles.clear();
	}
	AutoKnitMesh(ObjectMesh obj) {
		for (QVector3D v : obj.vertices) {
			vertices.push_back(glm::vec3(v.x(), v.y(), v.z()));
		}

		for (int i = 0; i < obj.indices.size(); i += 3) {
			triangles.push_back(glm::uvec3(obj.indices[i], obj.indices[i + 1], obj.indices[i + 2]));
		}
	}
	ObjectMesh toObjMesh() {
		ObjectMesh res;
		for (glm::vec3 v : vertices) {
			res.vertices.push_back(QVector3D(v.x, v.y, v.z));
		}
		for (int i = 0; i < triangles.size(); i++) {
			res.indices.push_back(triangles[i].x);
			res.indices.push_back(triangles[i].y);
			res.indices.push_back(triangles[i].z);
		}
		return res;
	}
};