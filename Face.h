#pragma once
#include <vector>
#include <GL/gl.h>


//honestly not needed anymore since I am using the indices vector from ObjectMesh
class Face
{
public:
	std::vector<GLuint> indices;

};