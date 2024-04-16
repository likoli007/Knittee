#pragma once
#include <vector>
#include <qvector2d.h>
#include <qvector3d.h>
#include <Qt3DRender>
#include "Face.h"





class ObjectMesh
{
public:
    std::vector<QVector3D> vertices;
    std::vector<GLuint> indices;

    //std::vector<QVector3D> normal_indices;
    //std::vector<QVector2D> uv_indices;
public:
    
};

