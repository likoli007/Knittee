#pragma once
#include <vector>
#include <qvector2d.h>
#include <qvector3d.h>
#include <Qt3DRender>

class ObjectMesh
{
public:
    std::vector<QVector3D> vertices;
    std::vector<QVector3D> normal_indices;
    std::vector<QVector2D> uv_indices;
public:
    
};

