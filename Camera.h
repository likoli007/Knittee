// Camera.h
#ifndef CAMERA_H
#define CAMERA_H

#include <QVector3D>
#include <QMatrix4x4>
#include <QMouseEvent>

class Camera {
public:
    Camera();

    void setAspectRatio(float aspectRatio);
    void zoom(float delta);
    void rotate(QPoint delta);
    void pan(QPoint delta);

    QMatrix4x4 getViewMatrix();

public:
    QVector3D position;
    QVector3D up;
    QVector3D front;

    float fov;
    float aspectRatio;
    float nearPlane;
    float farPlane;

    void updateVectors();
};

#endif // CAMERA_H