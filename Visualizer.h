#pragma once
#include<qopenglwidget.h>
#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include "ObjLoader.h"




class Visualizer :
    public QOpenGLWidget,
    protected QOpenGLFunctions
{
public:
    Visualizer(QWidget* parent = nullptr);
    ~Visualizer();
    void loadMesh(ObjectMesh loaded_mesh);

protected:
    QPoint lastMousePos;
    ObjectMesh mesh;
    QVector3D focusPoint, cameraLocation, up;
    QMatrix4x4 modelMatrix;
    bool isRotating = false;
    bool isZooming = false;
    bool isDragging = false;
    float translateX = 0.0f;
    float translateY = 0.0f;
    float translateZ = 0.0f;
    float m_rotationAngleX = 0.0f;
    float m_rotationAngleY = 0.0f;
    float zoomLevel = 1.0f;
    float rotationAngle = 0.0f;
    void initializeGL();
    void resizeGL(int w, int h);
    void paintGL();
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
};

