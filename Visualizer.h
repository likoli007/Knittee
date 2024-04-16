#pragma once
#include<qopenglwidget.h>
#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFunctions>
#include "ObjLoader.h"




class Visualizer :
    public QOpenGLWidget,
    protected QOpenGLExtraFunctions
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
    QMatrix4x4 viewMatrix;
    QMatrix4x4 projectionMatrix;
    QMatrix4x4 mvpMatrix;

    QOpenGLShaderProgram shaderProgram;

    GLuint vao, vbo, indexBuffer, shadingBuffer;

    QOpenGLFramebufferObject *fbo;


    std::vector<QVector3D> origins;
    std::vector<QVector3D> directions;
    QVector3D rayOrigin;
    QVector3D rayDirection;
    bool rayDebug = false;

    bool objectLoaded = false;
    bool isRotating = false;
    bool isZooming = false;
    bool isDragging = false;
    float translateX = 0.0f;
    float translateY = 0.0f;
    float translateZ = -5.0f;
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
    void computeRayFromMouse(QPoint currentPosition);
    void drawPickFrame(int, int);
};

