#pragma once
#include<qopenglwidget.h>
#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFunctions>
#include "ObjHandler.h"




struct Constraint
{
    std::vector<GLuint> vertices;
    float timeValue = 0.0f;
};

class Visualizer :
    public QOpenGLWidget,
    protected QOpenGLExtraFunctions
{
public:
    Visualizer(QWidget* parent = nullptr);
    ~Visualizer();
    void loadMesh(ObjectMesh loaded_mesh);
    void setConstraintsMode();
    std::vector<Constraint*> getConstraints();
    void setConstraints(std::vector<Constraint*>);

protected:
    int mouseX, mouseY;

    std::vector<Constraint*> constraints;
    Constraint* currentConstraint;

    QPoint lastMousePos;
    ObjectMesh mesh;
    QVector3D focusPoint, cameraLocation, up;
    QMatrix4x4 modelMatrix;
    QMatrix4x4 viewMatrix;
    QMatrix4x4 projectionMatrix;
    QMatrix4x4 modelViewMatrix;
    QMatrix4x4 mvpMatrix;

    QOpenGLShaderProgram shaderProgram;

    GLuint vao, vbo, indexBuffer, shadingBuffer;

    int selectedFace = -1;
    int chosenVertex = -1;
    QVector3D rayDir;
    QPointF wpos;

    std::vector<QVector3D> origins;
    std::vector<QVector3D> directions;
    QVector3D rayOrigin;
    QVector3D rayDirection;
    bool rayDebug = false;
    bool getFace = false;
    bool objectLoaded = false;
    bool isRotating = false;
    bool isZooming = false;
    bool isDragging = false;
    bool addingConstraints = false;
    bool showConstraints = false;            //both of these values should be set from the toolbar inside Knittee
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
    void paintConstraints();
    void paintPickFrame();
    void paintConstraintHighlight(QVector3D, QVector3D);
    void paintPickedFace(int);
    void paintMesh();



    void buildmvpMatrix();
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
    void keyPressEvent(QKeyEvent* event);
    void computeRayFromMouse(QPoint currentPosition);

    QVector3D getMouseCoordinatesInWorld();
    void addConstraintVertex();

    void pushConstraints();
    void pickFromMesh();



};

