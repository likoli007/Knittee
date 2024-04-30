#pragma once
#include<qopenglwidget.h>
#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFunctions>
#include "ObjHandler.h"

#include "EmbeddedVertex.h"
#include "Stitch.h"
#include "RowColGraph.h"


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
public slots:
    void meshInterpolated(ObjectMesh mesh, std::vector<float> values);
    void firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
        std::vector< std::vector< Stitch > >* active_stitches,
        RowColGraph* graph);
protected:
    int mouseX, mouseY;

    std::vector<Constraint*> constraints;
    Constraint* currentConstraint;

    ObjectMesh interpolatedMesh;
    std::vector<float> interpolatedValues;

    std::vector<std::array<float, 3>> colors;

    QPoint lastMousePos;
    ObjectMesh mesh;
    QVector3D focusPoint, cameraLocation, up;
    QMatrix4x4 modelMatrix;
    QMatrix4x4 viewMatrix;
    QMatrix4x4 projectionMatrix;
    QMatrix4x4 modelViewMatrix;
    QMatrix4x4 mvpMatrix;

    QOpenGLShaderProgram shaderProgram;
    QOpenGLShaderProgram timeShaderProgram;
    QOpenGLTexture* timeTexture;

    GLuint vao, vbo, indexBuffer, shadingBuffer;
    GLuint colorBuffer;
    GLuint coordsBuffer;

    int selectedFace = -1;
    int chosenVertex = -1;
    QVector3D rayDir;
    QPointF wpos;



    std::vector< std::vector< EmbeddedVertex > > activeChains;
    std::vector< std::vector< Stitch > > activeStitches;
    RowColGraph rowColGraph;






    QImage timeTexImage;

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
    

    bool firstChainsLoaded = false;
    bool interpolatedLoaded = false;
    
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
    void paintConstraintHighlight(QVector3D, QVector3D, float, float);
    void paintPickedFace(int);
    void paintMesh();
    void paintInterpolatedMesh();
    void paintFirstActiveChains();


    void buildmvpMatrix();
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
    void keyPressEvent(QKeyEvent* event);
    void computeRayFromMouse(QPoint currentPosition);

    void computeColors();

    void generateTimeTexture();

    QVector3D getMouseCoordinatesInWorld();
    void addConstraintVertex();

    void pushConstraints();
    void pickFromMesh();

    void increaseConstraintValue();
    void decreaseConstraintValue();
    int getClosestConstraint();
};

