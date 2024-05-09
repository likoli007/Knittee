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
#include "Link.h"
#include "FlatPoint.h"

struct Constraint
{
    std::vector<GLuint> vertices;
    float timeValue = 0.0f;
};

class Visualizer :
    public QOpenGLWidget,
    protected QOpenGLExtraFunctions
{
	Q_OBJECT
public:
    Visualizer(QWidget* parent = nullptr);
    ~Visualizer();
    void loadMesh(ObjectMesh loaded_mesh);
    void setConstraintsMode();
    std::vector<Constraint*> getConstraints();
    void setConstraints(std::vector<Constraint*>);
    void peelSliceDone(ObjectMesh*, std::vector< std::vector< uint32_t > >*, std::vector< std::vector< uint32_t > >*);
    void linkChainsDone(std::vector< std::vector< Stitch > >*, std::vector< Link >*);
    void nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex>>*);

    //could be getted and setted in the future
    int projectType = 0; //is the user operating on a 3D model (0) or a 2D sheet? (1)

public slots:
    void meshInterpolated(ObjectMesh mesh, std::vector<float> values);
    void firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
        std::vector< std::vector< Stitch > >* active_stitches,
        RowColGraph* graph);
    void sheetChanged(std::vector<std::vector<FlatPoint>>*);
    void showInterpolatedChanged(int state);
    void showGraphChanged(int state);
    void showTracedChanged(int state);
protected:
    std::vector<std::vector<FlatPoint>> sheet;
    
    
    void paintSheet();
    bool isMovingLoop = false;
    QPair<int, int> getLoopPos(int, int);
    QPair<int, int> from, to;
    
    
    
    int mouseX, mouseY;

   
    int timerID = 0;
    void setTimer();
    void timerEvent(QTimerEvent* event);
    void stopTimer();



    std::vector<Constraint*> constraints;
    Constraint* currentConstraint;

    ObjectMesh interpolatedMesh;
    ObjectMesh sliceMesh;
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

    GLuint originalVao, originalVbo, originalIndexBuffer;
    GLuint interpolatedVao, interpolatedVbo, interpolatedIndexBuffer;
    GLuint sliceVao, sliceVbo, sliceIndexBuffer;


    int selectedFace = -1;
    int chosenVertex = -1;
    QVector3D rayDir;
    QPointF wpos;



    std::vector< std::vector< EmbeddedVertex > > activeChains;
    std::vector< std::vector< Stitch > > activeStitches;
    RowColGraph rowColGraph;
    std::vector< std::vector< uint32_t > > sliceActiveChains;
    std::vector< std::vector< uint32_t > > sliceNextChains;
    std::vector< std::vector< Stitch > > nextStitches;
    std::vector< Link > links;
    std::vector<std::vector<EmbeddedVertex>> nextActiveChains;



    QImage timeTexImage;

    std::vector<QVector3D> origins;
    std::vector<QVector3D> directions;
    QVector3D rayOrigin;
    QVector3D rayDirection;
    bool rayDebug = false;
    bool getFace = false;
   
    bool isRotating = false;
    bool isZooming = false;
    bool isDragging = false;
    bool addingConstraints = false;
   
    
    bool showModel = false;
    bool showFirstChains = false;
    bool showInterpolated = false;
    bool showSlice = false;
    bool showConstraints = false;            //both of these values should be set from the toolbar inside Knittee
    bool showLinks = false;
    bool showNextChains = false;
    bool showSliceChains = false;
    bool showGraph = false;
    bool showTraced = false;

    float translateX = 0.0f;
    float translateY = 0.0f;
    float translateZ = -5.0f;
    float m_rotationAngleX = 0.0f;
    float m_rotationAngleY = 0.0f;
    float zoomLevel = 1.0f;
    float rotationAngle = 0.0f;
    void initializeGL();
    void resizeGL(int w, int h);



    void paintYarn(int x, int y, int ex, int ey);
    void paintGL();
    void paintConstraints();
    void paintPickFrame();
    void paintPath(QVector3D, QVector3D, float, float, float);
    void paintPickedFace(int);
    void paintMesh();
    void paintInterpolatedMesh();
    void paintFirstActiveChains();
    void paintSliceMesh();
    void paintSphere(QVector3D, float);
    void paintLinks();
    void paintNextChains();

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

