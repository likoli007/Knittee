#pragma once
#include<qopenglwidget.h>
#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLExtraFunctions>
#include <QOpenGLFunctions>
#include <GL/gl.h>
#include <QDebug>
#include <QToolTip>
#include <QApplication>
#include <QPair>
#include <QScreen>
#include <QGraphicsDropShadowEffect>

#include "ObjHandler.h"
#include "EmbeddedVertex.h"
#include "Stitch.h"
#include "RowColGraph.h"
#include "Link.h"
#include "FlatPoint.h"
#include "TracedStitch.h"

/*
*   This struct represents a set of vertices on which a constraint resides
*   timeValue is the value of the constraint in the range [-1, 1]
*/

struct Constraint
{
    std::vector<GLuint> vertices;
    float timeValue = 0.0f;
};


/*
*   Visualizer is in charge of using OpenGL to show all the different visualization methods that the program has
*       It is also in charge of visual based options (such as adding constraint), as well as keyboard functions
*/
class Visualizer :
    public QOpenGLWidget,
    protected QOpenGLExtraFunctions
{
	Q_OBJECT
public:
    Visualizer(QWidget* parent = nullptr);
    ~Visualizer();
    void loadMesh(ObjectMesh loaded_mesh);
    void reset();
    std::vector<Constraint*> getConstraints();
    void setConstraints(std::vector<Constraint*>);
    void linkChainsDone(std::vector< std::vector< Stitch > >*, std::vector< Link >*);
    void nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex>>*);
    void knitGraphCreated();
    //could be getted and setted in the future
    int projectType = 0; //is the user operating on a 3D model (0) or a 2D sheet? (1)

signals:
    void requestConstraintsSave();
    void moveLoop(QPair<int, int>, QPair<int, int>);

public slots:
    void peelSliceDone(ObjectMesh, std::vector< std::vector< uint32_t > >, std::vector< std::vector< uint32_t > >);
    void setConstraintsMode(bool);
    void meshInterpolated(ObjectMesh mesh, std::vector<float> values);
    void firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
        std::vector< std::vector< Stitch > >* active_stitches,
        RowColGraph* graph);
    void sheetChanged(std::vector<std::vector<FlatPoint>>*);
    void showInterpolatedChanged(int state);
    void showGraphChanged(int state);
    void showTracedChanged(int state);
    void showYarnChanged(int state);
    void knitGraphTraced(std::vector< TracedStitch >*);
protected:
    //Color information for various steps
    glm::vec4 clearColor = { 0.65f, 0.65f, 0.65f, 1.0f };
    glm::vec3 graphRowColor = { 0.3, 0.3, 0.3 };//{ 0.44, 0.16, 0.39 };
    glm::vec3 graphLinkColor = { 0.196, 0.710, 0.176 };
    glm::vec3 graphCollapseColor = { 1.0 , 0.34, 0.2 };
    glm::vec3 graphExpandColor = { 0.953, 0.953, 0.110 };
    glm::vec3 lowerBoundColor = { 0.090, 0.216, 0.753 };
    glm::vec3 upperBoundColor = { 0.780, 0.0, 0.224 };
    glm::vec3 yarnMeshColor = { 0.831, 0.110, 0.447 };
    std::vector<glm::vec3> ycolors = { glm::vec3(0.831, 0.110, 0.447), glm::vec3(0.098, 0.580, 0.129),
        glm::vec3(0.137, 0.275, 0.737), glm::vec3(0.745, 0.420, 0.749), glm::vec3(0.808, 0.851, 0.082),
        glm::vec3(0.082, 0.851, 0.835), glm::vec3(0.898, 0.886, 0.863) };
    std::vector<std::array<float, 3>> colors;
    //used for mouse events
    int mouseX, mouseY;
    QPoint lastMousePos;

    //timer information, going for 60fps visualization
    int timerID = 0;
    void setTimer();
    void timerEvent(QTimerEvent* event);
    void stopTimer();

    //2D sheet information
    std::vector<std::vector<FlatPoint>> sheet;
    void paintSheet();
    bool isMovingLoop = false;
    QPair<int, int> getLoopPos(int, int);
    QPair<int, int> from, to;
    
    //3D mesh information
    std::vector<Constraint*> constraints;
    Constraint* currentConstraint;
    ObjectMesh mesh;
    ObjectMesh interpolatedMesh;
    ObjectMesh sliceMesh;
    std::vector<float> interpolatedValues;
    std::vector< TracedStitch > tracedMesh;
    std::vector< std::vector< EmbeddedVertex > > activeChains;
    std::vector< std::vector< Stitch > > activeStitches;
    RowColGraph rowColGraph;
    std::vector< std::vector< uint32_t > > sliceActiveChains;
    std::vector< std::vector< uint32_t > > sliceNextChains;
    std::vector< std::vector< Stitch > > nextStitches;
    std::vector< Link > links;
    std::vector<std::vector<EmbeddedVertex>> nextActiveChains;

    QMatrix4x4 modelMatrix;
    QMatrix4x4 viewMatrix;
    QMatrix4x4 projectionMatrix;
    QMatrix4x4 mvpMatrix;

    QOpenGLShaderProgram shaderProgram;
    QOpenGLShaderProgram timeShaderProgram;
    QOpenGLTexture* timeTexture;

    GLuint vao, vbo, indexBuffer;
    GLuint coordsBuffer;

    int selectedFace = -1;
    int chosenVertex = -1;
    QImage timeTexImage;
    const uint TimeTexSize = 16;               //how many discrete colors are there? 16 seems to be enough, could be made more granular
   
    bool isRotating = false;
    bool isZooming = false;
    bool isDragging = false;
    bool addingConstraints = false;
   
    bool showLastChain = false;
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
    bool showYarn = false;

    float translateX = 0.0f;
    float translateY = 0.0f;
    float translateZ = -20.0f;
    float m_rotationAngleX = 0.0f;
    float m_rotationAngleY = 0.0f;
    float zoomLevel = 1.0f;
    float zoomFactor = 10.0f;

    QVector3D centroid;

    void initializeGL();
    void resizeGL(int w, int h);
    void computeBoundaries();

    void paintYarn(int x, int y, int ex, int ey);
    void paintGL();
    void paintConstraints();
    void paintPickFrame();
    void paintPath(QVector3D, QVector3D, float, float, float, float);
    void paintPickedFace(int);
    void paintOriginalMesh();
    void paintInterpolatedMesh();
    void paintFirstActiveChains();
    void paintSliceMesh();
    void paintLinks();
    void paintNextChains();
    void paintTraced();
    void paintCurve(QPainter& painter, const std::vector<QPoint>& controlPoints);
    void paintYarnMesh();
    void paint3DCurve(QVector3D start, QVector3D end, QVector3D cp1, QVector3D cp2);

    void buildmvpMatrix();
    QColor timeColor(float time);
    void mousePressEvent(QMouseEvent* event);
    void mouseReleaseEvent(QMouseEvent* event);
    void mouseMoveEvent(QMouseEvent* event);
    void wheelEvent(QWheelEvent* event);
    void keyPressEvent(QKeyEvent* event);
    void generateTimeTexture();
    void addConstraintVertex();
    void pushConstraints();
    void pickFromMesh();
    void increaseConstraintValue();
    void decreaseConstraintValue();
    int getClosestConstraint();
    void deleteConstraint();
    void loadInterpolated();

    std::vector< std::vector< QVector3D > > copyLocations(ObjectMesh const& model, std::vector< std::vector< uint32_t > > const& chains);
    QVector3D mmix(const QVector3D& wa, const QVector3D& wb, float m);
    std::vector< std::vector< QVector3D > > interpolateStitchLocations(std::vector< std::vector< QVector3D > > const& chains, std::vector< std::vector< Stitch > > const& stitches);
    std::vector< std::vector< QVector3D > > interpolateLocations(ObjectMesh const& model, std::vector< std::vector< EmbeddedVertex > > const& chains);


};

