#include "Visualizer.h"
#include <QWidget>
#include <QOpenGLWidget>
#include <QOpenGLFunctions>
#include <QOpenGLExtraFunctions>
#include <GL/gl.h>
#include "ObjLoader.h"


/*
*   Function: constructor, adds mouse tracking for the visualization subwindow
*/
Visualizer::Visualizer(QWidget* parent )
    : QOpenGLWidget(parent)
{
    setMouseTracking(true);

}

/*
*   Function: destructor
*/
Visualizer::~Visualizer()
{

}


/*
*   Function: initialize the OpenGL graphics API, set some default parameters
*   Return: no return values
*/
void Visualizer::initializeGL()
{
    initializeOpenGLFunctions();


    glMatrixMode(GL_PROJECTION);
    glEnable(GL_DEPTH_TEST);
    glLoadIdentity();
    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    glOrtho(-1 * aspectRatio, 1 * aspectRatio, -1, 1, -1000, 1000); 
    glMatrixMode(GL_MODELVIEW);
    
    glLoadIdentity();
}

/*
*   Function: inherited function that is used when the user resizes the application
*/
void Visualizer::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
}

/*
*   Function: paint the OpenGL scene, that is the 3D object selected by the user 
*   Return: no return values
*   TODO: quad having objects may need to be drawn in a different method (GL_TRIANGLES vs GL_QUADS?)
*       currently the color is just a simple gradient, add shading and different texture
*/
void Visualizer::paintGL()
{

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // apply zoom, rotation, drag transformations on the view matrix
    glTranslatef(translateX, translateY, translateZ);
    glRotatef(m_rotationAngleX, 1.0f, 0.0f, 0.0f);
    glRotatef(m_rotationAngleY, 0.0f, 1.0f, 0.0f);
    glScalef(1.0f / zoomLevel, 1.0f / zoomLevel, 1.0f / zoomLevel);

    // draw the object itself
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < mesh.vertices.size(); i++) {
        QVector3D vertex = mesh.vertices[i];
        glColor3f(1.0f * (i % 3 == 0), 1.0f * (i % 3 == 1), 1.0f * (i % 3 == 2));
        glVertex3f(vertex.x(), vertex.y(), vertex.z());
    }

    glEnd();
    glDisable(GL_DEPTH_TEST);
    update();
}

/*
*   Function: load the ObjectMesh object into the Visualizer module 
*   Return: no return values
*   TODO: a quad having object could set a flag for that in this function
*/
void Visualizer::loadMesh(ObjectMesh object) {
    mesh = object;
}

/*
*   Function: make sure a value of angles is within the 360 degree range
*   Return: a float value that corresponds to the rotation angle without excessive 360 rotations
*/
float loopAround(float val, float min, float max)
{
    if (val < min)
    {
        return max - fmod(min - val, max - min);
    }
    else if (val >= max)
    {
        return min + fmod(val - max, max - min);
    }
    return val;
}

/*
*   Function: handler of mouse press/drag/release/wheel events,
*       in essence the implementation of rotate/zoom/drag functions
*   Return: no return values
*/
void Visualizer::mousePressEvent(QMouseEvent* event)
{
    // left button: rotate the model
    if (event->button() == Qt::LeftButton) {
        isRotating = true;
        lastMousePos = event->pos();
    }
    // right button: drag and move the object around
    if (event->button() == Qt::RightButton) {
        isDragging = true;
        lastMousePos = event->pos();
    }
}
void Visualizer::mouseReleaseEvent(QMouseEvent* event)
{
    // left button: rotate the model
    if (event->button() == Qt::LeftButton) {
        isRotating = false;
    }
    // right button: drag and move the object around
    if (event->button() == Qt::RightButton) {
        isDragging = false;
    }
}
void Visualizer::mouseMoveEvent(QMouseEvent* event)
{
    if (isRotating) {
        QPoint currentPos = event->pos();

        // calculate the change in mouse position between previous and current positions
        float deltaX = currentPos.x() - lastMousePos.x();
        float deltaY = currentPos.y() - lastMousePos.y();

        // update the rotation angles based on the mouse movement
        // TODO: the 0.5 value could be a variable set by the user as a rotation speed 
        m_rotationAngleY += deltaX * 0.5f;
        m_rotationAngleX += deltaY * 0.5f;

        // make sure the rotation angles are within the 360 degrees range
        m_rotationAngleX = loopAround(m_rotationAngleX, 0.0f, 359.999f);
        m_rotationAngleY = loopAround(m_rotationAngleY, 0.0f, 359.999f);
      
        lastMousePos = currentPos;
        update();
    }
    if (isDragging) {
        QPoint currentPos = event->pos();
        QVector3D viewDirection;

        // get the view direction so that the object can be dragged in all axes not just in the X and Y direction
        viewDirection.setX(sin(qDegreesToRadians(m_rotationAngleY)) * cos(qDegreesToRadians(m_rotationAngleX)));
        viewDirection.setY(-sin(qDegreesToRadians(m_rotationAngleX)));
        viewDirection.setZ(-cos(qDegreesToRadians(m_rotationAngleY)) * cos(qDegreesToRadians(m_rotationAngleX)));

        // QVector3D movementDirection = viewDirection.normalized();
        // calculate the change in mouse position
        float deltaX = currentPos.x() - lastMousePos.x();
        float deltaY = currentPos.y() - lastMousePos.y();

        // Update the drag values based on the mouse movement
        // TODO: currently doing two divisions by zoomLevel to get a good dragging speed, could use a variable and have the user set it
        translateX += deltaX * 0.05f / zoomLevel / zoomLevel; //* movementDirection.x();
        translateY -= deltaY * 0.05f / zoomLevel / zoomLevel;//* movementDirection.y();
        translateZ += deltaX * 0.05f / zoomLevel / zoomLevel;//* movementDirection.z();

        lastMousePos = currentPos;
        update();
    }
}
void Visualizer::wheelEvent(QWheelEvent* event)
{
    // zoom in/out based on the scroll wheel delta
    QPoint angleDelta = event->angleDelta();
    if (!angleDelta.isNull()) {
        //qDebug() << zoomLevel;
        float zoomFactor = 0.3f;
        if (angleDelta.y() > 0) {
            zoomLevel -= zoomFactor;
        }
        else {
            zoomLevel += zoomFactor;
        }
        isZooming = true;
        update();
    }
}

