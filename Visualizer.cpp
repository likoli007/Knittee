#include "Visualizer.h"
#include <QWidget>
#include <QOpenGLWidget>
//#include <QOpenGLFunctions>
#include <QOpenGLExtraFunctions>
#include <GL/gl.h>
//#include <GL/glu.h>
//#include "ObjLoader.h"
#include <QDebug>
#include <QToolTip>
#include <QApplication>
#include <QScreen>
/*
*   Function: constructor, adds mouse tracking for the visualization subwindow
*/
Visualizer::Visualizer(QWidget* parent )
    : QOpenGLWidget(parent)
{
    setMouseTracking(true);
    //setFixedSize(800, 600); // Set the initial size
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

    glEnable(GL_DEPTH_TEST);
    glClearColor(0.0f, 0.0f, 1.0f, 1.0f); // Set the clear color to blue
    
   
    // Create and bind the VAO, VBO 
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    glGenBuffers(1, &vbo);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    // Set up the vertex buffer data with a null buffer
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);

    glGenBuffers(1, &indexBuffer);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    // Set up the index buffer data with a null buffer
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
    // Specify the format of the vertex data (position in this case)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), nullptr);
    glEnableVertexAttribArray(0);



    // Create and bind the  buffer
    glGenBuffers(1, &shadingBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, shadingBuffer);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), nullptr);
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    // Load the shader program
    shaderProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, "vertex_shader.glsl");
    shaderProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, "fragment_shader.glsl");
    //shaderProgram.addShaderFromSourceFile(QOpenGLShader::Geometry, "geometry_shader.glsl");
    shaderProgram.link();

    // Set the default matrices
    viewMatrix.setToIdentity();
    projectionMatrix.setToIdentity();
    modelMatrix.setToIdentity();
    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    projectionMatrix.perspective(45.0f, aspectRatio, 0.1f, 10000.0f); // fov?
}

/*
*   Function: inherited function that is used when the user resizes the application
*/
void Visualizer::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);

    projectionMatrix.setToIdentity();
    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    projectionMatrix.perspective(45.0f, aspectRatio, 0.1f, 10000.0f);
}

//not working yet...


template <typename T>
T clamp(const T& value, const T& min, const T& max) {
    return std::min(std::max(value, min), max);
}

void Visualizer::drawPickFrame() {
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
        GLuint id = i / 3; // ID corresponding to the triangle
        int red = (id >> 16) & 0xFF;
        int green = (id >> 8) & 0xFF;
        int blue = id & 0xFF;
        //qDebug() << "ID: " << id;
        QColor color(red, green, blue);

        QVector3D vertex1 = mesh.vertices[mesh.indices[i]];
        QVector3D vertex2 = mesh.vertices[mesh.indices[i + 1]];
        QVector3D vertex3 = mesh.vertices[mesh.indices[i + 2]];

        vertex1 = mvpMatrix.map(vertex1);
        vertex2 = mvpMatrix.map(vertex2);
        vertex3 = mvpMatrix.map(vertex3);


        float r = color.redF();
        float g = color.greenF();
        float b = color.blueF();


        glBegin(GL_TRIANGLES);
        glColor3f(r, g, b);
        glVertex3f(vertex1.x(), vertex1.y(), vertex1.z());
        glVertex3f(vertex2.x(), vertex2.y(), vertex2.z());
        glVertex3f(vertex3.x(), vertex3.y(), vertex3.z());
        glEnd();
    }
    glFinish();
    unsigned char pixel[3];
    glReadPixels(mouseX, mouseY, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, &pixel);


    //qDebug() << mouseX << mouseY;

    GLuint id = (pixel[0] << 16) | (pixel[1] << 8) | pixel[2];
    qDebug() << "within ID: " << id;
    selectedFace = id;
    glFinish();
}


/*
*   Function: paint the OpenGL scene, that is the 3D object selected by the user 
*   Return: no return values
*   TODO: quad having objects may need to be drawn in a different method (GL_TRIANGLES vs GL_QUADS?)
*       currently the color is just a simple gradient, add shading and different texture
*/
void Visualizer::paintGL()
{
    projectionMatrix.setToIdentity();
    viewMatrix.setToIdentity();
    modelMatrix.setToIdentity();

    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    projectionMatrix.perspective(45, aspectRatio, 0.1f, 10000.0f);

    viewMatrix.translate(0, 0, translateZ);
    viewMatrix.rotate(m_rotationAngleX, 1.0f, 0.0f, 0.0f);
    viewMatrix.rotate(m_rotationAngleY, 0.0f, 1.0f, 0.0f);
    
    modelMatrix.translate(translateX, translateY, 0.0f);
    modelMatrix.scale(1.0f / zoomLevel);

    //combine the matrices to form mvp matrix
    mvpMatrix = projectionMatrix * viewMatrix * modelMatrix;


    QMatrix4x4 invViewMatrix = viewMatrix.inverted();
    QVector3D cameraPosition = invViewMatrix.column(3).toVector3D();
    //qDebug() << cameraPosition;


    glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    if (objectLoaded) {
        drawPickFrame();
        glFinish();
        glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        

        shaderProgram.bind();
        shaderProgram.setUniformValue("mvpMatrix", mvpMatrix);

        
        //QVector3D cameraPosition = viewMatrix.inverted().column(3).toVector3D();
        QVector3D lightDirection = cameraPosition.normalized();
        shaderProgram.setUniformValue("lightDirection", lightDirection);

        glBindVertexArray(vao);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
        glFinish();
        for (int i = 0; i < mesh.indices.size(); i += 3) {
            
            QVector3D v0 = mesh.vertices[mesh.indices[i]].normalized();
            QVector3D v1 = mesh.vertices[mesh.indices[i + 1]].normalized();
            QVector3D v2 = mesh.vertices[mesh.indices[i + 2]].normalized();

            QVector3D edge1 = (v1 - v0).normalized();
            QVector3D edge2 = (v2 - v0).normalized();
            QVector3D faceNormal = QVector3D::crossProduct(edge1, edge2).normalized();
            float shading = qMax(0.0f, QVector3D::dotProduct(faceNormal, lightDirection ));

            int isSelected = 0;
            if (i / 3 == selectedFace) {
                isSelected = 1;


			}

            shaderProgram.setUniformValue("selectedFace", isSelected);

            //qDebug() << shading;
            shaderProgram.setUniformValue("shadingValue", shading);
            glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, reinterpret_cast<void*>(sizeof(GLuint) * i));
        }
    }
    shaderProgram.release();
    glBindVertexArray(0);
    //update();
}


void Visualizer::computeRayFromMouse(QPoint currentPosition) {


    viewMatrix.setToIdentity();
    viewMatrix.translate(0,0,translateZ);
    projectionMatrix.setToIdentity();
    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    projectionMatrix.perspective(45.0f, aspectRatio, 0.1f, 10000.0f);

    viewMatrix.translate(translateX, translateY, translateZ);
    viewMatrix.rotate(m_rotationAngleX, 1.0f, 0.0f, 0.0f);
    viewMatrix.rotate(m_rotationAngleY, 0.0f, 1.0f, 0.0f);


    float ndcX = (2.0f * static_cast<float>(currentPosition.x())) / width() - 1.0f;
    float ndcY = 1.0f - (2.0f * static_cast<float>(currentPosition.y())) / height();

    //qDebug() << ndcX << ndcY;

    // Create near and far points in NDC
    QVector4D rayNear(ndcX, ndcY, -1.0f, 1.0f);
    QVector4D rayFar(ndcX, ndcY, 1.0f, 1.0f);

    // Convert the points to world space


    QMatrix4x4 invProjectionMatrix = projectionMatrix.inverted();
    QVector4D rayEyeNear = invProjectionMatrix * rayNear;
    QVector4D rayEyeFar = invProjectionMatrix * rayFar;
    rayEyeNear = rayEyeNear.normalized();
    rayEyeFar = rayEyeFar.normalized();

    //qDebug() << rayEyeFar << rayEyeNear;

    QMatrix4x4 invViewMatrix = viewMatrix.inverted();
    QVector3D rayNearWorld = QVector3D(invViewMatrix * rayEyeNear);
    QVector3D rayFarWorld = QVector3D(invViewMatrix * rayEyeFar);

    // Calculate the ray direction
    rayDir = (rayFarWorld - rayNearWorld).normalized();
    //qDebug() << rayNear << rayFar << rayDir;
    update();
}
/*
*   Function: load the ObjectMesh object into the Visualizer module 
*   Return: no return values
*   TODO: a quad having object could set a flag for that in this function
*/
void Visualizer::loadMesh(ObjectMesh object) {
    mesh = object;

    qDebug() << "start load " << mesh.vertices.size() << mesh.indices.size();
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    // Set up the vertex buffer data
    glBufferData(GL_ARRAY_BUFFER, mesh.vertices.size() * sizeof(QVector3D), mesh.vertices.data(), GL_STATIC_DRAW);

    // Create and bind the index buffer
    
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.indices.size() * sizeof(GLuint), mesh.indices.data(), GL_STATIC_DRAW);
    GLenum error;
    while ((error = glGetError()) != GL_NO_ERROR) {
        qDebug() << "OpenGL error: " << error;
    }
    // Specify the format of the vertex data (position only)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), nullptr);
    glEnableVertexAttribArray(0);

    // Unbind the VAO
    glBindVertexArray(0);

    objectLoaded = true;
    qDebug() << "end load";
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
    QPoint currentPos = event->pos();
    // left button: rotate the model
    if (event->button() == Qt::LeftButton) {
        isRotating = true;
        lastMousePos = event->pos();


        mouseX = event->pos().x() * devicePixelRatio();
        mouseY = (height() - event->pos().y() - 1) * devicePixelRatio();

        update();

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
    QPoint currentPos = event->pos();

    //QString text;
   // text = QString("%1 X %2").arg(event->pos().x()).arg(event->pos().y());
    // Global coordinate-system, but widget coord. position
    // QToolTip::showText(event->pos(), text); 
    // this should fix it
    //QToolTip::showText(this->mapToGlobal(event->pos()), text);


    if (isRotating) {
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
        // TODO: currently doing divisions by zoomLevel to get a good dragging speed, could use a variable and have the user set it
        translateX += deltaX * 0.05f / zoomLevel; //* movementDirection.x();
        translateY -= deltaY * 0.05f / zoomLevel;//* movementDirection.y();

        lastMousePos = currentPos;
        update();
    }
    update();
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

