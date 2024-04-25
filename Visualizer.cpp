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
Visualizer::Visualizer(QWidget* parent)
    : QOpenGLWidget(parent)
{
    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
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
    // Set the clear color to blue
    glClearColor(0.0f, 0.0f, 1.0f, 1.0f);


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
    //shaderProgram.addShaderFromSourceFile(QOpenGLShader::Geometry, "geometry_shader.glsl"); //unneeded?
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
    projectionMatrix.perspective(45.0f, aspectRatio, 0.1f, 100.0f);
}

/*
* Function: build an mvp matrix from the different operations applied to the mesh, such as rotations..
*/
void Visualizer::buildmvpMatrix() {
    // Remove the original changes to the mvp matrix by setting each matrix to identity
    projectionMatrix.setToIdentity();
    viewMatrix.setToIdentity();
    modelMatrix.setToIdentity();

    // Apply the different operations
    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    projectionMatrix.perspective(45, aspectRatio, 0.1f, 100.0f);

    viewMatrix.translate(0, 0, translateZ);
    viewMatrix.rotate(m_rotationAngleX, 1.0f, 0.0f, 0.0f);
    viewMatrix.rotate(m_rotationAngleY, 0.0f, 1.0f, 0.0f);

    modelMatrix.translate(translateX, translateY, 0.0f);
    modelMatrix.scale(1.0f / zoomLevel);


    // Combine the matrices to form mvp matrix
    mvpMatrix = projectionMatrix * viewMatrix * modelMatrix;
}

/*
* Function: helper function that only draws the triangle currently under the cursor.
*   done for color picking of vertex, each vertex of the triangle emits an RGB value,
*   the pixel under cursor will have an rgb value with the closest vertex having the
*   highest value, which will be used to find it
*/
void Visualizer::paintPickedFace(int faceID) {
    // Set up
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    // Get coordinates of the face
    QVector3D vertex1 = mesh.vertices[mesh.indices[faceID * 3]];
    QVector3D vertex2 = mesh.vertices[mesh.indices[faceID * 3 + 1]];
    QVector3D vertex3 = mesh.vertices[mesh.indices[faceID * 3 + 2]];

    vertex1 = mvpMatrix.map(vertex1);
    vertex2 = mvpMatrix.map(vertex2);
    vertex3 = mvpMatrix.map(vertex3);

    // Draw the face in fixed mode, maybe change to programmable?
    glBegin(GL_TRIANGLES);
    glColor3f(1.0, 0.0, 0.0);
    glVertex3f(vertex1.x(), vertex1.y(), vertex1.z());
    glColor3f(0.0, 1.0, 0.0);
    glVertex3f(vertex2.x(), vertex2.y(), vertex2.z());
    glColor3f(0.0, 0.0, 1.0);
    glVertex3f(vertex3.x(), vertex3.y(), vertex3.z());
    glEnd();
    glFinish();
}

/*
*   Function: after the color picking mesh is drawn, this function gets the id of the face under the cursor,
*       also gets the closest vertex, both the face and vertex info is stored in Visualizer member variables
*/
void Visualizer::pickFromMesh() {
    unsigned char pixel[3];

    // Get pixel color under cursor and derive the ID from it
    glReadPixels(mouseX, mouseY, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, &pixel);
    GLuint id = (pixel[0] << 16) | (pixel[1] << 8) | pixel[2];

    // If the id is valid (the pixel value does not correspond to the white background)
    if (id >= 0 && id != 16777215 && addingConstraints)
    {
        paintPickedFace(id);
        glReadPixels(mouseX, mouseY, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, &pixel);

        // Get the maximum value in the RGB composition of the pixel, vertex1 is R, 2 is G, 3 is B...
        // The highest value corresponds to the closest vertex as it is linearly interpolated through the triangle
        int max = 0;
        int val = pixel[0];
        for (int i = 1; i < 3; i++) {
            if (pixel[i] > val) {
                max = i;
                val = pixel[max];
            }
        }
        chosenVertex = mesh.indices[id * 3 + max];
    }
    selectedFace = id;
}


/*
*   Function: paint the mesh with each face having a different color, which sets the screen up for the color picking function
*/
void Visualizer::paintPickFrame() {
    // Set up
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    // For each face, draw the triangle with a different color
    for (size_t i = 0; i < mesh.indices.size(); i += 3) {
        GLuint id = i / 3; // ID corresponding to the triangle
        int red = (id >> 16) & 0xFF;
        int green = (id >> 8) & 0xFF;
        int blue = id & 0xFF;

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

        // Paint the triangle
        glBegin(GL_TRIANGLES);
        glColor3f(r, g, b);
        glVertex3f(vertex1.x(), vertex1.y(), vertex1.z());
        glVertex3f(vertex2.x(), vertex2.y(), vertex2.z());
        glVertex3f(vertex3.x(), vertex3.y(), vertex3.z());
        glEnd();
    }
    glFinish();
}


/*
*   Function: paint the constraint set by the user, currently a thick line running from the selected constrained vertices
*/
void Visualizer::paintConstraintHighlight(QVector3D start, QVector3D end) {
    // Set up
    glLineWidth(5.0f);
    glColor3f(0.0f, 1.0f, 0.0f);

    glBegin(GL_LINES);
    glVertex3f(start.x(), start.y(), start.z());
    glVertex3f(end.x(), end.y(), end.z());

    glEnd();
}

/*
*   Function: paints all the constraints of different constraint sets iteratively
*/
void Visualizer::paintConstraints() {

    for (Constraint* c : constraints) {
        if (c->vertices.size() > 1) {
            for (int i = 0; i < c->vertices.size() - 1; i++) {
                // Get the 2 vertices of the constraint and perform the mvp matrix transformations on them
                QVector3D vertex1 = mesh.vertices[c->vertices[i]];
                vertex1 = mvpMatrix.map(vertex1);
                QVector3D vertex2 = mesh.vertices[c->vertices[i + 1]];
                vertex2 = mvpMatrix.map(vertex2);

                paintConstraintHighlight(vertex1, vertex2);
            }
        }
    }
}

/*
*   Function: allow the user to set constraints in the visualizer, in the future
*       will disable a different mode selected by the user previously
*/
void Visualizer::setConstraintsMode() {
    qDebug() << "Setting constraints mode...";
    addingConstraints = true;
    showConstraints = true;
}

/*
*   Function: paint the current 3D mesh with all transformations
*/
void Visualizer::paintMesh()
{
    // General OpenGL set up
    glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    shaderProgram.bind();
    shaderProgram.setUniformValue("mvpMatrix", mvpMatrix);

    // Shader lighting values computation set up
    QMatrix4x4 invViewMatrix = viewMatrix.inverted();
    QVector3D cameraPosition = invViewMatrix.column(3).toVector3D();
    QVector3D lightDirection = cameraPosition.normalized();
    //shaderProgram.setUniformValue("lightDirection", lightDirection);

    // Start the paint process
    glBindVertexArray(vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glFinish();

    // Need to check which face is selected, best done in the CPU rather than shader, as it will be used later on
    for (int i = 0; i < mesh.indices.size(); i += 3) {
        QVector3D v0 = mesh.vertices[mesh.indices[i]].normalized();
        QVector3D v1 = mesh.vertices[mesh.indices[i + 1]].normalized();
        QVector3D v2 = mesh.vertices[mesh.indices[i + 2]].normalized();

        QVector3D edge1 = (v1 - v0).normalized();
        QVector3D edge2 = (v2 - v0).normalized();
        QVector3D faceNormal = QVector3D::crossProduct(edge1, edge2).normalized();
        float shading = qMax(0.0f, QVector3D::dotProduct(faceNormal, lightDirection));

        int isSelected = 0;
        if (i / 3 == selectedFace) {
            isSelected = 1;
        }

        // Set the uniform values for the shader, selectedFace determines face color,
        // shadingValue determines the 'darkness' of the face about to be painted
        shaderProgram.setUniformValue("selectedFace", isSelected);
        shaderProgram.setUniformValue("shadingValue", shading);

        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, reinterpret_cast<void*>(sizeof(GLuint) * i));
    }
    shaderProgram.release();
    glBindVertexArray(0);
}

/*
*   Function: paint the OpenGL scene, including mesh, pick mesh, constraints etc.
*   Return: no return values
*   TODO: quad having objects may need to be drawn in a different method (GL_TRIANGLES vs GL_QUADS?)
*       currently the color is just a simple gradient, add shading and different texture
*/
void Visualizer::paintGL()
{
    // Set up
    buildmvpMatrix();
    glClearColor(0.0f, 0.0f, 0.5451f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    if (objectLoaded) {
        // Reset the chosen vertex so that it is not added to the constraints twice
        chosenVertex = -1;         // 0 for first in indexed face, 1 for 2nd etc..
        paintPickFrame();
        pickFromMesh();
        glFinish();
        paintMesh();
        paintConstraints();
    }
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
//TODO: more valid constraint checking
void Visualizer::addConstraintVertex() {
    if (chosenVertex != -1) {
        QVector3D vertex = mesh.vertices[chosenVertex];             //should constraints be vertices or indices?
        qDebug() << "Vertex: " << chosenVertex << ": " << vertex;

        int cnt = std::count(currentConstraint->vertices.begin(), currentConstraint->vertices.end(), chosenVertex);

        if (cnt > 0) {
            qDebug() << "Vertex already in constraint...";
        }
        currentConstraint->vertices.push_back(chosenVertex);
        qDebug() << "pushed back, size: " << currentConstraint->vertices.size();
    }
    else {
        qDebug() << "No vertex in selection...";
    }
}

std::vector<Constraint*> Visualizer::getConstraints()
{
    return constraints;
}

void Visualizer::setConstraints(std::vector<Constraint*> c)
{
    qDebug() << "size of inputted constraints" << c.size();
    constraints = c;
    currentConstraint = new Constraint();
    constraints.push_back(currentConstraint);
}

void Visualizer::mousePressEvent(QMouseEvent* event)
{
    QPoint currentPos = event->pos();
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

void Visualizer::pushConstraints() {
    currentConstraint = new Constraint();
    constraints.push_back(currentConstraint);
}

void Visualizer::keyPressEvent(QKeyEvent* event)
{
    qDebug() << "Key pressed: " << event->key();
    if (event->key() == Qt::Key_C && addingConstraints) {
        qDebug() << "Adding constraint vertex..";
        addConstraintVertex();
    }
    if (event->key() == Qt::Key_Return && addingConstraints) {
        qDebug() << "Pushing constraints..";
        pushConstraints();
    }
    if (event->modifiers() == Qt::ControlModifier && event->key() == Qt::Key_Z)
    {
        qDebug() << "Ctrl+Z pressed";
        // in the future will backtrack many options, for now only constraints
        if (currentConstraint->vertices.size() > 0) {
            currentConstraint->vertices.pop_back();
        }
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
    mouseX = event->pos().x() * devicePixelRatio();
    mouseY = (height() - event->pos().y() - 1) * devicePixelRatio();
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

