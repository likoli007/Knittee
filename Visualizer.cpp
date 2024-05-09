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
#include <QPair>
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


void Visualizer::sheetChanged(std::vector<std::vector<FlatPoint>>* sheet) {
	qDebug() << "sheet changed";
	this->sheet = *sheet;
	//update();
}

QColor timeColor(float time) {
    constexpr int size = 3;
    static const QColor grad[size] = {
        QColor(51, 51, 255),   // Blue
        QColor(204, 204, 204), // Gray
        QColor(204, 51, 51)    // Red
    };

    time *= (size - 1);
    int i = qBound(0, qFloor(time), size - 2);
    qreal f = qBound(0.0, time - i, 1.0);

    float r = grad[i].redF() * (1 - f) + grad[i + 1].redF() * f;
    float g = grad[i].greenF() * (1 - f) + grad[i + 1].greenF() * f;
    float b = grad[i].blueF() * (1 - f) + grad[i + 1].blueF() * f;

    int red = static_cast<int>(r * 255.0f);
    int green = static_cast<int>(g * 255.0f);
    int blue = static_cast<int>(b * 255.0f);

    // Interpolate between colors
    QColor color(red, green, blue);

    return color;
}
const uint TimeTexSize = 16;

void Visualizer::generateTimeTexture() {
    timeTexImage = QImage(TimeTexSize, 1, QImage::Format_RGB888);
    for (int i = 0; i < TimeTexSize; ++i) {
        QColor color = timeColor(i / float(TimeTexSize - 1));
        //QColor qColor = QColor::fromRgbF(color.r, color.g, color.b);
        timeTexImage.setPixel(i, 0, color.rgb());
    }

    timeTexture = new QOpenGLTexture(timeTexImage);

    timeTexture->setWrapMode(QOpenGLTexture::ClampToEdge);
    timeTexture->setMinificationFilter(QOpenGLTexture::Nearest);
    timeTexture->setMagnificationFilter(QOpenGLTexture::Linear);

    timeTexImage.save("stuffensie.jpg");

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
    glGenBuffers(1, &coordsBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, coordsBuffer);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(float)*2, nullptr);
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    // Load the shader program
    shaderProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, "vertex_shader.glsl");
    shaderProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, "fragment_shader.glsl");
    //shaderProgram.addShaderFromSourceFile(QOpenGLShader::Geometry, "geometry_shader.glsl"); //unneeded?
    shaderProgram.link();

    timeShaderProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, "time_vertex_shader.glsl");
    timeShaderProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, "time_fragment_shader.glsl");
    timeShaderProgram.link();


    // Set the default matrices
    viewMatrix.setToIdentity();
    projectionMatrix.setToIdentity();
    modelMatrix.setToIdentity();
    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    projectionMatrix.perspective(45.0f, aspectRatio, 0.1f, 10000.0f); // fov?

    setTimer();
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

void Visualizer::firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
    std::vector< std::vector< Stitch > >* active_stitches,
    RowColGraph* graph) {

    qDebug() << "FIRST ACTIVE CHAINS ARRIVED AT VISUALIZER";




    //assign the passed pointers to class vars
    this->activeChains = *active_chains;
    this->activeStitches = *active_stitches;
    this->rowColGraph = *graph;


    showFirstChains = true;
    showInterpolated = true;
    showLinks = false;
    showSlice = false;
    showSliceChains = false;

}

void Visualizer::meshInterpolated(ObjectMesh mesh, std::vector<float> values) {
    qDebug() << "visualizer loaded interpolated mesh!";
    showModel = false;




    interpolatedMesh = mesh;
    interpolatedValues = values;
    qDebug() << "interpolatedvalues size: " << interpolatedValues.size();

    generateTimeTexture();

    float min = -1.0f;
    float max = 1.0f;

    for (auto v : interpolatedValues) {
        min = std::min(min, v);
        max = std::max(max, v);
    }


    std::vector< float > texcoords;
    for (auto v : interpolatedValues) {
        texcoords.emplace_back(
            (((v - min) / (max - min)) * (TimeTexSize - 1) + 0.5f) / float(TimeTexSize));
        texcoords.emplace_back(0.5f);
    }


    qDebug() << texcoords.size() << interpolatedMesh.vertices.size() << interpolatedMesh.indices.size();

    qDebug() << "original mesh tris: " << mesh.indices.size() / 3 << " interpolated mesh tris: " << interpolatedMesh.indices.size() / 3;

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    // Set up the vertex buffer data
    glBufferData(GL_ARRAY_BUFFER, interpolatedMesh.vertices.size() * sizeof(QVector3D), interpolatedMesh.vertices.data(), GL_STATIC_DRAW);

    // Create and bind the index buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, interpolatedMesh.indices.size() * sizeof(GLuint), interpolatedMesh.indices.data(), GL_STATIC_DRAW);
    GLenum error;
    while ((error = glGetError()) != GL_NO_ERROR) {
        qDebug() << "OpenGL error: " << error;
    }
    // Specify the format of the vertex data (position only)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), nullptr);


    glBindBuffer(GL_ARRAY_BUFFER, coordsBuffer);
    glBufferData(GL_ARRAY_BUFFER, texcoords.size() * sizeof(float), texcoords.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glFinish();

    showInterpolated = true;
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
    //qDebug() << "Selected face: " << selectedFace;
    //qDebug() << "Chosen vertex: " << chosenVertex;
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
void Visualizer::paintPath(QVector3D start, QVector3D end, float r, float g, float b) {
    // Set up
    glLineWidth(5.0f);

    buildmvpMatrix();

    start = mvpMatrix.map(start);
    end = mvpMatrix.map(end);

    if (showConstraints) {
        if (r == 0.0f && b == 0.0f) {
            glColor3f(1.0f, 1.0f, 1.0f);
        }
        else {
            glColor3f(r, g, b);
        }
    }
    else {
        glColor3f(r, g, b);   
    }

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
            //select a color of the constraint based on the time value
            QColor color = timeColor(c->timeValue);
            int red = color.red();
            int green = color.green();
            int blue = color.blue();


            float r = static_cast<float>(red) / 255.0f;
            float g = static_cast<float>(green) / 255.0f;
            float b = static_cast<float>(blue) / 255.0f;

            
            for (int i = 0; i < c->vertices.size() - 1; i++) {
                //a constraint is a set of vertices, paint a line between each pair of vertices
                QVector3D vertex1 = mesh.vertices[c->vertices[i]];
                QVector3D vertex2 = mesh.vertices[c->vertices[i + 1]];

                paintPath(vertex1, vertex2, r, g, b);
            }
        }
    }
}

/*
*   Function: allow the user to set constraints in the visualizer, in the future
*       will disable a different mode selected by the user previously
*/
void Visualizer::setConstraintsMode(bool type) {
    addingConstraints = type;
    showConstraints = type;

    //if constraints mode is enabled, create a new constraint set
    if (type == true) {
        pushConstraints();
    }

    //if constraints mode is disabled, save them through Knittee
    if (type == false) {
        
        //traverse the constraint list and delete the ones that are empty or singular
        for (int i = 0; i < constraints.size(); i++) {
            if (constraints[i]->vertices.size() < 2) {
				constraints.erase(constraints.begin() + i);
				i--;
			}
		}

        emit requestConstraintsSave();
    }
}

/*
*   Function: paint the current 3D mesh with all transformations
*/
void Visualizer::paintOriginalMesh()
{
    // General OpenGL set up
    glClearColor(0.0f, 0.0f, 0.5451f, 1.0f);
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

void Visualizer::paintInterpolatedMesh() {
    glClearColor(0.0f, 0.0f, 0.5451f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

   
    //print out 2 texcord values per line in qdebug
    /*
    for (int i = 0; i < texcoords.size(); i += 2) {
		qDebug() << texcoords[i] << " " << texcoords[i + 1];
	}*/
    
    // For each face, draw the triangle with a different color
    glBindVertexArray(vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glFinish();

    //qDebug() << "painting interpolated mesh with (verts, tris) = " << interpolatedMesh.vertices.size() << ", " << interpolatedMesh.indices.size()/3;

   // qDebug() << "1";
    // Set the texture unit in the shader program
    timeShaderProgram.bind();
    timeTexture->bind(1);
    timeShaderProgram.setUniformValue("tex", 1);

    buildmvpMatrix();
    timeShaderProgram.setUniformValue("mvpMatrix", mvpMatrix);
    //qDebug() << "2";
    glDrawElements(GL_TRIANGLES, interpolatedMesh.indices.size(), GL_UNSIGNED_INT, nullptr);

    glFinish();
   // qDebug() << "3";
    timeShaderProgram.release();
    timeTexture->release();
    glBindVertexArray(0);

}

void Visualizer:: paintCube(QVector3D center, float sideLength) {
    float halfSide = sideLength;
    halfSide = 0.05f;
    buildmvpMatrix();
    center = mvpMatrix.map(center);
    
    qDebug() << sideLength;

    float x = center.x();
    float y = center.y();
    float z = center.z();
    
    glBegin(GL_QUADS);

    glColor3f(1.0f, 1.0f, 1.0f);
    // Front face
    glVertex3f(x - halfSide, y - halfSide, z + halfSide);
    glVertex3f(x + halfSide, y - halfSide, z + halfSide);
    glVertex3f(x + halfSide, y + halfSide, z + halfSide);
    glVertex3f(x - halfSide, y + halfSide, z + halfSide);

    // Back face
    glVertex3f(x - halfSide, y - halfSide, z - halfSide);
    glVertex3f(x - halfSide, y + halfSide, z - halfSide);
    glVertex3f(x + halfSide, y + halfSide, z - halfSide);
    glVertex3f(x + halfSide, y - halfSide, z - halfSide);

    // Top face
    glVertex3f(x - halfSide, y + halfSide, z - halfSide);
    glVertex3f(x - halfSide, y + halfSide, z + halfSide);
    glVertex3f(x + halfSide, y + halfSide, z + halfSide);
    glVertex3f(x + halfSide, y + halfSide, z - halfSide);

    // Bottom face
    glVertex3f(x - halfSide, y - halfSide, z - halfSide);
    glVertex3f(x + halfSide, y - halfSide, z - halfSide);
    glVertex3f(x + halfSide, y - halfSide, z + halfSide);
    glVertex3f(x - halfSide, y - halfSide, z + halfSide);

    // Right face
    glVertex3f(x + halfSide, y - halfSide, z - halfSide);
    glVertex3f(x + halfSide, y + halfSide, z - halfSide);
    glVertex3f(x + halfSide, y + halfSide, z + halfSide);
    glVertex3f(x + halfSide, y - halfSide, z + halfSide);

    // Left face
    glVertex3f(x - halfSide, y - halfSide, z - halfSide);
    glVertex3f(x - halfSide, y - halfSide, z + halfSide);
    glVertex3f(x - halfSide, y + halfSide, z + halfSide);
    glVertex3f(x - halfSide, y + halfSide, z - halfSide);

    glEnd();
}

void Visualizer::paintTraced() {
    std::vector<glm::vec3> colors = { glm::vec3(1.0f,0.0f,0.0f), glm::vec3(0.0f, 1.0f, 0.0f),
        glm::vec3(0.0f, 0.0f, 1.0f), glm::vec3(1.0f, 0.0f, 1.0f), glm::vec3(0.0f, 1.0f, 1.0f),
    glm::vec3(1.0f, 1.0f, 0.0f), glm::vec3(1.0f, 1.0f, 1.0f)};
    
    int yarnColor = 0;

    for (int i = 1; i < tracedMesh.size(); i++) {
        if (tracedMesh[i].yarn != tracedMesh[i - 1].yarn) {
            yarnColor = (yarnColor + 1) % colors.size();
            continue;
        }
        QVector3D start = QVector3D(tracedMesh[i - 1].at.x, tracedMesh[i - 1].at.y, tracedMesh[i - 1].at.z);
        QVector3D end = QVector3D(tracedMesh[i].at.x, tracedMesh[i].at.y, tracedMesh[i].at.z);
        paintPath(
            start, end,
            colors[yarnColor].x, colors[yarnColor].y, colors[yarnColor].z);
    }
    for (int i = 0; i < tracedMesh.size(); i++) {
        glm::vec3 const& at = tracedMesh[i].at;
        if (tracedMesh[i].ins[0] != -1U && tracedMesh[i].ins[1] != -1U) {
            glm::vec3 const& in0 = tracedMesh[tracedMesh[i].ins[0]].at;
            glm::vec3 const& in1 = tracedMesh[tracedMesh[i].ins[1]].at;
            
            QVector3D start = QVector3D(0.5f * (at.x + in0.x), 0.5f * (at.y + in0.y), 0.5f * (at.z + in0.z));
            QVector3D end = QVector3D(0.5f * (at.x + in1.x), 0.5f * (at.y + in1.y), 0.5f * (at.z + in1.z));

            paintPath(
				start, end,
				0.33f, 0.33f, 0.33f);
            paintPath(start, end, 0.8, 0.8, 0.8);
        }
        if (tracedMesh[i].outs[0] != -1U && tracedMesh[i].outs[1] != -1U) {
            glm::vec3 const& out0 = tracedMesh[tracedMesh[i].outs[0]].at;
            glm::vec3 const& out1 = tracedMesh[tracedMesh[i].outs[1]].at;

            QVector3D start = QVector3D(0.5f * (at.x + out0.x), 0.5f * (at.y + out0.y), 0.5f * (at.z + out0.z));
            QVector3D end = QVector3D(0.5f * (at.x + out1.x), 0.5f * (at.y + out1.y), 0.5f * (at.z + out1.z));

            paintPath(
                start, end,
                0.26f, 0.26f, 0.26f);
            paintPath(start, end, 0.73, 0.73, 0.73);
        }
        if (tracedMesh[i].ins[0] != -1U && tracedMesh[i].ins[1] == -1U) {
            glm::vec3 const& in0 = tracedMesh[tracedMesh[i].ins[0]].at;

            QVector3D start = QVector3D(0.5f * (at.x + in0.x), 0.5f * (at.y + in0.y), 0.5f * (at.z + in0.z));
            QVector3D end = QVector3D(at.x, at.y, at.z);

            paintPath(
				start, end,
				0.53f, 0.53f, 0.53f);

        }
        if (tracedMesh[i].outs[0] != -1U && tracedMesh[i].outs[1] == -1U) {
            glm::vec3 const& out0 = tracedMesh[tracedMesh[i].outs[0]].at;

            QVector3D start = QVector3D(0.5f * (at.x + out0.x), 0.5f * (at.y + out0.y), 0.5f * (at.z + out0.z));
            QVector3D end = QVector3D(at.x, at.y, at.z);

            paintPath(
				start, end,
				0.46f, 0.46f, 0.46f);
        }
    }
}

void Visualizer::knitGraphTraced(std::vector< TracedStitch >* t) {
    tracedMesh = *t;
    showTraced = true;
    showFirstChains = false;
    showInterpolated = false;
}

void Visualizer::paintFirstActiveChains() {
    std::vector< QVector3D > locations;
    locations.reserve(rowColGraph.vertices.size());
    for (auto const& v : rowColGraph.vertices) {
        locations.emplace_back(v.at.interpolate(interpolatedMesh.vertices));
    }

    
    //TODO: get parameters from knitgrapher
    float radius = 1.5f * 0.075f * 3.66f * 10.0f;
    for (uint32_t vi = 0; vi < rowColGraph.vertices.size(); ++vi) {
        auto const& v = rowColGraph.vertices[vi];


       /*
        paintSphere(
            locations[vi],
            radius
        );*/

        //paintCube(locations[vi], radius);

        if (v.row_in != -1U) {
            //qDebug() << "i have row_in";
            paintPath(
                locations[vi],
                locations[v.row_in],
                1.0f, 0.0f, 1.0f
            );
        }
        if (v.row_out != -1U) {
            //qDebug() << "i have row_out";
            paintPath(
                locations[vi],
                locations[v.row_out],
                1.0f, 0.0f, 1.0f
            );
        }
        if (v.col_in[0] != -1U) {
            paintPath(
                locations[vi],
                locations[v.col_in[0]],
                1.0f, 0.0f, 1.0f
            );
        }
        if (v.col_in[1] != -1U) {
            paintPath(
                locations[vi],
                locations[v.col_in[1]],
                1.0f, 0.0f, 1.0f
            );
        }

        if (v.col_out[0] != -1U) {
            paintPath(
                locations[vi],
                locations[v.col_out[0]],
                1.0f, 0.0f, 1.0f
            );
        }
        if (v.col_out[1] != -1U) {
            paintPath(
                locations[vi],
                locations[v.col_out[1]],
                1.0f, 0.0f, 1.0f
            );
        }
    }
}

void Visualizer::paintSphere(QVector3D center, float radius) {
    const int slices = 15;
    const int stacks = 15;

    radius = 0.04f / zoomLevel;
    //glClear(GL_DEPTH_BUFFER_BIT);
    
    //qDebug() << "drawing stuff";
    buildmvpMatrix();
    center = mvpMatrix.map(center);

    glBegin(GL_TRIANGLE_STRIP);
    
    //set white color
    glColor3f(0.0f, 0.57f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    for (int i = 0; i <= stacks; ++i)
    {
        float phi = M_PI * i / stacks;
        float sinPhi = sin(phi);
        float cosPhi = cos(phi);

        for (int j = 0; j <= slices; ++j)
        {
            float theta = 2 * M_PI * j / slices;
            float sinTheta = sin(theta);
            float cosTheta = cos(theta);

            float x = cosTheta * sinPhi;
            float y = cosPhi;
            float z = sinTheta * sinPhi;

            glVertex3f(center.x() + x * radius, center.y() + y * radius, center.z() + z * radius);
            x = cosTheta * sin(phi + M_PI / stacks);
            y = cos(phi + M_PI / stacks);
            z = sinTheta * sin(phi + M_PI / stacks);

            glVertex3f(center.x() + x * radius, center.y() + y * radius, center.z() + z * radius);
        }
    }
    glEnd();
}

void Visualizer::peelSliceDone(ObjectMesh* slice, std::vector< std::vector< uint32_t > >* slice_active_chains_, std::vector< std::vector< uint32_t > >* slice_next_chains_) {
    showInterpolated = false;
    qDebug() << "slice arrived at visualizer";
    sliceMesh = *slice;
    sliceActiveChains = *slice_active_chains_;
    sliceNextChains = *slice_next_chains_;



    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);

    // Set up the vertex buffer data
    glBufferData(GL_ARRAY_BUFFER, sliceMesh.vertices.size() * sizeof(QVector3D), sliceMesh.vertices.data(), GL_STATIC_DRAW);

    // Create and bind the index buffer
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sliceMesh.indices.size() * sizeof(GLuint), sliceMesh.indices.data(), GL_STATIC_DRAW);
    GLenum error;
    while ((error = glGetError()) != GL_NO_ERROR) {
        qDebug() << "OpenGL error: " << error;
    }
    // Specify the format of the vertex data (position only)
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), nullptr);


    glEnableVertexAttribArray(0);
    
    glFinish();

    showSlice = true;
    showSliceChains = true;
}


//paint both the slice mesh and the chains
void Visualizer::paintSliceMesh() {
    glClearColor(0.0f, 0.0f, 0.5451f, 1.0f);
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
    for (int i = 0; i < sliceMesh.indices.size(); i += 3) {
        QVector3D v0 = sliceMesh.vertices[sliceMesh.indices[i]].normalized();
        QVector3D v1 = sliceMesh.vertices[sliceMesh.indices[i + 1]].normalized();
        QVector3D v2 = sliceMesh.vertices[sliceMesh.indices[i + 2]].normalized();

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

    if (showSliceChains) {
        for (auto activeChain : sliceActiveChains) {
            for (int i = 0; i < activeChain.size() - 1; i++) {
                QVector3D vertex1 = sliceMesh.vertices[activeChain[i]];
                QVector3D vertex2 = sliceMesh.vertices[activeChain[i + 1]];
                paintPath(vertex1, vertex2, 1.0f, 0.0f, 0.0f);
            }
        }

        for (auto nextChain : sliceNextChains) {
            for (int i = 0; i < nextChain.size() - 1; i++) {
                QVector3D vertex1 = sliceMesh.vertices[nextChain[i]];
                QVector3D vertex2 = sliceMesh.vertices[nextChain[i + 1]];
                paintPath(vertex1, vertex2, 0.0f, 0.0f, 1.0f);
            }
        }
    }
}

//brazenly copied from knitgrapher, move to its own 'common functions' file later
QVector3D mmix(const QVector3D& wa, const QVector3D& wb, float m) {
    return wa + m * (wb - wa);
}

std::vector< std::vector< QVector3D > > copyLocations(ObjectMesh const& model, std::vector< std::vector< uint32_t > > const& chains) {
    std::vector< std::vector< QVector3D > > locations;
    locations.reserve(chains.size());
    for (auto const& chain : chains) {
        locations.emplace_back();
        locations.back().reserve(chain.size());
        for (auto v : chain) {
            locations.back().emplace_back(model.vertices[v]);
        }
    }
    return locations;
}
std::vector< std::vector< QVector3D > > interpolateStitchLocations(std::vector< std::vector< QVector3D > > const& chains, std::vector< std::vector< Stitch > > const& stitches) {
    assert(stitches.size() == chains.size());

    std::vector< std::vector< QVector3D > > locations;
    locations.reserve(locations.size());

    for (auto const& chain : chains) {
        locations.emplace_back();

        uint32_t ci = &chain - &chains[0];
        if (stitches[ci].empty()) continue;

        std::vector< float > lengths;
        lengths.reserve(chain.size());
        lengths.emplace_back(0.0f);
        for (uint32_t pi = 0; pi + 1 < chain.size(); ++pi) {
            QVector3D a = chain[pi];
            QVector3D b = chain[pi + 1];
            lengths.emplace_back(lengths.back() + (b - a).length());
        }
        assert(lengths.size() == chain.size());

        locations.back().reserve(stitches[ci].size());
        auto li = lengths.begin();
        for (auto const& s : stitches[ci]) {
            float l = s.t * lengths.back();
            while (li != lengths.end() && *li <= l) ++li;
            assert(li != lengths.end());
            assert(li != lengths.begin());
            float m = (l - *(li - 1)) / (*li - *(li - 1));
            uint32_t i = li - lengths.begin();
            locations.back().emplace_back(mmix(chain[i - 1], chain[i], m));
        }
    }

    return locations;
}

void Visualizer::nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >* next_chains) {
    nextActiveChains = *next_chains;

    showLinks = false;
    showSliceChains = false;
    showNextChains = true;
}

void Visualizer::paintLinks() {
    std::vector< std::vector< QVector3D > > from = interpolateStitchLocations(copyLocations(sliceMesh, sliceActiveChains), activeStitches);
    std::vector< std::vector< QVector3D > > to = interpolateStitchLocations(copyLocations(sliceMesh, sliceNextChains), nextStitches);



    for (auto link : links) {

        QVector3D const& vertex1 = from[link.from_chain][link.from_stitch];
        QVector3D const& vertex2 = to[link.to_chain][link.to_stitch];

		paintPath(vertex1, vertex2, 0.0f, 1.0f, 0.0f);
	}



}

void Visualizer::linkChainsDone(std::vector< std::vector< Stitch > >* next_stitches, std::vector< Link >* links_) {
    nextStitches = *next_stitches;
    links = *links_;

    qDebug() << "links arrived at visualizer";

    showLinks = true;


}

std::vector< std::vector< QVector3D > > interpolateLocations(ObjectMesh const& model, std::vector< std::vector< EmbeddedVertex > > const& chains) {
    std::vector< std::vector< QVector3D > > locations;
    locations.reserve(chains.size());
    for (auto const& chain : chains) {
        locations.emplace_back();
        locations.back().reserve(chain.size());
        for (auto const& ev : chain) {
            locations.back().emplace_back(ev.interpolate(model.vertices));
        }
    }
    return locations;
}

void Visualizer::paintNextChains() {

    std::vector< std::vector< QVector3D > > next_active_locations = interpolateLocations(interpolatedMesh, nextActiveChains);

    for (auto nextChain : next_active_locations) {
        for (int i = 0; i < nextChain.size() - 1; i++) {
            QVector3D vertex1 = nextChain[i];
            QVector3D vertex2 = nextChain[i + 1];
            paintPath(vertex1, vertex2, 1.0f, 0.0f, 0.0f);
        }
    }

}



void Visualizer::paintSheet() {

    glClearColor(0.7451f, 0.7451f, 0.7451f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // Set up orthogonal projection


    // +1 here so that loops are not on the edges
    int sheetWidth = sheet[0].size()+1;
    int sheetHeight = sheet.size()+1;

    int portWidth = this->width();
    int portHeight = this->height();
    
    int cellWidth = portWidth / sheetWidth;
    int cellHeight = portHeight / sheetHeight;
    // Calculate the offset to center the drawing
    int xOffset = (portWidth - (sheetWidth - 1) * cellWidth) / 2;
    int yOffset = (portHeight - (sheetHeight - 1) * cellHeight) / 2;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, portWidth, 0, portHeight, -1, 1);
    glMatrixMode(GL_MODELVIEW);



    
    //qDebug() << xOffset << yOffset;

    for (int i = 0; i < sheet.size(); i++) {
        for (int j = 0; j < sheet[i].size(); j++) {
            // Calculate the position of the rectangle
            int rectX = j * cellWidth + xOffset; // Left of the cell
            int rectY = i * cellHeight + yOffset; // Top of the cell

            // Draw the rectangle around the cell
            glColor3f(0.0f, 0.0f, 0.0f); // Black color
            glBegin(GL_LINE_LOOP);
            glVertex2i(rectX, rectY); // Top-left corner of the cell
            glVertex2i(rectX + cellWidth, rectY); // Top-right corner of the cell
            glVertex2i(rectX + cellWidth, rectY + cellHeight); // Bottom-right corner of the cell
            glVertex2i(rectX, rectY + cellHeight); // Bottom-left corner of the cell
            glEnd();

            // Calculate the position of the line
            int lineX = j * cellWidth + xOffset; // Center of the cell horizontally
            int lineY = i * cellHeight + yOffset + cellHeight / 6; // About 1/6th of the way up from the bottom of the cell
            int endX = lineX + cellWidth / 4; // End of the line
            int endY = lineY;
        

            paintYarn(lineX, lineY, endX, endY);
           
            lineX = endX;
            lineY = endY;

            endX = lineX + sheet[i][j].offset*cellWidth - cellWidth/8;
            endY = i * cellHeight + yOffset + cellHeight + 2*(cellHeight / 6); // About 1/6th of the way down from the top of the cell

            paintYarn(lineX, lineY, endX, endY);

            lineX = endX;
            lineY = endY;
            endX = lineX + cellWidth / 2 + cellWidth/8;
            endY = lineY;

            paintYarn(lineX, lineY, endX, endY);


            lineX = endX;
            lineY = endY;
            endX = lineX - sheet[i][j].offset * cellWidth - cellWidth / 8;
            endY =  i * cellHeight + yOffset + cellHeight / 6;

            paintYarn(lineX, lineY, endX, endY);

            lineX = endX;
            lineY = endY;
            endX = lineX+cellWidth/2;
            endY = lineY;
            paintYarn(lineX, lineY, endX, endY);

            glLineWidth(1.0f); // Reset line width to default
        }
    }

}

void Visualizer::paintYarn(int ax, int ay, int bx, int by) {

    glColor3f(1.0f, 0.0f, 0.0f); // Red color
    glLineWidth(5.0f); // Set line width to make it thicker
    glBegin(GL_LINES);
    glVertex2i(ax, ay); // Start from the left side of the rectangle
    glVertex2i(bx, by); // End at the right side of the rectangle
    glEnd();

    // Draw the thinner line (black) on top for the outline effect
    glColor3f(0.0f, 0.0f, 0.0f); // Black color
    glLineWidth(10.0f); // Reset line width to default
    glBegin(GL_LINES);
    glVertex2i(ax, ay); // Start from the left side of the rectangle
    glVertex2i(bx, by); // End at the right side of the rectangle
    glEnd();
}

/*
*   Function: paint the OpenGL scene, including mesh, pick mesh, constraints etc.
*   Return: no return values
*   TODO: quad having objects may need to be drawn in a different method (GL_TRIANGLES vs GL_QUADS?)
*       currently the color is just a simple gradient, add shading and different texture
*/
void Visualizer::paintGL()
{
    if (projectType == 0) {
        // Set up
        buildmvpMatrix();
        glClearColor(0.0f, 0.0f, 0.5451f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);

        if (showModel) {
            // Reset the chosen vertex so that it is not added to the constraints twice

            //WILL IT BREAK? WHO KNOWS? I DON'T
            //chosenVertex = -1;         // 0 for first in indexed face, 1 for 2nd etc..
            paintPickFrame();
            pickFromMesh();
            glFinish();
            paintOriginalMesh();
            paintConstraints();
        }
        if (showInterpolated) {
            //qDebug() << "drawing interpolatedmesh";
            paintInterpolatedMesh();
        }
        if (showFirstChains) {
            //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            //glEnable(GL_DEPTH_TEST);
            //paintInterpolatedMesh();
            //qDebug() << "drawing first active chains";
            paintFirstActiveChains();
        }
        if (showSlice) {
            paintSliceMesh();
        }
        if (showLinks) {
            paintLinks();
        }
        if (showNextChains) {
            paintNextChains();
        }
        if (showTraced) {
			paintTraced();
		}
    }
    else {
        paintSheet();
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

    showModel = true;
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


void Visualizer::pushConstraints() {
    currentConstraint = new Constraint();
    constraints.push_back(currentConstraint);
}

int Visualizer::getClosestConstraint() {
    if (chosenVertex != -1) {
        for (int i = 0; i < constraints.size(); i++) {
            if (constraints[i]->vertices.size() > 0) {
                for (int j = 0; j < constraints[i]->vertices.size(); j++) {
                    if (constraints[i]->vertices[j] == chosenVertex) {
                        return i;
                    }
                }
            }
		}
	}
	return -1;
}

void Visualizer::increaseConstraintValue() {

    qDebug() << "getting closest constraint, chosenVertex == " << chosenVertex;
    int closestConstraint = getClosestConstraint();

    qDebug() << "Closest constraint: " << closestConstraint;
    if (closestConstraint != -1 && constraints[closestConstraint]->timeValue < 1.0) {
		constraints[closestConstraint]->timeValue += 0.1f;
	}
}

void Visualizer::decreaseConstraintValue() {
	int closestConstraint = getClosestConstraint();
    if (closestConstraint != -1 && constraints[closestConstraint]->timeValue > -1.0) {
		constraints[closestConstraint]->timeValue -= 0.1f;
	}
}

void Visualizer::deleteConstraint() {
    int closestConstraint = getClosestConstraint();
    if (closestConstraint != -1) {
		constraints.erase(constraints.begin() + closestConstraint);
    }
}

void Visualizer::keyPressEvent(QKeyEvent* event)
{
    //if 'c', user wants to add a vertex
    if (event->key() == Qt::Key_C && addingConstraints) {
        qDebug() << "Adding constraint vertex..";
        addConstraintVertex();
    }
    //if 'enter', user wants to finish adding constraint vertices
    if (event->key() == Qt::Key_Return && addingConstraints) {
        qDebug() << "Pushing constraints..";
        pushConstraints();
    }
    //if '+', user wants to increase time value
    if (event->key() == Qt::Key_Plus ) {
		qDebug() << "Increasing constraint value..";
		increaseConstraintValue();
	}
    //if '-', user wants to delete time value
    if (event->key() == Qt::Key_Minus) {
        decreaseConstraintValue();
    }
    //if 'delete', user wants to delete the constraint
    if (event->key() == Qt::Key_Delete) {
        deleteConstraint();
    }

    //if 'ctrl-z', user wants to delete the last constraint vertex
    if (event->modifiers() == Qt::ControlModifier && event->key() == Qt::Key_Z)
    {
        qDebug() << "Ctrl+Z pressed";
        // in the future will backtrack many options, for now only constraints
        if (currentConstraint->vertices.size() > 0) {
            currentConstraint->vertices.pop_back();
        }
    }
}


QPair<int, int> Visualizer::getLoopPos(int mX, int mY) {
    QPair<int, int> res(-1, -1);

    int sheetWidth = sheet[0].size() + 1;
    int sheetHeight = sheet.size() + 1;

    int portWidth = this->width();
    int portHeight = this->height();

    int cellWidth = portWidth / sheetWidth;
    int cellHeight = portHeight / sheetHeight;
    // Calculate the offset to center the drawing
    int xOffset = (portWidth - (sheetWidth - 1) * cellWidth) / 2;
    int yOffset = (portHeight - (sheetHeight - 1) * cellHeight) / 2;
    for (int i = 0; i < sheet.size(); i++) {
        for (int j = 0; j < sheet[i].size(); j++) {
            int rectX = j * cellWidth + xOffset; // Left of the cell
            int rectY = i * cellHeight + yOffset; // Top of the cell
            int rectWidth = cellWidth;
            int rectHeight = cellHeight;

            if (mX > rectX && mX < rectX + rectWidth && mY > rectY && mY < rectY + rectHeight) {
                qDebug() << "Mouse is in cell: " << i << ", " << j;
                res.first = i;
                res.second = j;
                //chosenVertex = i * sheet[i].size() + j;
                //qDebug() << "Chosen vertex: " << chosenVertex;
            }
        }
    }
    return res;
}

void Visualizer::mousePressEvent(QMouseEvent* event)
{
    QPoint currentPos = event->pos();
    // left button: rotate the model
    if (event->button() == Qt::LeftButton) {
        if (projectType == 0) isRotating = true;
        else {
            isMovingLoop = true;
            mouseX = event->pos().x();
            mouseY = (height() - event->pos().y() - 1);
            from = getLoopPos(mouseX, mouseY);

        }
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
        isMovingLoop = false;
        to.first = -1;
        to.second = -1;
        from.first = -1;
        from.second = -1;
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
    if (isMovingLoop) {
        mouseX = event->pos().x();
        mouseY = (height() - event->pos().y() - 1);
        
        to = getLoopPos(mouseX, mouseY);

        if (from.first != -1 && to.first != -1) {
            sheet[from.first][from.second].offset = to.second - from.second;
        }

    }

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
        //update();
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
        //update();
    }
    //update();
}
void Visualizer::wheelEvent(QWheelEvent* event)
{
    // zoom in/out based on the scroll wheel delta
    QPoint angleDelta = event->angleDelta();
    if (!angleDelta.isNull()) {
        //qDebug() << zoomLevel;
        float zoomFactor = 0.03f;
        if (angleDelta.y() > 0) {
            zoomLevel -= zoomFactor;
        }
        else {
            zoomLevel += zoomFactor;
        }
        isZooming = true;
        //update();
    }
}

void Visualizer::stopTimer() {
    killTimer(timerID);
    timerID = 0;
}

void Visualizer::setTimer()
{
	timerID = startTimer(16);  //16ms should be 60fps

}

void Visualizer::timerEvent(QTimerEvent* event) 
{
    if (event->timerId() == timerID)
    {
        // Trigger a repaint event to update the scene
        update();
    }
}

void Visualizer::showInterpolatedChanged(int state) {
    if (state == Qt::Checked) {
        showInterpolated = true;
    }
    else if (state == Qt::Unchecked) {
        showInterpolated = false;
    }
}

void Visualizer::showGraphChanged(int state) {
    if (state == Qt::Checked) {
		showGraph = true;
	}
    else if (state == Qt::Unchecked) {
		showGraph = false;
	}
}

void Visualizer::showTracedChanged(int state) {
    if (state == Qt::Checked) {
		showTraced = true;
	}
    else if (state == Qt::Unchecked) {
		showTraced = false;
	}
}

