#include "Visualizer.h"

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
*   Function: get suitable color for a corresponding time value, creates the nice gradient for interpolated mesh
*/
QColor Visualizer::timeColor(float time) 
{
    constexpr int size = 3;
    static const QColor grad[size] = {
        QColor(51, 51, 255),   // Blue
        QColor(204, 204, 204), // Grey
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

    QColor color(red, green, blue);

    return color;
}

/*
*  Function: generate a time texture which will be used to pick the color of a vertex in the interpolated mesh
*/
void Visualizer::generateTimeTexture() 
{
    timeTexImage = QImage(TimeTexSize, 1, QImage::Format_RGB888);
    for (int i = 0; i < TimeTexSize; ++i) {
        QColor color = timeColor(i / float(TimeTexSize - 1));
        timeTexImage.setPixel(i, 0, color.rgb());
    }

    timeTexture = new QOpenGLTexture(timeTexImage);

    timeTexture->setWrapMode(QOpenGLTexture::ClampToEdge);
    timeTexture->setMinificationFilter(QOpenGLTexture::Nearest);
    timeTexture->setMagnificationFilter(QOpenGLTexture::Linear);

    //debug - show the texture
    //timeTexImage.save("stuffensie.jpg");
}

/*
*   Function: initialize the OpenGL graphics API, set some default parameters
*   Return: no return values
*/
void Visualizer::initializeGL()
{
    initializeOpenGLFunctions();

    glEnable(GL_DEPTH_TEST);

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

    // Create and bind the buffer
    glGenBuffers(1, &coordsBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, coordsBuffer);
    glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_STATIC_DRAW);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, sizeof(float)*2, nullptr);
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    // Load the shader program
    shaderProgram.addShaderFromSourceFile(QOpenGLShader::Vertex, "vertex_shader.glsl");
    shaderProgram.addShaderFromSourceFile(QOpenGLShader::Fragment, "fragment_shader.glsl");
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
void Visualizer::buildmvpMatrix() 
{
    // Remove the original changes to the mvp matrix by setting each matrix to identity
    projectionMatrix.setToIdentity();
    viewMatrix.setToIdentity();
    modelMatrix.setToIdentity();

    // Apply the different operations
    float aspectRatio = static_cast<float>(width()) / static_cast<float>(height());
    projectionMatrix.perspective(45, aspectRatio, 0.1f, 100.0f);

    viewMatrix.translate(-centroid.x(), -centroid.y(), -centroid.z());
    
    viewMatrix.translate(0.0, 0.0, -zoomLevel*zoomFactor);
    viewMatrix.translate(translateX, translateY, 0.0f);

    viewMatrix.rotate(m_rotationAngleX, 1.0f, 0.0f, 0.0f);
    viewMatrix.rotate(m_rotationAngleY, 0.0f, 1.0f, 0.0f);

    // Combine the matrices to form mvp matrix
    mvpMatrix = projectionMatrix * viewMatrix * modelMatrix;
}


/*
*  Function: change the visualized sheet as it was changed in the toolbar (usually the dimensions changed)
*/
void Visualizer::sheetChanged(std::vector<std::vector<FlatPoint>>* sheet)
{
    this->sheet = *sheet;
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
    RowColGraph* graph) 
{
    //assign the passed pointers to class vars
    this->activeChains = *active_chains;
    this->activeStitches = *active_stitches;
    this->rowColGraph = *graph;

    showFirstChains = true;
    showConstraints = false;
    showInterpolated = false;
    showLinks = false;
    showSlice = false;
    showSliceChains = false;
}

/*
*   Function: slot function that sets up the visualization of the remeshed and interpolated mesh
*/
void Visualizer::meshInterpolated(ObjectMesh mesh, std::vector<float> values) 
{
    showModel = false;
    interpolatedMesh = mesh;
    interpolatedValues = values;

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

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, interpolatedMesh.vertices.size() * sizeof(QVector3D), interpolatedMesh.vertices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, interpolatedMesh.indices.size() * sizeof(GLuint), interpolatedMesh.indices.data(), GL_STATIC_DRAW);
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
void Visualizer::pickFromMesh() 
{
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
void Visualizer::paintPickFrame() 
{
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
*   Function: paint the a path between two points, done for constraints, links, chains etc..
*/
void Visualizer::paintPath(QVector3D start, QVector3D end, float r, float g, float b, float radius = 5.0f) 
{
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(radius);

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

    glDisable(GL_LINE_SMOOTH);
    glEnd();
}

/*
*   Function: reset all parameters of the visualizer, used when the user wants to start a new project or resets the current one
*/
void Visualizer::reset() 
{
    interpolatedMesh.clear();
    interpolatedValues.clear();
    activeChains.clear();
    activeStitches.clear();
    rowColGraph.clear();
    showFirstChains = false;
    showInterpolated = false;
    showLinks = false;
    showSlice = false;
    showSliceChains = false;
    sliceMesh.clear();
    sliceActiveChains.clear();
    sliceNextChains.clear();
    nextStitches.clear();
    links.clear();
    nextActiveChains.clear();
    showNextChains = false;
    showGraph = false;
    showTraced = false;
    tracedMesh.clear();
    addingConstraints = false;
    showLastChain = false;
    showModel = true;
    showConstraints= true;

    translateX = 0.0f;
    translateY = 0.0f;
    translateZ = -5.0f;
    m_rotationAngleX = 0.0f;
    m_rotationAngleY = 0.0f;
    zoomLevel = 1.0f;

    loadMesh(mesh);
}

/*
*   Function: paints all the constraints of different constraint sets iteratively
*/
void Visualizer::paintConstraints() 
{
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
void Visualizer::setConstraintsMode(bool type) 
{
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
    glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    shaderProgram.bind();
    shaderProgram.setUniformValue("mvpMatrix", mvpMatrix);

    // Shader lighting values computation set up
    QMatrix4x4 invViewMatrix = viewMatrix.inverted();
    QVector3D cameraPosition = invViewMatrix.column(3).toVector3D();
    QVector3D lightDirection = cameraPosition.normalized();

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
        // in the future pass it by buffer?
        shaderProgram.setUniformValue("selectedFace", isSelected);
        shaderProgram.setUniformValue("shadingValue", shading);

        glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, reinterpret_cast<void*>(sizeof(GLuint) * i));
    }
    shaderProgram.release();
    glBindVertexArray(0);
}

/*
*  Function: load the interpolated mesh into the OpenGL buffers
*/
void Visualizer::loadInterpolated() 
{
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

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, interpolatedMesh.vertices.size() * sizeof(QVector3D), interpolatedMesh.vertices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, interpolatedMesh.indices.size() * sizeof(GLuint), interpolatedMesh.indices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), nullptr);
    glBindBuffer(GL_ARRAY_BUFFER, coordsBuffer);
    glBufferData(GL_ARRAY_BUFFER, texcoords.size() * sizeof(float), texcoords.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 0, nullptr);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glFinish();
}

/// the following are self-explanatory paint functions for the different visualizations
void Visualizer::paintInterpolatedMesh() 
{
    loadInterpolated();

    glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);

    // For each face, draw the triangle with a different color
    glBindVertexArray(vao);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glFinish();

    // Set the texture unit in the shader program
    timeShaderProgram.bind();
    timeTexture->bind(1);
    timeShaderProgram.setUniformValue("tex", 1);
    buildmvpMatrix();
    timeShaderProgram.setUniformValue("mvpMatrix", mvpMatrix);
    glDrawElements(GL_TRIANGLES, interpolatedMesh.indices.size(), GL_UNSIGNED_INT, nullptr);
    glFinish();
    timeShaderProgram.release();
    timeTexture->release();
    glBindVertexArray(0);

}
void Visualizer::paint3DCurve(QVector3D start, QVector3D end, QVector3D cp1, QVector3D cp2) 
{
    glColor3f(yarnMeshColor.x, yarnMeshColor.y, yarnMeshColor.z); 
  
    buildmvpMatrix();
    start = mvpMatrix.map(start);
	end = mvpMatrix.map(end);
    cp1 = mvpMatrix.map(cp1);
    cp2 = mvpMatrix.map(cp2);

    int numSegments = 16;

    //Draw a cubic Bezier curve
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i <= numSegments; ++i) {
        float t = static_cast<float>(i) / numSegments;
        float oneMinusT = 1.0f - t;

        //Bezier curve equation
        QVector3D point = oneMinusT * oneMinusT * oneMinusT * start +
            3 * oneMinusT * oneMinusT * t * cp1 +
            3 * oneMinusT * t * t * cp2 +
            t * t * t * end;

        glVertex3f(point.x(), point.y(), point.z());
    }
    glEnd();
}
void Visualizer::paintYarnMesh() 
{
    QVector3D cp1;
    QVector3D cp2;
    QVector3D direction;
    QVector3D currToTopDir;
    float radius = 1.0f;
    for (int i = 1; i < tracedMesh.size(); i++) {
        if (tracedMesh[i].yarn != tracedMesh[i - 1].yarn) {
            continue;
        }

        QVector3D previousPoint = QVector3D(tracedMesh[i - 1].at.x, tracedMesh[i - 1].at.y, tracedMesh[i - 1].at.z);
        QVector3D currentPoint = QVector3D(tracedMesh[i].at.x, tracedMesh[i].at.y, tracedMesh[i].at.z);

        // Calculate the direction vector from start to end
        direction = currentPoint - previousPoint;

        // Calculate the 2/3rds point along the direction vector
        QVector3D twoThirdsPoint = previousPoint + (2.0 / 3.0) * direction;

        QVector3D oneThirdPoint = previousPoint + (1.0 / 3.0) * direction;

        paintPath(
            oneThirdPoint, twoThirdsPoint,
            yarnMeshColor.x, yarnMeshColor.y, yarnMeshColor.z, radius);

        if (tracedMesh[i].outs[0] != -1U) {

            //painting yarn going 'up' to the stitch above
            QVector3D loopstart = currentPoint;
            QVector3D up = { tracedMesh[tracedMesh[i].outs[0]].at.x, tracedMesh[tracedMesh[i].outs[0]].at.y, tracedMesh[tracedMesh[i].outs[0]].at.z };
            QVector3D previousToUp = { tracedMesh[tracedMesh[i].outs[0] - 1].at.x, tracedMesh[tracedMesh[i].outs[0] - 1].at.y, tracedMesh[tracedMesh[i].outs[0] - 1].at.z };
            QVector3D nextToUp;
            QVector3D upperTwoThirdsPoint;
            if (tracedMesh[tracedMesh[i].outs[0] - 1].yarn == tracedMesh[tracedMesh[i].outs[0]].yarn) {

                direction = up - previousToUp;

                upperTwoThirdsPoint = previousToUp + (2.0 / 3.0) * direction;

                cp1 = currentPoint;
                cp2 = currentPoint;
                paint3DCurve(twoThirdsPoint, upperTwoThirdsPoint, cp1, cp2);
            }

            //if the previous nodes yarn is different, it necessarily implies you are at the start of a new yarn
            //therefore there must be at least one more point after right???
            else {
                nextToUp = { tracedMesh[(tracedMesh[i].outs[0]) + 1].at.x, tracedMesh[(tracedMesh[i].outs[0]) + 1].at.y, tracedMesh[(tracedMesh[i].outs[0]) + 1].at.z };
                direction = nextToUp - up;
                upperTwoThirdsPoint = up + (2.0 / 3.0) * direction;
                cp1 = currentPoint;
                cp2 = currentPoint;
                paint3DCurve(twoThirdsPoint, upperTwoThirdsPoint, cp1, cp2);
            }

            //painting the top 'loop cap' of a stitch
            QVector3D nextTwoThirdsPoint;   //last loop will have to do something special
            if (tracedMesh[i].outs[0] + 1 < tracedMesh.size()) {
                if (tracedMesh[tracedMesh[i].outs[0] + 1].yarn == tracedMesh[tracedMesh[i].outs[0]].yarn) {
                    nextToUp = { tracedMesh[(tracedMesh[i].outs[0]) + 1].at.x, tracedMesh[(tracedMesh[i].outs[0]) + 1].at.y, tracedMesh[(tracedMesh[i].outs[0]) + 1].at.z };
                    direction = nextToUp - up;
                    nextTwoThirdsPoint = up + (1.0 / 3.0) * direction;
                }
                else {
                    direction = up - previousToUp;
                    nextTwoThirdsPoint = up + (1.0 / 3.0) * direction;
                }
            }
            else {
                direction = up - previousToUp;
                nextTwoThirdsPoint = up + (1 / 3.0) * direction;
            }
            direction = up - currentPoint;
            QVector3D upperFourThirdsPoint = currentPoint + (5.0 / 3.0) * direction;

            cp1 = upperFourThirdsPoint;
            cp2 = upperFourThirdsPoint;

            paint3DCurve(upperTwoThirdsPoint, nextTwoThirdsPoint,cp1, cp2);

            //painting the 'down' yarn to the stitch below
            //going from nextTwoThirdsPoint to currentOneThirdPoint

            QVector3D start = nextTwoThirdsPoint;
            QVector3D nextToCurr = { tracedMesh[i+1].at.x, tracedMesh[i+1].at.y, tracedMesh[i+1].at.z };
            direction = nextToCurr - currentPoint;
            QVector3D nextOneThirdsPoint = currentPoint + (1.0 / 3.0) * direction;

            cp1 = currentPoint;
            cp2 = currentPoint;

            paint3DCurve(start, nextOneThirdsPoint, cp1, cp2);

        }
        
        if (tracedMesh[i].outs[1] != -1U) {
            //painting yarn going 'up' to the stitch above
            QVector3D loopstart = currentPoint;
            QVector3D up = { tracedMesh[tracedMesh[i].outs[1]].at.x, tracedMesh[tracedMesh[i].outs[1]].at.y, tracedMesh[tracedMesh[i].outs[1]].at.z };
            QVector3D previousToUp = { tracedMesh[tracedMesh[i].outs[1] - 1].at.x, tracedMesh[tracedMesh[i].outs[1] - 1].at.y, tracedMesh[tracedMesh[i].outs[1] - 1].at.z };

            direction = up - previousToUp;

            QVector3D upperTwoThirdsPoint = previousToUp + (2.0 / 3.0) * direction;

            cp1 = currentPoint;
            cp2 = currentPoint;
            paint3DCurve(twoThirdsPoint, upperTwoThirdsPoint, cp1, cp2);

            //painting the top 'loop cap' of a stitch
            QVector3D nextTwoThirdsPoint;   //last loop will have to do something special
            if (tracedMesh[i].outs[1] + 1 < tracedMesh.size()) {
                QVector3D nextToUp = { tracedMesh[(tracedMesh[i].outs[1]) + 1].at.x, tracedMesh[(tracedMesh[i].outs[1]) + 1].at.y, tracedMesh[(tracedMesh[i].outs[1]) + 1].at.z };
                direction = nextToUp - up;
                nextTwoThirdsPoint = up + (1.0 / 3.0) * direction;
            }
            else {
                nextTwoThirdsPoint = { 0,0,0 };
            }
            direction = up - currentPoint;
            QVector3D upperFourThirdsPoint = currentPoint + (5.0 / 3.0) * direction;

            cp1 = upperFourThirdsPoint;
            cp2 = upperFourThirdsPoint;

            paint3DCurve(upperTwoThirdsPoint, nextTwoThirdsPoint, cp1, cp2);

            //painting the 'down' yarn to the stitch below

            QVector3D start = nextTwoThirdsPoint;
            QVector3D nextToCurr = { tracedMesh[i + 1].at.x, tracedMesh[i + 1].at.y, tracedMesh[i + 1].at.z };
            direction = nextToCurr - currentPoint;
            QVector3D nextOneThirdsPoint = currentPoint + (1.0 / 3.0) * direction;

            cp1 = currentPoint;
            cp2 = currentPoint;

            paint3DCurve(start, nextOneThirdsPoint, cp1, cp2);
        }
        //else: it is the topmost row, draw according to extrapolated data from the previous row
        else if (tracedMesh[i].outs[1] == -1U && tracedMesh[i].outs[0] == -1U) {
            //paint the yarn going 'up' to the stitch above
            QVector3D current = QVector3D(tracedMesh[i].at.x, tracedMesh[i].at.y, tracedMesh[i].at.z);
            QVector3D previous = QVector3D(tracedMesh[i - 1].at.x, tracedMesh[i - 1].at.y, tracedMesh[i - 1].at.z);
            QVector3D oneDown = QVector3D(tracedMesh[tracedMesh[i].ins[0]].at.x, tracedMesh[tracedMesh[i].ins[0]].at.y, tracedMesh[tracedMesh[i].ins[0]].at.z);
            QVector3D oneDownPrev = QVector3D(tracedMesh[tracedMesh[i].ins[0] - 1].at.x, tracedMesh[tracedMesh[i].ins[0] - 1].at.y, tracedMesh[tracedMesh[i].ins[0] - 1].at.z);


            direction = current - oneDown;
            QVector3D extrapolatedNextPoint = oneDown + 2.0 * direction; 
            QVector3D extrapolatedTopPoint = oneDown + (2.5 * direction);

            direction = (oneDown - oneDownPrev);


            QVector3D extrapolatedPrevDirection = extrapolatedNextPoint + (1.0 / 3.0) * direction;
            QVector3D extrapolatedNextDirection = extrapolatedNextPoint - (1.0 / 3.0) * direction;

            //borrow the first ever yarn point as well as its prev and curr to get artificial stitch width

            cp1 = current;
            cp2 = current;

            paint3DCurve(twoThirdsPoint, extrapolatedNextDirection, cp1, cp2);
        
            //paint the top 'loop cap' of a stitch
            QVector3D start = extrapolatedNextPoint;
            cp1 = extrapolatedTopPoint;
            cp2 = extrapolatedTopPoint;
            paint3DCurve(extrapolatedNextDirection, extrapolatedPrevDirection, cp1, cp2);

            QVector3D nextOneThirdsPoint;
            cp1 = current;
            cp2 = current;
            //paint the 'down' yarn to the stitch below
            if (i < tracedMesh.size() - 1 && tracedMesh[i+1].yarn == tracedMesh[i].yarn) {
                QVector3D nextToCurr = { tracedMesh[i + 1].at.x, tracedMesh[i + 1].at.y, tracedMesh[i + 1].at.z };
                direction = nextToCurr - currentPoint;
                nextOneThirdsPoint = currentPoint + (1.0 / 3.0) * direction;
                
            }
            else {
                direction = (oneDown - oneDownPrev);
                nextOneThirdsPoint = currentPoint + (1.0 / 3.0) * direction;

            }
            paint3DCurve(extrapolatedPrevDirection, nextOneThirdsPoint, cp1, cp2);
        }
    }
}
void Visualizer::paintTraced() 
{
    int yarnColor = 0;

    for (int i = 1; i < tracedMesh.size(); i++) {
        if (tracedMesh[i].yarn != tracedMesh[i - 1].yarn) {
            yarnColor = (yarnColor + 1) % ycolors.size();
            continue;
        }
        QVector3D start = QVector3D(tracedMesh[i - 1].at.x, tracedMesh[i - 1].at.y, tracedMesh[i - 1].at.z);
        QVector3D end = QVector3D(tracedMesh[i].at.x, tracedMesh[i].at.y, tracedMesh[i].at.z);
        paintPath(start, end, ycolors[yarnColor].x, ycolors[yarnColor].y, ycolors[yarnColor].z);
    }
    for (int i = 0; i < tracedMesh.size(); i++) {
        glm::vec3 const& at = tracedMesh[i].at;
        if (tracedMesh[i].ins[0] != -1U && tracedMesh[i].ins[1] != -1U) {
            glm::vec3 const& in0 = tracedMesh[tracedMesh[i].ins[0]].at;
            glm::vec3 const& in1 = tracedMesh[tracedMesh[i].ins[1]].at;
            
            QVector3D start = QVector3D(0.5f * (at.x + in0.x), 0.5f * (at.y + in0.y), 0.5f * (at.z + in0.z));
            QVector3D end = QVector3D(0.5f * (at.x + in1.x), 0.5f * (at.y + in1.y), 0.5f * (at.z + in1.z));

            QVector3D atStart = QVector3D(at.x, at.y, at.z);

            paintPath(atStart, end, 0.33f, 0.33f, 0.33f);
            paintPath(atStart, start, 0.8, 0.8, 0.8);
        }
        if (tracedMesh[i].outs[0] != -1U && tracedMesh[i].outs[1] != -1U) {
            glm::vec3 const& out0 = tracedMesh[tracedMesh[i].outs[0]].at;
            glm::vec3 const& out1 = tracedMesh[tracedMesh[i].outs[1]].at;

            QVector3D start = QVector3D(0.5f * (at.x + out0.x), 0.5f * (at.y + out0.y), 0.5f * (at.z + out0.z));
            QVector3D end = QVector3D(0.5f * (at.x + out1.x), 0.5f * (at.y + out1.y), 0.5f * (at.z + out1.z));
            QVector3D atStart = QVector3D(at.x, at.y, at.z);

            paintPath(atStart, end, 0.26f, 0.26f, 0.26f);
            paintPath(atStart, start, 0.8, 0.8, 0.8);
        }
        if (tracedMesh[i].ins[0] != -1U && tracedMesh[i].ins[1] == -1U) {
            glm::vec3 const& in0 = tracedMesh[tracedMesh[i].ins[0]].at;

            QVector3D start = QVector3D(0.5f * (at.x + in0.x), 0.5f * (at.y + in0.y), 0.5f * (at.z + in0.z));
            QVector3D end = QVector3D(at.x, at.y, at.z);

            paintPath(start, end, 0.53f, 0.53f, 0.53f);
        }
        if (tracedMesh[i].outs[0] != -1U && tracedMesh[i].outs[1] == -1U) {
            glm::vec3 const& out0 = tracedMesh[tracedMesh[i].outs[0]].at;

            QVector3D start = QVector3D(0.5f * (at.x + out0.x), 0.5f * (at.y + out0.y), 0.5f * (at.z + out0.z));
            QVector3D end = QVector3D(at.x, at.y, at.z);

            paintPath(start, end, 0.46f, 0.46f, 0.46f);
        }
    }
}
void Visualizer::paintRowColGraph() 
{
    std::vector< QVector3D > locations;
    locations.reserve(rowColGraph.vertices.size());
    for (auto const& v : rowColGraph.vertices) {
        locations.emplace_back(v.at.interpolate(interpolatedMesh.vertices));
    }

    float r, g, b;
    float radius = 1.5f * 0.075f * 3.66f * 10.0f;
    for (uint32_t vi = 0; vi < rowColGraph.vertices.size(); ++vi) {
        auto& v = rowColGraph.vertices[vi];
        if (v.row_in != -1U) {
            if (showLastChain || (showNextChains && v.col_out[0] != -1U))
                paintPath(locations[vi],locations[v.row_in],
                    graphRowColor.x, graphRowColor.y, graphRowColor.z);
        }
        if (v.row_out != -1U) {
            if (!showNextChains || (showNextChains && v.col_out[0] != -1U))
                paintPath(locations[vi],locations[v.row_out],
                    graphRowColor.x, graphRowColor.y, graphRowColor.z);
        }
        if (v.col_in[0] != -1U) {
            r = v.col_in[1] == -1U ? graphLinkColor.x : graphCollapseColor.x;
            g = v.col_in[1] == -1U ? graphLinkColor.y : graphCollapseColor.y;
            b = v.col_in[1] == -1U ? graphLinkColor.z : graphCollapseColor.z;

            paintPath(locations[vi],locations[v.col_in[0]],r, g, b);
        }
        if (v.col_in[1] != -1U) {
            paintPath(locations[vi],locations[v.col_in[1]],
                graphCollapseColor.x, graphCollapseColor.y, graphCollapseColor.z);
        }

        if (v.col_out[0] != -1U && rowColGraph.vertices[v.col_out[0]].col_in[1] == -1) {
            r = v.col_out[1] == -1U ? graphLinkColor.x : graphExpandColor.x;
            g = v.col_out[1] == -1U ? graphLinkColor.y : graphExpandColor.y;
            b = v.col_out[1] == -1U ? graphLinkColor.z : graphExpandColor.z;

            paintPath(locations[vi],locations[v.col_out[0]], r, g, b);
        }
        if (v.col_out[1] != -1U && rowColGraph.vertices[v.col_out[1]].col_in[1] == -1) {
            paintPath(locations[vi], locations[v.col_out[1]],
                graphExpandColor.x, graphExpandColor.y, graphExpandColor.z);
        }
    }
}
void Visualizer::paintSliceMesh() 
{
    glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);
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
                paintPath(vertex1, vertex2, lowerBoundColor.x, lowerBoundColor.y, lowerBoundColor.z);
            }
        }

        for (auto nextChain : sliceNextChains) {
            for (int i = 0; i < nextChain.size() - 1; i++) {
                QVector3D vertex1 = sliceMesh.vertices[nextChain[i]];
                QVector3D vertex2 = sliceMesh.vertices[nextChain[i + 1]];
                paintPath(vertex1, vertex2, upperBoundColor.x, upperBoundColor.y, upperBoundColor.z);
            }
        }
    }
}
void Visualizer::paintLinks() 
{
    std::vector< std::vector< QVector3D > > from = interpolateStitchLocations(copyLocations(sliceMesh, sliceActiveChains), activeStitches);
    std::vector< std::vector< QVector3D > > to = interpolateStitchLocations(copyLocations(sliceMesh, sliceNextChains), nextStitches);

    for (auto link : links) {
        QVector3D const& vertex1 = from[link.from_chain][link.from_stitch];
        QVector3D const& vertex2 = to[link.to_chain][link.to_stitch];

        paintPath(vertex1, vertex2, graphLinkColor.x, graphLinkColor.y, graphLinkColor.z);
    }
}
void Visualizer::paintNextChains() 
{
    std::vector< std::vector< QVector3D > > next_active_locations = interpolateLocations(interpolatedMesh, nextActiveChains);

    for (auto nextChain : next_active_locations) {
        for (int i = 0; i < nextChain.size() - 1; i++) {
            QVector3D vertex1 = nextChain[i];
            QVector3D vertex2 = nextChain[i + 1];
            paintPath(vertex1, vertex2, lowerBoundColor.x, lowerBoundColor.y, lowerBoundColor.z);
        }
    }
}

void Visualizer::paintCurve(QPainter& painter, const std::vector<QPoint>& controlPoints) 
{
    QPainterPath path;

    for (int i = 0; i < controlPoints.size(); i += 4) {
        QPoint p0 = controlPoints[i];
        QPoint p1 = controlPoints[i + 1];
        QPoint c0 = controlPoints[i + 2];
        QPoint c1 = controlPoints[i + 3];

        path.moveTo(p0);
        path.cubicTo(c0.x(), c0.y(), c1.x(), c1.y(), p1.x(), p1.y());

    }

    // Top curve

    QPen blackPen(Qt::black, 6);
    blackPen.setCapStyle(Qt::RoundCap);

    painter.setPen(blackPen);
    painter.drawPath(path);

    QPen pen;
    pen.setWidth(3);
    pen.setColor(QColor(255, 25, 25));
    pen.setCapStyle(Qt::RoundCap);
    painter.setPen(pen);
    painter.drawPath(path);
}
void Visualizer::paintSheet() 
{
    glClearColor(0.85f, 0.85f, 0.85f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    QPainter painter(this);
    painter.setRenderHint(QPainter::Antialiasing);
    // +1 here so that loops are not on the edges
    int sheetWidth = sheet[0].size() + 1;
    int sheetHeight = sheet.size() + 1;

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

    QPoint p0;
    QPoint p1;
    QPoint c0;
    QPoint c1;
    std::vector<QPoint> controlPoints;

    for (int i = 0; i < sheet.size(); i++) {
        if (i == sheet.size() - 1) {
            p0 = QPoint(0, (i + 1) * cellHeight + yOffset);
            p1 = QPoint(xOffset, (i + 1) * cellHeight + yOffset);
            c0 = QPoint(0, (i + 1) * cellHeight + yOffset);
            c1 = QPoint(xOffset, (i + 1) * cellHeight + yOffset);
            controlPoints.push_back(p0);
            controlPoints.push_back(p1);
            controlPoints.push_back(c0);
            controlPoints.push_back(c1);
        }

        for (int j = 0; j < sheet[i].size(); j++) {
            // Calculate the position of the rectangle
            int rectX = j * cellWidth + xOffset; // Left of the cell
            int rectY = (sheet.size() - 1 - i) * cellHeight + yOffset; // Top of the cell

            //paint in yarn
            p0 = QPoint(rectX, rectY + cellHeight);
            p1 = QPoint(rectX + cellWidth / 2.4, rectY + (cellHeight / 1.25));
            c0 = QPoint(rectX + cellWidth / 4, rectY + (cellHeight));
            c1 = QPoint(rectX + cellWidth / 2, rectY + (cellHeight));
            controlPoints.push_back(p0);
            controlPoints.push_back(p1);
            controlPoints.push_back(c0);
            controlPoints.push_back(c1);

            int x = cellWidth * sheet[i][j].offset;
            p0 = p1;
            p1 = QPoint(rectX + cellWidth / 4 + x, rectY - cellHeight / 3);
            c0 = QPoint(rectX + (cellWidth / 2.6), rectY + (cellHeight / 1.15));
            c1 = QPoint(rectX + cellWidth / 5 + x, rectY - cellHeight / 4);
            controlPoints.push_back(p0);
            controlPoints.push_back(p1);
            controlPoints.push_back(c0);
            controlPoints.push_back(c1);

            p0 = QPoint(rectX + cellWidth / 4 + x, rectY - cellHeight / 3);
            p1 = QPoint(rectX + cellWidth - (cellWidth / 4) + x, rectY - cellHeight / 3);
            c0 = QPoint(rectX + cellWidth / 2.7 + x, rectY - (cellHeight / 2));
            c1 = QPoint(rectX + cellWidth - (cellWidth / 2.7) + x, rectY - (cellHeight / 2));
            controlPoints.push_back(p0);
            controlPoints.push_back(p1);
            controlPoints.push_back(c0);
            controlPoints.push_back(c1);


            p0 = QPoint(rectX + cellWidth - (cellWidth / 4) + x, rectY - cellHeight / 3);
            p1 = QPoint(rectX + cellWidth - (cellWidth / 2.4), rectY + (cellHeight / 1.25));
            c0 = QPoint(rectX + cellWidth - (cellWidth / 5) + x, rectY - cellHeight / 4);
            c1 = QPoint(rectX + cellWidth - (cellWidth / 2.6), rectY + (cellHeight / 1.15));
            controlPoints.push_back(p0);
            controlPoints.push_back(p1);
            controlPoints.push_back(c0);
            controlPoints.push_back(c1);

            p0 = QPoint(rectX + cellWidth, rectY + cellHeight);
            p1 = QPoint(rectX + cellWidth - (cellWidth / 2.4), rectY + (cellHeight / 1.25));
            c0 = QPoint(rectX + cellWidth - (cellWidth / 4), rectY + (cellHeight));
            c1 = QPoint(rectX + cellWidth - (cellWidth / 2), rectY + (cellHeight));
            controlPoints.push_back(p0);
            controlPoints.push_back(p1);
            controlPoints.push_back(c0);
            controlPoints.push_back(c1);
        }
        if (i != 0) {
            if ((i % 2 == 1 && sheet.size() % 2 == 0) || (i % 2 == 0 && sheet.size() % 2 == 1)) {
                p0 = QPoint(sheet[i].size() * cellWidth + xOffset, i * cellHeight + yOffset);
                p1 = QPoint(sheet[i].size() * cellWidth + xOffset, (i + 1) * cellHeight + yOffset);
                c0 = QPoint(sheet[i].size() * cellWidth + xOffset + (cellWidth / 5), i * cellHeight + yOffset + (cellHeight - cellHeight / 1.7));
                c1 = QPoint(sheet[i].size() * cellWidth + xOffset + (cellWidth / 5), i * cellHeight + yOffset + cellHeight / 1.7);
                controlPoints.push_back(p0);
                controlPoints.push_back(p1);
                controlPoints.push_back(c0);
                controlPoints.push_back(c1);
            }
            else {
                p0 = QPoint(xOffset, i * cellHeight + yOffset);
                p1 = QPoint(xOffset, (i + 1) * cellHeight + yOffset);
                c0 = QPoint(xOffset - (cellWidth / 5), i * cellHeight + yOffset + (cellHeight - cellHeight / 1.7));
                c1 = QPoint(xOffset - (cellWidth / 5), i * cellHeight + yOffset + cellHeight / 1.7);
                controlPoints.push_back(p0);
                controlPoints.push_back(p1);
                controlPoints.push_back(c0);
                controlPoints.push_back(c1);
            }
        }
        else {
            if (sheet.size() % 2 == 1) {
                p0 = QPoint(sheet[i].size() * cellWidth + xOffset, (i + 1) * cellHeight + yOffset);
                p1 = QPoint(sheet[i].size() * cellWidth + xOffset + xOffset, (i + 1) * cellHeight + yOffset);
                c0 = QPoint(sheet[i].size() * cellWidth + xOffset, (i + 1) * cellHeight + yOffset);
                c1 = QPoint(sheet[i].size() * cellWidth + xOffset, (i + 1) * cellHeight + yOffset);
                controlPoints.push_back(p0);
                controlPoints.push_back(p1);
                controlPoints.push_back(c0);
                controlPoints.push_back(c1);
            }
            else {
                p0 = QPoint(0, (i + 1) * cellHeight + yOffset);
                p1 = QPoint(xOffset, (i + 1) * cellHeight + yOffset);
                c0 = QPoint(0, (i + 1) * cellHeight + yOffset);
                c1 = QPoint(xOffset, (i + 1) * cellHeight + yOffset);
                controlPoints.push_back(p0);
                controlPoints.push_back(p1);
                controlPoints.push_back(c0);
                controlPoints.push_back(c1);
            }
        }
    }
    paintCurve(painter, controlPoints);
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
        glClearColor(clearColor.r, clearColor.g, clearColor.b, clearColor.a);
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
            paintInterpolatedMesh();
        }
        if (showFirstChains) {
            if (showSlice) {
                paintSliceMesh();
            }
            paintRowColGraph();
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
        if (showYarn) {
            paintYarnMesh();
        }
    }
    else {
        paintSheet();
    }
}

//Slot functions that set the relevant visualization parameters and flags
void Visualizer::knitGraphCreated() 
{
    showLastChain = true;
}
void Visualizer::knitGraphTraced(std::vector< TracedStitch >* t) 
{
    tracedMesh = *t;
    showTraced = true;
    showFirstChains = false;
    showInterpolated = false;
}
void Visualizer::peelSliceDone(ObjectMesh slice, std::vector< std::vector< uint32_t > > slice_active_chains_, std::vector< std::vector< uint32_t > > slice_next_chains_) 
{
    showInterpolated = false;
    sliceMesh = slice;
    sliceActiveChains = slice_active_chains_;
    sliceNextChains = slice_next_chains_;

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sliceMesh.vertices.size() * sizeof(QVector3D), sliceMesh.vertices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sliceMesh.indices.size() * sizeof(GLuint), sliceMesh.indices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), nullptr);
    glEnableVertexAttribArray(0);
    glFinish();

    showSlice = true;
    showSliceChains = true;
}
void Visualizer::nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >* next_chains) 
{
    nextActiveChains = *next_chains;

    showLinks = false;
    showSliceChains = false;
    showNextChains = true;
}
void Visualizer::linkChainsDone(std::vector< std::vector< Stitch > >* next_stitches, std::vector< Link >* links_) 
{
    nextStitches = *next_stitches;
    links = *links_;

    showLinks = true;
}

//helper function that returns the position of the loop in the sheet
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
                //qDebug() << "Mouse is in cell: " << i << ", " << j;
                res.first = i;
                res.second = j;
                //chosenVertex = i * sheet[i].size() + j;
                //qDebug() << "Chosen vertex: " << chosenVertex;
            }
        }
    }
    return res;
}

//no longer brazenly copied, helper function that linearly interpolated between 2 vertices
QVector3D Visualizer::mmix(const QVector3D& wa, const QVector3D& wb, float m) 
{
    return wa + m * (wb - wa);
}
//helper function that copied the locations of the vertices of the mesh
std::vector< std::vector< QVector3D > > Visualizer::copyLocations(ObjectMesh const& model, std::vector< std::vector< uint32_t > > const& chains) {
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
//helper function that interpolated between stitch locations
std::vector< std::vector< QVector3D > > Visualizer::interpolateStitchLocations(std::vector< std::vector< QVector3D > > const& chains, std::vector< std::vector< Stitch > > const& stitches) 
{
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
//same as previous function, but done before the stitch locations are calculated
std::vector< std::vector< QVector3D > > Visualizer::interpolateLocations(ObjectMesh const& model, std::vector< std::vector< EmbeddedVertex > > const& chains) {
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

/*
*   Function:compute the object boundaries so that the camera and mesh rotation can be set up correctly
*/
void Visualizer::computeBoundaries() {
    float sumX, sumY, sumZ;
    float maxDistSquared;

    sumX = sumY = sumZ = 0.0f;
    maxDistSquared = 0.0f;

    for (int i = 0; i < mesh.vertices.size(); i++) {
        sumX += mesh.vertices[i].x();
        sumY += mesh.vertices[i].y();
        sumZ += mesh.vertices[i].z();
    }

    int num_vertices = mesh.vertices.size();
    float centroidX = sumX / num_vertices;
    float centroidY = sumY / num_vertices;
    float centroidZ = sumZ / num_vertices;

    centroid = QVector3D(centroidX, centroidY, centroidZ);
}
/*
*   Function: load the ObjectMesh object into the Visualizer module
*/
void Visualizer::loadMesh(ObjectMesh object) {
    mesh = object;

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, mesh.vertices.size() * sizeof(QVector3D), mesh.vertices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexBuffer);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, mesh.indices.size() * sizeof(GLuint), mesh.indices.data(), GL_STATIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(QVector3D), nullptr);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);


    computeBoundaries();
    
    showModel = true;
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
        //int cnt = std::count(currentConstraint->vertices.begin(), currentConstraint->vertices.end(), chosenVertex);

        //if (cnt > 0) {
            //qDebug() << "Vertex already in constraint...";
            //return;                                                          //WHY WAS THIS NOT HERE BUT IT STILL WORKED?
        //}
        currentConstraint->vertices.push_back(chosenVertex);
    }
    else {
        qDebug() << "No vertex in selection...";
    }
}
//public function that gives the constraints created in visualizer to the caller (Knittee)
std::vector<Constraint*> Visualizer::getConstraints()
{
    return constraints;
}
//public function that sets the constraints from Knittee to Visualizer
void Visualizer::setConstraints(std::vector<Constraint*> c)
{
    constraints = c;
    currentConstraint = new Constraint();
    constraints.push_back(currentConstraint);
}
//push the current constraints one back
void Visualizer::pushConstraints() {
    currentConstraint = new Constraint();
    constraints.push_back(currentConstraint);
}
//get closest constraint to vertex
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
//increase the time value of the closest constraint
void Visualizer::increaseConstraintValue() {
    int closestConstraint = getClosestConstraint();
    if (closestConstraint != -1 && constraints[closestConstraint]->timeValue < 1.0) {
		constraints[closestConstraint]->timeValue += 0.1f;
	}
}
//decrease the time value of the closest constraint
void Visualizer::decreaseConstraintValue() {
	int closestConstraint = getClosestConstraint();
    if (closestConstraint != -1 && constraints[closestConstraint]->timeValue > -1.0) {
		constraints[closestConstraint]->timeValue -= 0.1f;
	}
}
//delete the closest constraint
void Visualizer::deleteConstraint() {
    int closestConstraint = getClosestConstraint();
    if (closestConstraint != -1) {
		constraints.erase(constraints.begin() + closestConstraint);
    }
}


/*
* Handlers of different events follow:
*/
void Visualizer::keyPressEvent(QKeyEvent* event)
{
    //if 'c', user wants to add a vertex
    if (event->key() == Qt::Key_C && addingConstraints) {
        addConstraintVertex();
    }
    //if 'enter', user wants to finish adding constraint vertices
    if (event->key() == Qt::Key_Return && addingConstraints) {
        pushConstraints();
    }
    //if '+', user wants to increase time value
    if (event->key() == Qt::Key_Plus ) {
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
        // in the future will backtrack many options, for now only constraints
        if (currentConstraint->vertices.size() > 0) {
            currentConstraint->vertices.pop_back();
        }
    }
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

    if (isMovingLoop) {
        mouseX = event->pos().x();
        mouseY = (height() - event->pos().y() - 1);
        
        to = getLoopPos(mouseX, mouseY);

        if (from.first != -1 && to.first != -1) {
            emit moveLoop(from, to);
        }
    }

    if (isRotating) {
        // calculate the change in mouse position between previous and current positions
        float deltaX = currentPos.x() - lastMousePos.x();
        float deltaY = currentPos.y() - lastMousePos.y();

        // update the rotation angles based on the mouse movement
        // the 0.5 value could be a variable set by the user as a rotation speed 
        m_rotationAngleY += deltaX * 0.5f;
        m_rotationAngleX += deltaY * 0.5f;

        // make sure the rotation angles are within the 360 degrees range
        m_rotationAngleX = loopAround(m_rotationAngleX, 0.0f, 359.999f);
        m_rotationAngleY = loopAround(m_rotationAngleY, 0.0f, 359.999f);

        lastMousePos = currentPos;
    }
    if (isDragging) {
        QPoint currentPos = event->pos();
        QVector3D viewDirection;

        // get the view direction so that the object can be dragged in all axes not just in the X and Y direction
        viewDirection.setX(sin(qDegreesToRadians(m_rotationAngleY)) * cos(qDegreesToRadians(m_rotationAngleX)));
        viewDirection.setY(-sin(qDegreesToRadians(m_rotationAngleX)));
        viewDirection.setZ(-cos(qDegreesToRadians(m_rotationAngleY)) * cos(qDegreesToRadians(m_rotationAngleX)));

        // calculate the change in mouse position
        float deltaX = currentPos.x() - lastMousePos.x();
        float deltaY = currentPos.y() - lastMousePos.y();

        // Update the drag values based on the mouse movement
        // currently doing divisions by zoomLevel to get a good dragging speed, could use a variable and have the user set it
        translateX += deltaX * 0.05f / zoomLevel; 
        translateY -= deltaY * 0.05f / zoomLevel;

        lastMousePos = currentPos;
    }
}
void Visualizer::wheelEvent(QWheelEvent* event)
{
    // zoom in/out based on the scroll wheel delta
    QPoint angleDelta = event->angleDelta();
    if (!angleDelta.isNull()) {
        float zoomFactor = 0.03f;
        if (angleDelta.y() > 0) {
            zoomLevel -= zoomFactor;
        }
        else {
            zoomLevel += zoomFactor;
        }
        isZooming = true;
    }
}

//timer functions for 60fps follow:
void Visualizer::stopTimer() 
{
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

//slot functions that allow the user to control which mesh is shown
void Visualizer::showInterpolatedChanged(int state) 
{
    if (state == Qt::Checked) {
        showInterpolated = true;
    }
    else if (state == Qt::Unchecked) {
        showInterpolated = false;
    }
}
void Visualizer::showGraphChanged(int state) 
{
    if (state == Qt::Checked) {
		showGraph = true;
        showFirstChains = true;
	}
    else if (state == Qt::Unchecked) {
		showGraph = false;
        showFirstChains = false;
	}
}

void Visualizer::showTracedChanged(int state) 
{
    if (state == Qt::Checked) {
		showTraced = true;
	}
    else if (state == Qt::Unchecked) {
		showTraced = false;
	}
}
void Visualizer::showYarnChanged(int state) 
{
    if (state == Qt::Checked) {
        showYarn = true;
    }
    else if (state == Qt::Unchecked) {
        showYarn = false;
    }
}
