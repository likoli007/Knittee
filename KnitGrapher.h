#pragma once

#include "ObjectMesh.h"
#include "Visualizer.h"
//#include "myQVectors.h"

#include <glm/glm.hpp>
#include <QVector3D>
#include "Stitch.h"
#include "RowColGraph.h"
#include "EmbeddedVertex.h"
//#include "Constraint.h"
#include "EmbeddedPlanarMap.h"
#include "RowColGraph.h"



/*
* This class will be used to construct the KnitGraph from the parameters given by the user
* A lot of this code will be implemented by refering to the algorithms described by AUTOKNIT
*/
class KnitGrapher : public QObject
{
	Q_OBJECT


private:
	ObjectMesh originalMesh, newMesh;
	float stitchWidth;
	float stitchHeight;		//both float values in mm
	float modelUnitLength;	//the length of a unit in the model in mm

	const float minEdgeRatio = 0.3f; //'smallest allowed smallest-to-largest edge ratio in a triangle' 
	const float minEdgeRatioSquared = minEdgeRatio * minEdgeRatio;
	//maybe let them be set by the user?

	std::vector<Constraint*> constraints;
	std::vector<float> constrained_values;

	std::vector<glm::uvec3> oldTriangles;
	std::vector<glm::uvec3> newTriangles;


	//////////////////////class member variables used in peeling///////////////////////////////
	std::vector< std::vector< EmbeddedVertex > >activeChains;
	std::vector< std::vector< Stitch > >activeStitches;
	RowColGraph graph;
	std::vector<EmbeddedVertex> sampledChain;

	ObjectMesh slice;
	std::vector< EmbeddedVertex > sliceOnModel;
	std::vector< std::vector< uint32_t > > sliceActiveChains;
	std::vector< std::vector< uint32_t > > sliceNextChains;
	std::vector< bool > usedBoundary;



	float getChainSampleSpacing() const {
		return 0.25f * stitchWidth / modelUnitLength;
	}

	void sampleChain(
		float spacing,
		std::vector< EmbeddedVertex > const& chain //in: chain to be sampled
		);


	void peelSlice();

	int stepCount = 0;

	void quad(std::vector<glm::uvec3>& new_tris, std::vector<QVector3D>& verts, GLuint a, GLuint b, GLuint c, GLuint d);
	void divide(QSet<QPoint>& marked,
		std::vector<QVector3D>& verts,
		std::vector<std::vector<GLuint>>& paths,
		std::vector<glm::uvec3>& tris);
	void remesh();
	float getMaxEdgeLength();
	bool degenerateCheck(std::vector<glm::uvec3> tris);

	std::vector<GLuint> toIntArray(std::vector<glm::uvec3>);

	void generateTriangles();
	void interpolateValues();
	void printConstrainedValues();

	void findFirstActiveChains();

	//so. much. passing. by. reference. need to make some variables class members.... TODO!!!
	void unfold(GLuint depth, GLuint root, QVector2D const& flat_root,
		GLuint ai, QVector2D const& flat_a,
		GLuint bi, QVector2D const& flat_b,
		QVector2D const& limit_a, QVector2D const& limit_b, QHash< QPoint, GLuint > const& opposite,
		std::vector<QVector3D> const& newVertices, QHash< QPoint, float >& min_dis);

	GLuint lookup(GLuint a, GLuint b, QHash<QPoint, GLuint>& marked_verts);
public:
	KnitGrapher(QObject* parent = nullptr);
	void constructKnitGraph(std::vector<Constraint*> constraints);
public slots:
	void setStitchWidth(float width);
	void setStitchHeight(float height);
	void setModelUnitLength(float length);
	void setOriginalMesh(ObjectMesh mesh);
	void stepButtonClicked();
signals:
	void knitGraphInterpolated(ObjectMesh mesh, std::vector<float> values);

};

