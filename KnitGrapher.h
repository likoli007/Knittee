#pragma once

#include "ObjectMesh.h"
#include "Visualizer.h"
/*
* This class will be used to construct the KnitGraph from the parameters given by the user
* A lot of this code will be implemented by refering to the algorithms described by AUTOKNIT
*/
class KnitGrapher : public QObject
{
	Q_OBJECT


private:
	ObjectMesh originalMesh;
	float stitchWidth;
	float stitchHeight;		//both float values in mm
	float modelUnitLength;	//the length of a unit in the model in mm

	const float minEdgeRatio = 0.3f; //'smallest allowed smallest-to-largest edge ratio in a triangle' 
	const float minEdgeRatioSquared = minEdgeRatio * minEdgeRatio;
	//maybe let them be set by the user?
	/*
	void quad(std::vector<GLuint>& new_tris, std::vector<QVector3D>& verts, uint32_t a, uint32_t b, uint32_t c, uint32_t d);
	void divide(QSet<QPoint>& marked, std::vector<QVector3D>& verts, std::vector<std::vector<int>>& paths, std::vector<GLuint>& tris);
	ObjectMesh remesh(std::vector<Constraint*> constraints);
	float getMaxEdgeLength();
	bool degenerateCheck();

	//so. much. passing. by. reference. need to make some variables class members.... TODO!!!
	void unfold(GLuint depth, GLuint root, QVector2D const& flat_root,
		GLuint ai, QVector2D const& flat_a,
		GLuint bi, QVector2D const& flat_b,
		QVector2D const& limit_a, QVector2D const& limit_b, QHash< QPoint, GLuint > const& opposite,
		std::vector<QVector3D> const& newVertices, QHash< QPoint, float >& min_dis);

	int lookup(int a, int b, QHash<QPoint, int> &marked_verts);*/
public:
	KnitGrapher(QObject* parent = nullptr);
	void constructKnitGraph(std::vector<Constraint*> constraints);
public slots:
	void setStitchWidth(float width);
	void setStitchHeight(float height);
	void setModelUnitLength(float length);
	void setOriginalMesh(ObjectMesh mesh);

};

