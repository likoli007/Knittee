#pragma once
#include "KnitGrapher.h"
//#include "EmbeddedVertex.h"
//#include "EmbeddedPlanarMap.h"
#include <vector>
#include <QSet>
#include <QPair>

KnitGrapher::KnitGrapher(QObject* parent) : QObject(parent)
{
	stitchWidth = 5;
	stitchHeight = 5;
	modelUnitLength = 1;

}

void KnitGrapher::setOriginalMesh(ObjectMesh mesh)
{
	originalMesh = mesh;
}

QPair<int, int> qMinMax(int a, int b)
{
	if (a < b)
		return QPair<int, int>(a, b);
	else
		return QPair<int, int>(b, a);
}



void KnitGrapher::constructKnitGraph(std::vector<Constraint*> constraints)
{
	qDebug() << "Constructing Knit Graph, sizes:" << stitchWidth << stitchHeight << modelUnitLength;
	// Step 1: create a new mesh that conforms to the constraints and stitch size given by the user
	qDebug() << "Constraints:" << constraints.size();

	//ObjectMesh remeshedMesh = remesh(constraints);

}
void KnitGrapher::setStitchWidth(float width)
{
	qDebug() << "Setting stitch width to " << width;
	stitchWidth = width;
}
void KnitGrapher::setStitchHeight(float height)
{
	stitchHeight = height;
}
void KnitGrapher::setModelUnitLength(float length)
{
	modelUnitLength = length;
}