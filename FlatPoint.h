#pragma once

#include <QString>
#include <QDataStream>

/*
*   FlatPoint represents a single loop on a 2D mesh, a vector of flat points is a sheet
*/

class FlatPoint
{
public:
    FlatPoint();

	QString stitch;
	int offset;
	bool first;

	//is saved as a binary file to save space, maybe JSONs? but theyre too slow
    friend QDataStream& operator<<(QDataStream& out, const FlatPoint& obj);
    friend QDataStream& operator>>(QDataStream& in, FlatPoint& obj);
};

