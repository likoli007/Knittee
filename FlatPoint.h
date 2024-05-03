#pragma once

#include <QString>
#include <QDataStream>
class FlatPoint
{
public:
    FlatPoint();

	QString stitch;
	int offset;
	bool first;

    // Declare operator<< as a friend function
    friend QDataStream& operator<<(QDataStream& out, const FlatPoint& obj);
    // Declare operator>> as a friend function
    friend QDataStream& operator>>(QDataStream& in, FlatPoint& obj);
};

