#include "FlatPoint.h"



FlatPoint::FlatPoint()
{
	stitch = "k";
	offset = 0;
	first = false;
}


QDataStream& operator<<(QDataStream& out, const FlatPoint& obj)
{
	out << obj.stitch << obj.offset << obj.first;
	return out;
}

QDataStream& operator>>(QDataStream& in, FlatPoint& obj)
{
	in >> obj.stitch >> obj.offset >> obj.first;
	return in;
}