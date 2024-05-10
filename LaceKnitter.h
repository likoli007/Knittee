#pragma once
#include <vector>
#include "FlatPoint.h"

#include <QObject>
//this class is in charge of creating the sheet of stitches for a 2D project, manipulating it and storing it

class LaceKnitter : public QObject
{
	Q_OBJECT

	std::vector<std::vector<FlatPoint>> sheet;
	int maxOffset = 10;

	int sheetWidth, sheetHeight;
public:
	LaceKnitter(QObject* parent = nullptr);


	//LaceKnitter();
	//~LaceKnitter();
	
	void createSheet(int w, int h);
	void setDimensions(int w, int h);


	void saveToFile(QString);
	void loadFromFile(QString);


signals:
	void sheetChanged(std::vector<std::vector<FlatPoint>>* sheet);


};

