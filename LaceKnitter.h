#pragma once
#include <vector>
#include "FlatPoint.h"


#include <QObject>
//this class is in charge of creating the sheet of stitches for a 2D project, manipulating it and storing it

class LaceKnitter : public QObject
{
	Q_OBJECT

	std::vector<std::vector<FlatPoint>> sheet;
	int maxOffset = 3;

	int sheetWidth, sheetHeight;
public:
	LaceKnitter(QObject* parent = nullptr);


	//LaceKnitter();
	//~LaceKnitter();
	
	void createSheet(int w, int h, int racking);
	void setDimensions(int w, int h, int racking);


	void saveToFile(QString);
	void loadFromFile(QString);

	int getWidth();
	int getHeight();
	int getRacking();

signals:
	void sheetChanged(std::vector<std::vector<FlatPoint>>* sheet);
	void sheetDimensionsChanged(int, int, int);
public slots:
	void moveLoop(QPair<int, int>, QPair<int, int>);
	void rackingChanged(int);

	void widthChanged(int, int);
	void heightChanged(int, int);
};

