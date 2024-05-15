#pragma once
#include <vector>
#include "FlatPoint.h"


#include <QObject>
#include <QMap>
#include <QList>
#include <QFile>
#include <QDataStream>
#include <QDebug>
#include <qregularexpression.h>
#include <algorithm>
#include <cmath>
#include <set>
#include <QHash>
#include <queue>
#include <climits>
#include <regex>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include <QJsonValue>
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

	////////////////KNITOUT///////////////////////////
	void generateKnitout(int algo);
	int racking = 0;
	QStringList k;
	std::vector<int> offsets;
	QMap<QString, QVector<int>> needles;
	std::vector<bool> firsts;

	void xfer(QString& fromBed, int fromIndex, QString& toBed, int toIndex);
	int xferredStacks = 0;
	int xferredEmpty = 0;
	bool ignoreFirsts = false;
	bool ignoreStacks = false;
	bool ignoreEmpty = false;
	bool strict_order = false;
	int chosenAlgorithm = 0;  //0 = sb, 1=cse, 2=?




	void schoolbus(const std::vector<int>& offsets, int minRacking, int maxRacking);
	void cse(int maxRacking);

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

