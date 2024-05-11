#include "LaceKnitter.h"

#include <QFile>
#include <QDataStream>
#include <QDebug>

LaceKnitter::LaceKnitter(QObject* parent) : QObject(parent) {
	qDebug() << "loaded LaceKnitter";
}

void LaceKnitter::createSheet(int w, int h, int racking)
{
	sheet.clear();
	
	sheetWidth = w;
	sheetHeight = h;
	maxOffset = racking;
	for (int i = 0; i < h; i++)
	{
		std::vector<FlatPoint> row;
		for (int j = 0; j < w; j++)
		{
			FlatPoint fp;
			row.push_back(fp);
		}
		sheet.push_back(row);
	}
}

void LaceKnitter::setDimensions(int w, int h, int racking)
{
	sheetWidth = w;
	sheetHeight = h;
	maxOffset = racking;
}

void LaceKnitter::moveLoop(QPair<int, int> from, QPair<int, int> to)
{

	if (abs(to.second-from.second) > maxOffset){
		qDebug() << "cannot rack to this value!";
		return;
	}

	sheet[from.first][from.second].offset = to.second - from.second;

	emit sheetChanged(&sheet);
}

void LaceKnitter::rackingChanged(int r) {
	maxOffset = r;
}

void LaceKnitter::widthChanged(int increase, int left) {


	if (sheetWidth + increase < 2) {
		qDebug() << "cannot have a width of 0!";
		return;
	}

	sheetWidth += increase;

	if (increase == 1) {
		if (left == 1) {
			for (int i = 0; i < sheet.size(); i++) {
				FlatPoint fp;
				sheet[i].insert(sheet[i].begin(), fp);
			}
		}
		else {
			for (int i = 0; i < sheet.size(); i++) {
				FlatPoint fp;
				sheet[i].push_back(fp);
			}
		}
	}
	else {
		if (left == 1) {
			for (int i = 0; i < sheet.size(); i++) {
				sheet[i].erase(sheet[i].begin());
			}
		}
		else {
			for (int i = 0; i < sheet.size(); i++) {
				sheet[i].pop_back();
			}
		}
	}
	emit(sheetChanged(&sheet));
}

void LaceKnitter::heightChanged(int increase, int top) {

	//qDebug() << "arrrivaaa" << increase;
	if (sheetHeight + increase < 2) {
		qDebug() << "cannot have a height of 0!";
		return;
	}

	sheetHeight += increase;

	if (increase == 1) {
		if (top == 1) {
			std::vector<FlatPoint> row;
			for (int i = 0; i < sheetWidth; i++) {
				FlatPoint fp;
				row.push_back(fp);
			}
			sheet.insert(sheet.begin(), row);
		}
		else {
			std::vector<FlatPoint> row;
			for (int i = 0; i < sheetWidth; i++) {
				FlatPoint fp;
				row.push_back(fp);
			}
			sheet.push_back(row);
		}
	}
	else {
		if (top == 1) {
			sheet.erase(sheet.begin());
		}
		else {
			sheet.pop_back();
		}
	}
	emit sheetChanged(&sheet);
}

void LaceKnitter::loadFromFile(QString filename)
{
	sheet.clear();

	qDebug() << "loading from file";

	filename += "/sheet.dat";
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly))
	{
		qDebug() << "could not open file";
		return;
	}

	QDataStream in(&file);

	while (!in.atEnd())
	{
		
		for (int i = 0; i < sheetHeight; i++)
		{
			std::vector<FlatPoint> row;
			for (int j = 0; j < sheetWidth; j++) {
				FlatPoint fp;
				in >> fp;
				row.push_back(fp);
			}
			sheet.push_back(row);
		}
		
	}

	file.close();
	qDebug() << "loaded from file";
	emit sheetChanged(&sheet);
}

int LaceKnitter::getWidth()
{
	return sheetWidth;
}
int LaceKnitter::getHeight()
{
	return sheetHeight;
}
int LaceKnitter::getRacking()
{
	return maxOffset;
}

void LaceKnitter::saveToFile(QString filename)
{
	qDebug() << "saving to file";

	filename += "/sheet.dat";
	QFile file(filename);

	if (!file.open(QIODevice::WriteOnly))
	{
		qDebug() << "could not open file";
		return;
	}

	QDataStream out(&file);

	for (int i = 0; i < sheet.size(); i++)
	{
		for (int j = 0; j < sheet[i].size(); j++)
		{
			out << sheet[i][j];
		}
	}

	file.close();
	//qDebug() << "saved to file";
}