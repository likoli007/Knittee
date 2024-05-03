#include "LaceKnitter.h"

#include <QFile>
#include <QDataStream>
#include <QDebug>

LaceKnitter::LaceKnitter(QObject* parent) : QObject(parent) {
	qDebug() << "loaded LaceKnitter";
}

void LaceKnitter::createSheet(int w, int h)
{
	qDebug() << "creating sheet";
	sheet.clear();

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


	qDebug() << "sheet created";
}

void LaceKnitter::setDimensions(int w, int h)
{
	sheetWidth = w;
	sheetHeight = h;
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
	qDebug() << "saved to file";
}