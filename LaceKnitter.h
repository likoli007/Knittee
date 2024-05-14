#pragma once
#include <vector>
#include "FlatPoint.h"


#include <QObject>
#include <QMap>
#include <QList>
//this class is in charge of creating the sheet of stitches for a 2D project, manipulating it and storing it


struct LaceTransfer {
	QString from_bed;
	int from;
	QString to_bed;
	int to;
};

struct State {
	std::vector<int> current;
	std::vector<int> prev;
	std::vector<int> offsets;
	int Do;
	std::vector<LaceTransfer> path;
	int l;
	int r;
	std::vector<State*> chain;
	int rack;
	int penalty = -1;
};
//overload == operator that compares all fields


inline bool operator==(const State& lhs, const State& rhs) {
if (lhs.current != rhs.current) return false;
	if (lhs.prev != rhs.prev) return false;
	if (lhs.offsets != rhs.offsets) return false;
	if (lhs.l != rhs.l) return false;
	if (lhs.r != rhs.r) return false;
	if (lhs.rack != rhs.rack) return false;
	return true;
}




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

	void xfer( QString& fromBed, int fromIndex, QString& toBed, int toIndex);
	int xferredStacks = 0;
	int xferredEmpty = 0;
	bool ignoreFirsts = false;
	bool ignoreStacks = false;
	bool ignoreEmpty = false;
	bool strict_order = false;
	int chosenAlgorithm = 0;  //0 = sb, 1=cse, 2=?

	std::vector<State> priority_order(std::vector<State> states);
	bool has_state(State& state, std::vector<State*>& visited);
	void visit_state(State& state, std::vector<State*>& visited);
	int penalty(const State& state);

	int top = 0;

	constexpr int parent(int i) {
		return ((i + 1) >> 1) - 1;
	}

	constexpr int left(int i) {
		return (i << 1) + 1;
	}

	constexpr int right(int i) {
		return (i + 1) << 1;
	}



	void schoolbus(const std::vector<int>& offsets, int minRacking, int maxRacking);
	void cse(int maxRacking);
	void generate_transfers(std::vector<LaceTransfer>& path, bool log = false);
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

