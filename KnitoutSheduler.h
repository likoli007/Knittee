#pragma once
#include <cstdint>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <functional>
#include <string>
#include <qdebug.h>
#include <QFile>
#include <deque>
#include <QDebug>
#include <qregularexpression.h>


#include "KnitOutStitch.h"
#include "Shape.h"
#include "Typeset.h"
#include "ScheduleCost.h"
#include "EmbedDAG.h"
#include "PlanTransfers.h"
#include "ScheduleCost.h"
#include "KnitOutStitch.h"


typedef ScheduleCost DAGCost;

typedef uint32_t DAGNodeIndex;
typedef uint32_t DAGEdgeIndex;

struct DAGEdge {
	DAGNodeIndex from = -1U;
	DAGNodeIndex to = -1U;
	uint32_t from_shapes = 0;
	uint32_t to_shapes = 0;
	std::vector< DAGCost > costs;
	//If from != -1U && to != -1U: cost[f * to_shapes + t] --> cost given from/to shape
	//If from == -1U: cost[t] is min-cost given to shape
	//If to == -1U: cost[f] is min-cost given from shape
};

struct DAGOption {
	DAGCost cost;
	std::vector< DAGEdgeIndex > in_order;
	std::vector< DAGEdgeIndex > out_order;
	std::vector< uint32_t > in_shapes;
	std::vector< uint32_t > out_shapes;
};

struct DAGNode {
	std::vector< DAGOption > options;
};






//Ported from javascript
struct PortedBedNeedle {
	QString bed;
	int needle;
};

//Loop held on a needle:
struct Loop {
	constexpr Loop(uint32_t stitch_, uint32_t idx_) : stitch(stitch_), idx(idx_) {
	}
	uint32_t stitch;
	uint32_t idx;
	bool operator==(Loop const& o) const {
		return (stitch == o.stitch && idx == o.idx);
	}
	bool operator!=(Loop const& o) const {
		return (stitch != o.stitch || idx != o.idx);
	}
	bool operator<(Loop const& o) const {
		if (stitch != o.stitch) return stitch < o.stitch;
		else return idx < o.idx;
	}
	std::string to_string() const {
		if (stitch == -1U && idx == -1U) return "GAP";
		else return std::to_string(stitch) + "_" + std::to_string(idx);
	}
};
constexpr const Loop INVALID_LOOP = Loop(-1U, -1U);


struct CycleIndex {
	CycleIndex(uint32_t cycle_, uint32_t index_) : cycle(cycle_), index(index_) { }
	uint32_t cycle = -1U;
	uint32_t index = -1U;
	bool operator!=(CycleIndex const& o) const {
		return cycle != o.cycle || index != o.index;
	}
	std::string to_string() const {
		if (cycle == -1U && index == -1U) return ".";
		return std::to_string(int32_t(cycle)) + "-" + std::to_string(int32_t(index));
	}
	std::string to_string_simple() const {
		if (cycle == -1U && index == -1U) return ".";
		if (cycle < 26) {
			if (index == 0) {
				return std::string() + char('A' + cycle);
			}
			else if (index == 1) {
				return std::string() + '+';
			}
			else {
				return std::string() + char('a' + cycle);
			}
		}
		else {
			if (index == 0) {
				return std::string() + '*';
			}
			else {
				return std::string() + 'x';
			}
		}
	}
};

class KnitoutScheduler : public QObject
{
	Q_OBJECT

	uint32_t maxRacking = 3;					//TODO: make a parameter settable by the user

	std::vector< KnitOutStitch > stitches;
	std::vector< std::pair<BedNeedle, BedNeedle>> stitch_locations, input_locations;

	std::vector<QString> scheduledInstructions;

	std::vector<QString> knitoutInstructions;

public:
	KnitoutScheduler(QObject* parent = nullptr);
	void setFilePath(QString path);
	void schedule();
signals:
	void instructionsCreated(std::vector<std::string>);
	void knitoutGenerated(std::vector<QString>);
	void helpBoxCommunication(QString message);
private:
	//NOTE: edges must have from < to (i.e., nodes should be topologically sorted)
	bool embed_DAG(
		std::vector< DAGNode > const& nodes,
		std::vector< DAGEdge > const& edges,
		std::vector< uint32_t >* node_options,
		std::vector< int32_t >* node_positions, //positions give total left-to-right order of edges/nodes
		std::vector< int32_t >* edge_positions
	);

	QString filePath;

////////Helper functions AND VARS ported from javascript/////////////////
	QString Carrier = "5";
	int StartRows = 10;
	int EndRows = 10;

	int YarnInStitchNumber = 67;
	int CastOnStitchNumber = 67;
	int PlainStitchNumber = 68;

	std::vector<QString> knitout;
	std::map<QString, bool> inYarns;
	std::map<QString, int> yarnUsedCount;
	std::map<QString, bool> needsReleaseHook;
	int racking = 0;
	bool needIn = true;

	PortedBedNeedle parseBedNeedle(QString& bedNeedle);
	QString bnToHalf(QString& line);
	void startTube(QString& dir, QStringList& bns);
	void bringIn(QString& carrier);
	void knit(QString& d, QString& bn);
	void drop(QString& bn);
	void tuck(QString& d, QString& bn);
	void miss(QString& d, QString& bn);
	void decrease(QString& d, QString& bn);
	void increase(QString& d0, QString& bn0, QString& d1, QString& bn1);
	void xfer(QString& from, QString& to);
	void setRacking(QString& from_str, QString& to_str);
	void xferCycle(QString& opts, QString& from, QString& to, QList<QPair<QString, QString>>& xfers);
	void kStash(QStringList& from, QStringList& to);
	void kUnstash(QStringList& from, QStringList& to);
	void endTube(QString& dir, QStringList& bns);

	void out (QString& str) {
		knitout.push_back(str);
	}


};

