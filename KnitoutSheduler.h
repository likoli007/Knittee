#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "KnitOutStitch.h"
#include "Shape.h"
#include "Typeset.h"
#include "ScheduleCost.h"
#include "EmbedDAG.h"
#include "PlanTransfers.h"



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


public:
	KnitoutScheduler(QObject* parent = nullptr);
};

