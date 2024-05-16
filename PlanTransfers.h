#pragma once

#include <limits>
#include <vector>
#include <string>
#include <cassert>

/*
*	This file contains helper structs and functions used to calculate an optimal transfer plan for a set of scheduled stitches.
*		BedNeedle: represents a needle on a knitting machine. The needle is identified by its bed (front, front sliders, back sliders, back) and its number.
*		Constraints: represents the constraints for the transfer plan. The constraints are the minimum and maximum number of free needles, and the maximum racking.
* 		Transfer: represents a transfer operation. It contains the source and destination needles, and a string explaining the reason for the transfer.
*		Slack: represents the slack between a stitch and its counterclockwise neighbor.
*		plan_transfers(): calculates an optimal transfer plan  
*	BedNeedle was renamed since the previous class name was already defined in the project.
*/


struct BedNeedle {
	enum Bed : char {
		Front = 'f',
		FrontSliders = 'F',
		BackSliders = 'B',
		Back = 'b'
	} bed;
	int32_t needle;
	bool operator==(BedNeedle const& o) const {
		return bed == o.bed && needle == o.needle;
	}
	bool operator!=(BedNeedle const& o) const {
		return !(*this == o);
	}

	BedNeedle(Bed bed_ = Front, int32_t needle_ = 0) : bed(bed_), needle(needle_) { }
	std::string to_string() const {
		if (bed == Front)             return "f" + std::to_string(needle);
		else if (bed == FrontSliders) return "fs" + std::to_string(needle);
		else if (bed == BackSliders)  return "bs" + std::to_string(needle);
		else if (bed == Back)         return "b" + std::to_string(needle);
		else {
			assert(0 && "Should have handled all cases");
			return "?" + std::to_string(needle); //never executed
		}
	}
};

template< >
struct std::hash< BedNeedle > {
	size_t operator()(BedNeedle const& bn) const {
		return (size_t(bn.bed) << (8 * (sizeof(size_t) - sizeof(bn.bed)))) ^ size_t(bn.needle);
	};
};

struct Constraints {
	int32_t min_free = std::numeric_limits< int32_t >::min();
	int32_t max_free = std::numeric_limits< int32_t >::max();
	uint32_t max_racking = 4;
};

struct Transfer {
	Transfer() = default;
	Transfer(BedNeedle const& from_, BedNeedle const& to_, std::string const& why_ = "") : from(from_), to(to_), why(why_) { }
	BedNeedle from;
	BedNeedle to;
	std::string why;
	std::string to_string() const {
		std::string ret = from.to_string() + " -> " + to.to_string();
		if (why != "") ret += " (" + why + ")";
		return ret;
	}
};

typedef int32_t Slack;
constexpr const Slack SlackForNoYarn = std::numeric_limits< Slack >::max();

bool plan_transfers(
	Constraints const& constraints,
	std::vector< BedNeedle > const& from_ccw,
	std::vector< BedNeedle > const& to_ccw,
	std::vector< Slack > const& slack, //slack between stitch and its ccw-ward neighbor
	std::vector< Transfer>* transfers,
	std::string* error = nullptr
);


