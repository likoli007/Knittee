#pragma once

#include "PlanTransfers.h"

#include <cassert>
#include <iostream>
#include <algorithm>

struct NeedleRollGoal {
	NeedleRollGoal() = default;
	NeedleRollGoal(int32_t needle_, int32_t roll_, int32_t goal_, int32_t left_slack_, int32_t right_slack_) : needle(needle_), roll(roll_), goal(goal_), left_slack(left_slack_), right_slack(right_slack_) {
	}
	int32_t needle = 0;
	int32_t roll = 0;
	int32_t goal = 0;

	//slack available on both sides of the stitch:
	Slack left_slack = 0;
	Slack right_slack = 0;

	//can stitch be stacked to the left/right?
	bool can_stack_left = false;
	bool can_stack_right = false;

	NeedleRollGoal after_offset_and_roll(int32_t offset, int32_t roll) const {
		NeedleRollGoal ret = *this;
		ret.needle += offset;
		//adjust target roll based on transfer's roll:
		ret.roll = ret.roll - roll;
		//adjust roll/slack if on the other bed:
		if (roll % 2) {
			ret.roll = -ret.roll;
			std::swap(ret.left_slack, ret.right_slack);
			std::swap(ret.can_stack_left, ret.can_stack_right);
		}
		return ret;
	}

	bool has_same_goal_as(NeedleRollGoal const& o) const {
		return needle == o.needle && roll == o.roll && goal == o.goal;
	}

	//compare 'true' even if the other goal has 2*N more rolling involved.
	bool has_same_real_goal_as(NeedleRollGoal const& o) const {
		return needle == o.needle && (roll - o.roll) % 2 == 0 && goal == o.goal;
	}

	uint32_t penalty(int32_t min_free, int32_t max_free) const {
		uint32_t dir = penalty_dir(min_free, max_free);
		uint32_t rec = penalty_rec(min_free, max_free);
		//std::cout << "[n,r,g] = [" << needle << "," << roll << "," << goal << "] -> dir: " << dir << " vs rec: " << rec << std::endl;
		assert(dir == rec);
		return dir;
	}
	//These should be the same, but for now both are called and checked to shake out any bugs:
	uint32_t penalty_dir(int32_t min_free, int32_t max_free) const {
		if (roll == 0) {
			return std::abs(goal - needle);
		}
		else if (roll < 0) {
			return (needle - min_free)
				+ 2 * -roll
				+ (-roll - 1) * (max_free - min_free)
				+ (roll % 2 == 0 ? max_free - goal : goal - min_free);
		}
		else { //(roll > 0)
			return (max_free - needle)
				+ 2 * roll
				+ (roll - 1) * (max_free - min_free)
				+ (roll % 2 == 0 ? goal - min_free : max_free - goal);
		}
	}
	uint32_t penalty_rec(int32_t min_free, int32_t max_free) const {
		if (roll == 0) {
			return std::abs(goal - needle);
		}
		else if (roll < 0) {
			NeedleRollGoal after = *this;
			after.needle = min_free;
			after.roll = -(roll + 1);
			return (needle - min_free) + 2 + after.penalty_rec(min_free, max_free);
		}
		else { //(roll > 0)
			NeedleRollGoal after = *this;
			after.needle = max_free;
			after.roll = -(roll - 1);
			return (max_free - needle) + 2 + after.penalty_rec(min_free, max_free);
		}
	}
};

/* This might be useful but might also be overkill:

//the code generally takes two views of a transfer problem, lists of from/to stitches and flattened "top/bottom" arrays.
//notably, it's much easier to simulate transfers on from/to lists, but easier to actually run dyanmic programming on flattened top/bottom arrays.

struct FromTo;

//TopBottom stores a transfer problem in terms of "NeedleRollGoal" instances, which express the current needle of each stitch as well as its target bed (in terms of roll) and needle (goal).
struct TopBottom {
	std::vector< NeedleRollGoal > top; //top stitches, sorted left-to-right
	std::vector< NeedleRollGoal > bottom; //bottom stitches, sorted left-to-right
	BedNeedle::Bed top_bed;
	BedNeedle::Bed bottom_bed;

	TopBottom() = default;
	TopBottom(BedNeedle::Bed top_bed, BedNeedle::Bed bottom_bed, FromTo const &source); //construct from a FromTo instance
};

//FromTo stores a cycle in terms of source and destination needles in ccw order:
struct FromTo {
	std::vector< BedNeedle > from_ccw;
	std::vector< BedNeedle > to_ccw;
	std::vector< Slack > slack;

	FromTo() = default;
	FromTo(TopBottom const &source); //construct from a TopBottom instance
};
*/


//helper for constructing roll values:
void minimize_winding(
	std::vector< int32_t >* winding
);

//helper for drawing an (ASCII / ANSI) picture of beds:
void draw_beds(
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom
);


//helper for constructing output beds:
//NOTE: *recomputes* roll values [only cares about roll parity], cleans up doubled stitches.
void run_transfers(
	Constraints const& constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom,
	std::vector< Transfer > const& plan,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal >* to_top,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal >* to_bottom
);

//NOTE: all the functions below expect *no overlapped stitches* in the top/bottom input arrays, but may overlap stitches in the output arrays:

//collapse 'top' onto 'bottom':
void best_collapse(
	Constraints const& constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal >* to_top,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal >* to_bottom,
	std::vector< Transfer >* plan
);

//shift top/bottom and move to other bed:
void best_shift(
	Constraints const& constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal >* to_top,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal >* to_bottom,
	std::vector< Transfer >* plan
);

//expand bottom to other bed:
void best_expand(
	Constraints const& constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal >* to_top,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal >* to_bottom,
	std::vector< Transfer >* plan
);

