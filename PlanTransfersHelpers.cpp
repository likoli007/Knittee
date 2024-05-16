#include "PlanTransfersHelpers.h"
#include <sstream>
#include <algorithm>
#include <iostream>

#include <cassert>
#include <map>
#include <unordered_map>


/*
*	Implementation for the helper functions defined in the header file, used only by PlanTransfers.h
*/

void best_collapse(
	Constraints const& constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal >* to_top_,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal >* to_bottom_,
	std::vector< Transfer >* plan_
) {
	//Collapse won't change the bottom bed's location, but will change the top's:
	assert(top_bed != to_top_bed);
	assert(bottom_bed == to_bottom_bed);

	//make sure all output arrays exist:
	assert(to_top_);
	auto& to_top = *to_top_;
	assert(to_bottom_);
	auto& to_bottom = *to_bottom_;
	assert(plan_);
	auto& plan = *plan_;

	//Clear output:
	to_top.clear();
	to_bottom.clear();
	//NOTE: don't clear plan, only append.

	//------------
	//collapse moves stitches from top to bottom, from the edges in.

	//if no stitches on top, nothing to move, so done:
	if (top.empty()) {
		to_top.clear();
		to_bottom = bottom;
		return;
	}

	//Otherwise, Dijkstra-style search for optimal sequence of moves:

	//Search state:
#pragma pack(push,1)
	struct State {
		int32_t l; //index of left stitch on top
		int32_t r; //index of right stitch on top
		int32_t l_prev_needle; //needle of stitch to the left of 'l'
		int32_t r_next_needle; //needle of stitch to the right of 'r'
		enum : int8_t {
			LRollInvalid = -10,
			LRoll2 = -2,
			LRoll1 = -1,
			Roll0 = 0,
			RRoll1 = 1,
			RRoll2 = 2,
			RRollInvalid = 10
		};
		int8_t l_prev_roll; //{-10, -2, -1, 0} adjacent stitch is on the bottom bed to the left?
		int8_t r_next_roll; //{10, 2, 1, 0} adjacent stitch is on the bottom bed to the right?
		bool operator==(State const& o) const {
			return l == o.l
				&& r == o.r
				&& l_prev_needle == o.l_prev_needle
				&& r_next_needle == o.r_next_needle
				&& l_prev_roll == o.l_prev_roll
				&& r_next_roll == o.r_next_roll;
		};
		std::string to_string() const {
			return std::to_string(l_prev_needle) + "r" + std::to_string(int32_t(l_prev_roll))
				+ " [" + std::to_string(l) + "," + std::to_string(r) + "] "
				+ std::to_string(r_next_needle) + "r" + std::to_string(int32_t(r_next_roll))
				;
		}
	};
#pragma pack(pop)
	static_assert(sizeof(State) == 4 * 4 + 2 * 1, "collapse's State is packed");

	struct HashState {
		size_t operator()(State const& state) const {
			static std::hash< std::string > hash;
			return hash(std::string(reinterpret_cast<char const*>(&state), reinterpret_cast<char const*>(&state) + sizeof(state)));
		}
	};

	//Let's do this in terms of the actions that can be applied:
	struct Action {
		enum Type : uint8_t {
			None,
			MoveLeft,
			MoveRight,
			RollLeft,
			RollRight,
			Roll2Left,
			Roll2Right,
		} type;
		int32_t needle;

		Action(Type type_, int32_t needle_) : type(type_), needle(needle_) { }

		std::string to_string() const {
			if (type == None) return "CNone";
			else if (type == MoveLeft) return "CMoveLeft to " + std::to_string(needle);
			else if (type == MoveRight) return "CMoveRight to " + std::to_string(needle);
			else if (type == RollLeft) return "CRollLeft to " + std::to_string(needle);
			else if (type == RollRight) return "CRollRight to " + std::to_string(needle);
			else if (type == Roll2Left) return "CRoll2Left to " + std::to_string(needle);
			else if (type == Roll2Right) return "CRoll2Right to " + std::to_string(needle);
			else {
				assert(0 && "invalid move type");
				return "!"; //never reached
			}
		}
	};

	struct Cost {
		uint32_t penalty = 0;
		bool operator<(Cost const& o) const {
			return penalty < o.penalty;
		}
		bool operator==(Cost const& o) const {
			return penalty == o.penalty;
		}
	};

	struct StateInfo {
		Cost cost;
		State const* source;
		Action action;

		StateInfo(Cost const& cost_, State const* source_, Action const& action_) : cost(cost_), source(source_), action(action_) { }
	};

	std::multimap< Cost, const State* > todo;
	std::unordered_map< State, StateInfo, HashState > best_source;

	auto queue_state = [&](State const& state, Cost const& cost, State const* from, Action const& action) {
		auto ret = best_source.insert(std::make_pair(state, StateInfo(cost, from, action)));
		if (cost < ret.first->second.cost) {
			ret.first->second = StateInfo(cost, from, action);
			ret.second = true;
		}
		if (ret.second) {
			todo.insert(std::make_pair(cost, &ret.first->first));
		}
	};


	auto apply_action = [&queue_state, &top, &bottom, &constraints](Action const& action, State const& state, Cost const& cost) {
		//std::cout << "  doing '" << action.to_string() << "'" << std::endl; //DEBUG
		State next_state = state;
		Cost next_cost = cost;
		if (action.type == Action::MoveLeft) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));

			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, 0).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::Roll0;

			//if this is the first stitch, track it:
			if (state.r_next_roll == State::RRollInvalid) {
				assert(bottom.empty() && state.r + 1 == int32_t(top.size()));
				next_state.r_next_needle = action.needle;
				next_state.r_next_roll = State::RRoll2;
			}
		}
		else if (action.type == Action::MoveRight) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));

			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, 0).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::Roll0;

			//if this is the first stitch, track it:
			if (state.l_prev_roll == State::LRollInvalid) {
				assert(bottom.empty() && state.l == 0);
				next_state.l_prev_needle = action.needle;
				next_state.l_prev_roll = State::LRoll2;
			}

		}
		else if (action.type == Action::RollLeft) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));

			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, -1).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::LRoll1;

			//if this is the first stitch, track it on the other side:
			if (state.r_next_roll == State::RRollInvalid) {
				assert(bottom.empty() && state.r + 1 == int32_t(top.size()));
				next_state.r_next_needle = action.needle;
				next_state.r_next_roll = State::RRoll1;
			}
		}
		else if (action.type == Action::RollRight) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));

			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, +1).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::RRoll1;

			//if this is the first stitch, track it on the other side:
			if (state.l_prev_roll == State::LRollInvalid) {
				assert(bottom.empty() && state.l == 0);
				next_state.l_prev_needle = action.needle;
				next_state.l_prev_roll = State::LRoll1;
			}
		}
		else if (action.type == Action::Roll2Left) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			assert(bottom.empty());
			assert(state.l_prev_roll == State::LRoll2 || state.l_prev_roll == State::LRollInvalid);

			//do we add to penalty when stacking? I guess so.
			next_cost.penalty += top[state.l].after_offset_and_roll(action.needle - top[state.l].needle, -2).penalty(constraints.min_free, constraints.max_free);

			next_state.l += 1;
			next_state.l_prev_needle = action.needle;
			next_state.l_prev_roll = State::LRoll2;

			if (state.r_next_roll == State::RRollInvalid) {
				next_state.r_next_needle = action.needle;
				next_state.r_next_roll = State::Roll0;
			}

		}
		else if (action.type == Action::Roll2Right) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			assert(bottom.empty());
			assert(state.r_next_roll == State::RRoll2 || state.r_next_roll == State::RRollInvalid);

			//do we add to penalty when stacking? I guess so.
			next_cost.penalty += top[state.r].after_offset_and_roll(action.needle - top[state.r].needle, 2).penalty(constraints.min_free, constraints.max_free);

			next_state.r -= 1;
			next_state.r_next_needle = action.needle;
			next_state.r_next_roll = State::RRoll2;

			if (state.l_prev_roll == State::RRollInvalid) {
				next_state.l_prev_needle = action.needle;
				next_state.l_prev_roll = State::Roll0;
			}


			//shouldn't ever need to inform other side.
		}
		else {
			assert(0 && "Unhandled action type.");
		}

		queue_state(next_state, next_cost, &state, action);
	};

	auto expand_state = [&](State const& state, Cost const& cost) {
		assert(state.l <= state.r);
		assert(state.r < int32_t(top.size()));

		assert(state.l_prev_roll <= State::Roll0);
		assert(state.r_next_roll >= State::Roll0);

		std::vector< std::pair< int32_t, Action > > offset_action;

		//First, and most important range: what do the current bridges, constraints, and slack allow in terms of racking?
		int32_t min_ofs = -int32_t(constraints.max_racking);
		int32_t max_ofs = int32_t(constraints.max_racking);

		if (bottom.empty() && state.l == 0 && state.r + 1 == int32_t(top.size())) {
			//no bridges to worry about!
			assert(state.l_prev_roll == State::LRollInvalid && state.r_next_roll == State::RRollInvalid);
		}
		else {
			//can't have | ofs + top[l].needle - state.l_prev_needle | > top[l].left_slack
			//want -top[l].left_slack <= ofs + top[l].needle - state.l_prev_needle <= top[l].left_slack
			// -top[l].left_slack - (top[l].needle - state.l_prev_needle) <= ofs <= top[l].left_slack - (top[l].needle - state.l_prev_needle)
			if (top[state.l].left_slack != SlackForNoYarn) {
				min_ofs = std::max(min_ofs, -top[state.l].left_slack - (top[state.l].needle - state.l_prev_needle));
				max_ofs = std::min(max_ofs, top[state.l].left_slack - (top[state.l].needle - state.l_prev_needle));
			}
			if (top[state.r].right_slack != SlackForNoYarn) {
				min_ofs = std::max(min_ofs, -top[state.r].right_slack - (top[state.r].needle - state.r_next_needle));
				max_ofs = std::min(max_ofs, top[state.r].right_slack - (top[state.r].needle - state.r_next_needle));
			}
		}

		{ //"roll2" moves for left stitch:
			int32_t min = min_ofs + top[state.l].needle;
			int32_t max = max_ofs + top[state.l].needle;
			if (state.l_prev_roll == State::LRollInvalid) {
				//can roll2 to ~anywhere~
			}
			else if (state.l_prev_roll == State::LRoll2) {
				min = std::max(min, state.l_prev_needle + (top[state.l].can_stack_left ? 0 : +1));
				max = std::min(max, state.r_next_needle);
			}
			else {
				min = std::numeric_limits< int32_t >::max();
				max = std::numeric_limits< int32_t >::min();
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.l].needle, Action(Action::Roll2Left, needle));
				//apply_action(Action(Action::Roll2Left, needle), state, cost);
			}
		}

		{ //"roll2" moves for right stitch:
			int32_t min = min_ofs + top[state.r].needle;
			int32_t max = max_ofs + top[state.r].needle;
			if (state.r_next_roll == State::LRollInvalid) {
				//can roll2 to ~anywhere~
			}
			else if (state.r_next_roll == State::RRoll2) {
				min = std::max(min, state.l_prev_needle);
				max = std::min(max, state.r_next_needle + (top[state.r].can_stack_right ? 0 : -1));
			}
			else {
				min = std::numeric_limits< int32_t >::max();
				max = std::numeric_limits< int32_t >::min();
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.r].needle, Action(Action::Roll2Right, needle));
				//apply_action(Action(Action::Roll2Right, needle), state, cost);
			}
		}

		{ //"roll" moves for left stitch:
			int32_t min = min_ofs + top[state.l].needle;
			int32_t max = max_ofs + top[state.l].needle;

			//limit based on left stitches:
			if (state.l_prev_roll == State::LRollInvalid) {
				//nothing on the other bed, do whatever!
			}
			else if (state.l_prev_roll == State::LRoll2) {
				//must arrive to the left of l_prev_needle, as it's on the front:
				max = std::min(max, state.l_prev_needle - 1);
			}
			else if (state.l_prev_roll == State::LRoll1) {
				//can arrive to the left of l_prev_needle or stack:
				max = std::min(max, state.l_prev_needle + (top[state.l].can_stack_left ? 0 : -1));
			}
			else {
				assert(state.l_prev_roll == State::Roll0);
				//have already moved a stitch, so can't continue to roll:
				min = std::numeric_limits< int32_t >::max();
				max = std::numeric_limits< int32_t >::min();
			}
			//limit based on right stitches:
			if (state.r_next_roll == State::Roll0) {
				//must be to the left of the top-bed r_prev_needle:
				max = std::min(max, state.r_next_needle - 1);
			}
			else if (state.r_next_roll == State::RRoll2) {
				assert(state.l_prev_roll == State::Roll0); //would need to limit, but already on front bed
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.l].needle, Action(Action::RollLeft, needle));
				//apply_action(Action(Action::RollLeft, needle), state, cost);
			}
		}

		{ //"roll" moves for right stitch:
			int32_t min = min_ofs + top[state.r].needle;
			int32_t max = max_ofs + top[state.r].needle;

			//limit based on right stitches:
			if (state.r_next_roll == State::RRollInvalid) {
				//nothing on the other bed, do whatever!
			}
			else if (state.r_next_roll == State::RRoll2) {
				//must arrive to the right of r_next_needle, as it's on the front:
				min = std::max(min, state.r_next_needle + 1);
			}
			else if (state.r_next_roll == State::RRoll1) {
				//can arrive to the right of r_prev_needle or stack:
				min = std::max(min, state.r_next_needle + (top[state.r].can_stack_right ? 0 : 1));
			}
			else {
				assert(state.r_next_roll == State::Roll0);
				//have already moved a stitch, so can't continue to roll:
				min = std::numeric_limits< int32_t >::max();
				max = std::numeric_limits< int32_t >::min();
			}

			//limit based on left stitches:
			if (state.l_prev_roll == State::Roll0) {
				//must be to the right of the top-bed l_prev_needle:
				min = std::max(min, state.l_prev_needle + 1);
			}
			else if (state.l_prev_roll == State::LRoll2) {
				assert(state.r_next_roll == State::Roll0); //would need to limit, but already on front bed
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.r].needle, Action(Action::RollRight, needle));
				//apply_action(Action(Action::RollRight, needle), state, cost);
			}
		}

		{ //"move" moves for left stitch:
			int32_t min = min_ofs + top[state.l].needle;
			int32_t max = max_ofs + top[state.l].needle;

			//limit based on left stitches:
			if (state.l_prev_roll == State::LRollInvalid) {
				//nothing on the other bed, do whatever!
			}
			else if (state.l_prev_roll == State::LRoll2) {
				//must arrive to the left of l_prev_needle, as it's on the front; but r_next_needle case below will deal with this
			}
			else if (state.l_prev_roll == State::LRoll1) {
				//nothing to interfere on the left, do whatever!
			}
			else {
				assert(state.l_prev_roll == State::Roll0);
				//must place to the right of existing front bed stuff:
				min = std::max(min, state.l_prev_needle + (top[state.l].can_stack_left ? 0 : 1));
			}
			//limit based on right stitches:
			if (state.r_next_roll == State::Roll0) {
				//must be to the left of the top-bed r_prev_needle:
				max = std::min(max, state.r_next_needle - 1);
			}
			else if (state.r_next_roll == State::RRoll2) {
				assert(state.l_prev_roll == State::Roll0); //would need to limit, but already on front bed
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.l].needle, Action(Action::MoveLeft, needle));
				//apply_action(Action(Action::MoveLeft, needle), state, cost);
			}
		}

		{ //"move" moves for right stitch:
			int32_t min = min_ofs + top[state.r].needle;
			int32_t max = max_ofs + top[state.r].needle;

			//limit based on left stitches:
			if (state.r_next_roll == State::RRollInvalid) {
				//nothing on the other bed, do whatever!
			}
			else if (state.r_next_roll == State::RRoll2) {
				//must arrive to the right of r_next_needle, as it's on the front; but l_prev_needle case below will deal with this
			}
			else if (state.r_next_roll == State::RRoll1) {
				//nothing to interfere on the left, do whatever!
			}
			else {
				assert(state.r_next_roll == State::Roll0);
				//must place to the left of existing front bed stuff:
				max = std::min(max, state.r_next_needle + (top[state.r].can_stack_right ? 0 : -1));
			}
			//limit based on left stitches:
			if (state.l_prev_roll == State::Roll0) {
				//must be to the right of the top-bed l_prev_needle:
				min = std::max(min, state.l_prev_needle + 1);
			}
			else if (state.l_prev_roll == State::LRoll2) {
				assert(state.r_next_roll == State::Roll0); //would need to limit, but already on front bed
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				offset_action.emplace_back(needle - top[state.r].needle, Action(Action::MoveRight, needle));
				//apply_action(Action(Action::MoveRight, needle), state, cost);
			}
		}

		std::stable_sort(offset_action.begin(), offset_action.end(), [](std::pair< int32_t, Action > const& a, std::pair< int32_t, Action > const& b) -> bool {
			return std::abs(a.first) < std::abs(b.first);
			});

		for (auto const& oa : offset_action) {
			apply_action(oa.second, state, cost);
		}
	};

	//Initial state(s):
	{
		//NOTE: if bottom is empty, the prev_needle / next_needle fields don't matter until after first xfer
		State init;
		init.l = 0;
		init.l_prev_needle = (bottom.empty() ? top[0].needle : bottom[0].needle);
		init.l_prev_roll = (bottom.empty() ? State::LRollInvalid : State::LRoll1);
		init.r = top.size() - 1;
		init.r_next_needle = (bottom.empty() ? top.back().needle : bottom.back().needle);
		init.r_next_roll = (bottom.empty() ? State::RRollInvalid : State::RRoll1);
		Cost cost;
		cost.penalty = 0;
		queue_state(init, cost, nullptr, Action(Action::None, 0));
	}

	//Actual search:
	const State* best = nullptr;
	while (!todo.empty()) {
		Cost cost = todo.begin()->first;
		const State* state = todo.begin()->second;
		todo.erase(todo.begin());
		{ //see if this is the first time the state has been expanded:
			auto f = best_source.find(*state);
			assert(f != best_source.end());
			assert(&(f->first) == state);
			if (f->second.cost < cost) continue;
			assert(f->second.cost == cost);
		}
		//std::cout << "Considering " << state->to_string() << std::endl; //DEBUG

		//if this is an ending state, end:
		if (state->l > state->r) {
			best = state;
			break;
		}
		//otherwise, expand:
		expand_state(*state, cost);
	}
	assert(best && "Must have gotten to some ending state.");

	//read back operations from best:
	std::vector< Transfer > ops;

	while (best) {
		auto f = best_source.find(*best);
		assert(f != best_source.end());
		if (f->second.source == nullptr) break;
		State const& state = *f->second.source;
		Action const& action = f->second.action;

		if (action.type == Action::MoveLeft) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.l].needle), BedNeedle(to_top_bed, action.needle));
		}
		else if (action.type == Action::MoveRight) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.r].needle), BedNeedle(to_top_bed, action.needle));
		}
		else if (action.type == Action::RollLeft) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.l].needle), BedNeedle(to_bottom_bed, action.needle));
		}
		else if (action.type == Action::RollRight) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.r].needle), BedNeedle(to_bottom_bed, action.needle));
		}
		else if (action.type == Action::Roll2Left) {
			assert(state.l >= 0 && state.l < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.l].needle), BedNeedle(to_top_bed, action.needle));
		}
		else if (action.type == Action::Roll2Right) {
			assert(state.r >= 0 && state.r < int32_t(top.size()));
			ops.emplace_back(BedNeedle(top_bed, top[state.r].needle), BedNeedle(to_top_bed, action.needle));
		}
		else {
			assert(0 && "Invalid action type.");
		}
		ops.back().why = state.to_string() + "; " + action.to_string();

		best = f->second.source;
	}

	std::reverse(ops.begin(), ops.end());

	/*
		std::cout << "  Final plan:\n"; //DEBUG
		for (auto const &op : ops) {
			std::cout << "    " << op.to_string() << '\n';
		}
		std::cout.flush(); //DEBUG
		*/

		/*
			std::cout << "Before Collapse:\n"; //DEBUG
			draw_beds(top_bed, top, bottom_bed, bottom); //DEBUG
		*/
	run_transfers(constraints,
		top_bed, top,
		bottom_bed, bottom,
		ops,
		to_top_bed, &to_top,
		to_bottom_bed, &to_bottom);

	/*
		std::cout << "After Collapse:\n"; //DEBUG
		draw_beds(to_top_bed, to_top, to_bottom_bed, to_bottom); //DEBUG
	*/

	plan.insert(plan.end(), ops.begin(), ops.end());

}


void best_expand(
	Constraints const& constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal >* to_top_,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal >* to_bottom_,
	std::vector< Transfer >* plan_
) {
	//Expand won't change the top bed's location, but will change the bottom bed's:
	assert(top_bed == to_top_bed);
	assert(bottom_bed != to_bottom_bed);

	//make sure all output arrays exist:
	assert(to_top_);
	auto& to_top = *to_top_;
	assert(to_bottom_);
	auto& to_bottom = *to_bottom_;
	assert(plan_);
	auto& plan = *plan_;

	//Clear output:
	to_top.clear();
	to_bottom.clear();
	//NOTE: don't clear plan, only append.

	//------------

	//if no stitches on bottom, nothing to move, so done:
	if (bottom.empty()) {
		to_top = top;
		to_bottom.clear();
		return;
	}

	//Otherwise, Dijkstra-style search for optimal sequence of moves:

	//Moves in 'expand' move a stitch from the bottom bed (or roll a stitch from the top bed) to the to_bottom bed

	//Search state:
#pragma pack(push,1)
	struct State {
		//if 'l' or 'r' is < 0 or >= bottom.size(), this means that some top stitches are getting rolled
		int32_t l; //index of left stitch
		int32_t r; //index of right stitch
		int32_t l_next_needle; //needle that the stitch at index l+1 'l' was moved to
		int32_t r_prev_needle; //needle that the stitch at index r-1 was moved to
		bool operator==(State const& o) const {
			return l == o.l
				&& r == o.r
				&& l_next_needle == o.l_next_needle
				&& r_prev_needle == o.r_prev_needle
				;
		};
		std::string to_string() const {
			return std::to_string(l) + " " + std::to_string(l_next_needle) + "." + std::to_string(r_prev_needle) + " " + std::to_string(r);
		}
	};
#pragma pack(pop)
	static_assert(sizeof(State) == 4 * 4, "expand's State is packed");

	struct HashState {
		size_t operator()(State const& state) const {
			static std::hash< std::string > hash;
			return hash(std::string(reinterpret_cast<char const*>(&state), reinterpret_cast<char const*>(&state) + sizeof(state)));
		}
	};

	//Let's do this in terms of the actions that can be applied:
	struct Action {
		enum Type : uint8_t {
			None,
			MoveLeft,
			MoveRight,
			Finish //charge for all the remaining top stitches
		} type;
		int32_t needle;

		Action(Type type_, int32_t needle_) : type(type_), needle(needle_) { }

		std::string to_string() const {
			if (type == None) return "ENone";
			else if (type == MoveLeft) return "EMoveLeft to " + std::to_string(needle);
			else if (type == MoveRight) return "EMoveRight to " + std::to_string(needle);
			else if (type == Finish) return "EFinish";
			else {
				assert(0 && "invalid move type");
				return "!"; //never reached
			}
		}
	};


	struct Cost {
		uint32_t penalty = 0;
		bool operator<(Cost const& o) const {
			return penalty < o.penalty;
		}
		bool operator==(Cost const& o) const {
			return penalty == o.penalty;
		}
	};

	struct StateInfo {
		Cost cost;
		State const* source;
		Action action;

		StateInfo(Cost const& cost_, State const* source_, Action const& action_) : cost(cost_), source(source_), action(action_) { }
	};

	std::multimap< Cost, const State* > todo;
	std::unordered_map< State, StateInfo, HashState > best_source;

	auto queue_state = [&](State const& state, Cost const& cost, State const* from, Action const& action) {
		auto ret = best_source.insert(std::make_pair(state, StateInfo(cost, from, action)));
		if (cost < ret.first->second.cost) {
			ret.first->second = StateInfo(cost, from, action);
			ret.second = true;
		}
		if (ret.second) {
			todo.insert(std::make_pair(cost, &ret.first->first));
		}
	};


	auto apply_action = [&queue_state, &top, &bottom, &constraints](Action const& action, State const& state, Cost const& cost) {
		//std::cout << "  doing '" << action.to_string() << "'" << std::endl; //DEBUG
		State next_state = state;
		Cost next_cost = cost;

		int32_t top_l = -(state.l + 1);
		int32_t top_r = int32_t(top.size()) - 1 - (state.r - int32_t(bottom.size()));

		assert(state.l <= state.r);
		assert(top_l <= top_r);

		if (action.type == Action::MoveLeft) {
			if (state.l >= 0) {
				//l is on bottom bed
				next_cost.penalty += bottom[state.l].after_offset_and_roll(action.needle - bottom[state.l].needle, 0).penalty(constraints.min_free, constraints.max_free);

				next_state.l -= 1;
				next_state.l_next_needle = action.needle;

				if (state.r == state.l) {
					next_state.r += 1;
					next_state.r_prev_needle = action.needle;
				}
			}
			else {
				assert(state.l < 0);
				//l has wrapped to top bed
				assert(top_l >= 0 && uint32_t(top_l) < top.size());

				next_cost.penalty += top[top_l].after_offset_and_roll(action.needle - top[top_l].needle, -1).penalty(constraints.min_free, constraints.max_free);

				next_state.l -= 1;
				next_state.l_next_needle = action.needle;

				//Ignore top_l == top_r case.
			}
		}
		else if (action.type == Action::MoveRight) {
			if (state.r < int32_t(bottom.size())) {
				//r is on bottom bed
				assert(state.r >= 0);
				next_cost.penalty += bottom[state.r].after_offset_and_roll(action.needle - bottom[state.r].needle, 0).penalty(constraints.min_free, constraints.max_free);

				next_state.r += 1;
				next_state.r_prev_needle = action.needle;

				if (state.l == state.r) {
					next_state.l -= 1;
					next_state.l_next_needle = action.needle;
				}
			}
			else {
				assert(state.r >= int32_t(bottom.size()));
				//r has wrapped to top bed
				assert(top_r >= 0 && uint32_t(top_r) < top.size());

				next_cost.penalty += top[top_r].after_offset_and_roll(action.needle - top[top_r].needle, 1).penalty(constraints.min_free, constraints.max_free);

				next_state.r += 1;
				next_state.r_prev_needle = action.needle;

				//Ignore top_l == top_r case; as it's already a final state.
			}
		}
		else if (action.type == Action::Finish) {
			//got to have done all bottom stitches:
			assert(top_l >= 0 && top_r < int32_t(top.size()));

			//charge for all the not-done top stitches:
			for (int32_t i = top_l; i <= top_r; ++i) {
				next_cost.penalty += top[i].penalty(constraints.min_free, constraints.max_free);
			}

			//cross the indices:
			next_state.l = -1 - int32_t(top.size());
			next_state.r = bottom.size() + top.size();

			assert(-(next_state.l + 1) == int32_t(top.size()));
			assert(int32_t(top.size()) - 1 - (next_state.r - int32_t(bottom.size())) == -1);
			assert(next_state.r - next_state.l >= int32_t(top.size() + bottom.size()));

		}
		else {
			assert(0 && "Invalid action type");
		}
		//std::cout << "     penalty " << next_cost.penalty  << " from " << cost.penalty << std::endl; //DEBUG
		assert(next_cost.penalty >= cost.penalty);

		queue_state(next_state, next_cost, &state, action);
	};

	auto expand_state = [&](State const& state, Cost const& cost) {
		assert(state.l <= state.r);

		int32_t top_l = -(state.l + 1);
		int32_t top_r = int32_t(top.size()) - 1 - (state.r - int32_t(bottom.size()));

		assert(top_l <= top_r); //would have been flagged a final state otherwise

		//First, and most important range: what do the current bridges, constraints, and slack allow in terms of racking?
		int32_t min_ofs = -int32_t(constraints.max_racking);
		int32_t max_ofs = int32_t(constraints.max_racking);

		if (state.l == state.r) {
			//no bridges to worry about!
		}
		else {
			//can't have | ofs + top[l].needle - state.l_prev_needle | > top[l].left_slack
			//want -top[l].left_slack <= ofs + top[l].needle - state.l_prev_needle <= top[l].left_slack
			// -top[l].left_slack - (top[l].needle - state.l_prev_needle) <= ofs <= top[l].left_slack - (top[l].needle - state.l_prev_needle)

			int32_t l_needle = 0;
			Slack l_slack = SlackForNoYarn;
			if (state.l >= 0) {
				assert(state.l >= 0 && state.l < int32_t(bottom.size()));
				l_needle = bottom[state.l].needle;
				l_slack = bottom[state.l].right_slack;
			}
			else if (top_l < int32_t(top.size())) {
				assert(top_l >= 0 && top_l < int32_t(top.size()));
				l_needle = top[top_l].needle;
				l_slack = top[top_l].left_slack;
			}
			else if (top_l == int32_t(top.size())) {
				assert(!bottom.empty());
				l_needle = bottom.back().needle;
				l_slack = bottom.back().right_slack;
			}

			int32_t r_needle = 0;
			Slack r_slack = SlackForNoYarn;
			if (state.r < int32_t(bottom.size())) {
				assert(state.r >= 0 && state.r < int32_t(bottom.size()));
				r_needle = bottom[state.r].needle;
				r_slack = bottom[state.r].left_slack;
			}
			else if (top_r >= 0) {
				assert(top_r >= 0 && top_r < int32_t(top.size()));
				r_needle = top[top_r].needle;
				r_slack = top[top_r].right_slack;
			}
			else if (top_r == -1) {
				assert(!bottom.empty());
				r_needle = bottom[0].needle;
				r_slack = bottom[0].left_slack;
			}

			if (l_slack != SlackForNoYarn) {
				min_ofs = std::max(min_ofs, -l_slack - (l_needle - state.l_next_needle));
				max_ofs = std::min(max_ofs, l_slack - (l_needle - state.l_next_needle));
			}
			if (r_slack != SlackForNoYarn) {
				min_ofs = std::max(min_ofs, -r_slack - (r_needle - state.r_prev_needle));
				max_ofs = std::min(max_ofs, r_slack - (r_needle - state.r_prev_needle));
			}
		}

		//if it is possible to finish, try the finishing move:
		if (state.l < 0 && state.r >= int32_t(bottom.size())) {
			assert(top_l >= 0 && top_r < int32_t(top.size()));
			//only allow finish if zero offset is valid:
			if (min_ofs <= 0 && 0 <= max_ofs) {
				apply_action(Action(Action::Finish, 0), state, cost);
			}
		}

		{ //left moves:
			int32_t min, max;
			if (state.l >= 0) {
				assert(state.l < int32_t(bottom.size()));
				//where can bottom stitch be dropped?

				min = bottom[state.l].needle + min_ofs;
				max = bottom[state.l].needle + max_ofs;
				if (state.l == state.r) {
					//first move, no limit.
				}
				else {
					//NOTE: mis-use of 'can_stack_' in asymmetric case; probably will use stacking priorities for this anyway,however.
					//must place to the left of last placed stitch:
					max = std::min(max, state.l_next_needle + (bottom[state.l].can_stack_right ? 0 : -1));
				}
			}
			else if (top_l >= int32_t(top.size())) {
				//if we've somehow rolled *all* of the top stitches, not much to be done:
				max = std::numeric_limits< int32_t >::min();
				min = std::numeric_limits< int32_t >::max();
			}
			else {
				assert(top_l >= 0 && top_l < int32_t(top.size()));

				min = top[top_l].needle + min_ofs;
				max = top[top_l].needle + max_ofs;

				//must place left of last-placed stitch:
				//NOTE: mis-use of 'can_stack_' in asymmetric case; probably will use stacking priorities for this anyway,however.
				max = std::min(max, state.l_next_needle + (top[top_l].can_stack_left ? 0 : -1));

				if (state.r < int32_t(bottom.size())) {
					//can't move at all if pinned by bottom stitch:
					bool pinned = false;
					//    l --- o
					//  r ...
					if (bottom[state.r].needle <= top[top_l].needle) {
						pinned = true;
					}

					//  l --- o
					//     r ...
					if (top_l + 1 < int32_t(top.size()) && top[top_l + 1].needle > bottom[state.r].needle) {
						pinned = true;
					}

					//  l ------.
					//     r -- o
					if (top_l + 1 == int32_t(top.size()) && state.r + 1 < int32_t(bottom.size())) {
						pinned = true;
					}

					if (pinned) {
						max = std::numeric_limits< int32_t >::min();
						min = std::numeric_limits< int32_t >::max();
					}
				}
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				apply_action(Action(Action::MoveLeft, needle), state, cost);
			}

		}

		{ //right moves:
			int32_t min, max;
			if (state.r < int32_t(bottom.size())) {
				assert(state.r >= 0);
				//where can bottom stitch be dropped?

				min = bottom[state.r].needle + min_ofs;
				max = bottom[state.r].needle + max_ofs;
				if (state.l == state.r) {
					//first move, no limit.
				}
				else {
					//NOTE: mis-use of 'can_stack_' in asymmetric case; probably will use stacking priorities for this anyway,however.
					//must place to the right of last placed stitch:
					min = std::max(min, state.r_prev_needle + (bottom[state.r].can_stack_left ? 0 : +1));
				}
			}
			else if (top_r < 0) {
				//if we've somehow rolled *all* of the top stitches, not much to be done
				max = std::numeric_limits< int32_t >::min();
				min = std::numeric_limits< int32_t >::max();
			}
			else {
				assert(top_r >= 0 && top_r < int32_t(top.size()));

				min = top[top_r].needle + min_ofs;
				max = top[top_r].needle + max_ofs;

				//must place right of last-placed stitch:
				//NOTE: mis-use of 'can_stack_' in asymmetric case; probably will use stacking priorities for this anyway,however.
				min = std::max(min, state.r_prev_needle + (top[top_r].can_stack_right ? 0 : +1));

				if (state.l >= 0) {
					//can't move at all if pinned by bottom stitch:
					bool pinned = false;
					//  o --- r
					//     .... l
					if (bottom[state.l].needle >= top[top_r].needle) {
						pinned = true;
					}

					//  o --- r
					// ... l
					if (top_r - 1 >= 0 && top[top_r - 1].needle < bottom[state.l].needle) {
						pinned = true;
					}
					// /--- r
					// o ... l
					if (top_r == 0 && state.l > 0) {
						pinned = true;
					}

					if (pinned) {
						max = std::numeric_limits< int32_t >::min();
						min = std::numeric_limits< int32_t >::max();
					}
				}
			}

			min = std::max(min, constraints.min_free);
			max = std::min(max, constraints.max_free);

			for (int32_t needle = min; needle <= max; ++needle) {
				apply_action(Action(Action::MoveRight, needle), state, cost);
			}

		}

	};


	//Initial states:
	for (int32_t i = 0; uint32_t(i) < bottom.size(); ++i) {
		State init;
		init.l = i;
		init.l_next_needle = 0;
		init.r = i;
		init.r_prev_needle = 0;
		Cost cost;
		cost.penalty = 0;
		queue_state(init, cost, nullptr, Action(Action::None, 0));
	}
	//TODO: consider states that start on *top* needles (seems like a rare and unneeded case)

	//Actual search:
	const State* best = nullptr;
	while (!todo.empty()) {
		Cost cost = todo.begin()->first;
		const State* state = todo.begin()->second;
		todo.erase(todo.begin());
		{ //see if this is the first time the state has been expanded:
			auto f = best_source.find(*state);
			assert(f != best_source.end());
			assert(&(f->first) == state);
			if (f->second.cost < cost) continue;
			assert(f->second.cost == cost);
		}
		//std::cout << "Considering " << state->to_string() << "  [penalty: " << cost.penalty << "]" << std::endl; //DEBUG

		//if this is an ending state, end:
		if (state->l < 0 && state->r >= int32_t(bottom.size()) && state->r - state->l > int32_t(top.size() + bottom.size())) {
			best = state;
			break;
		}
		//otherwise, expand:
		expand_state(*state, cost);
	}
	assert(best && "Must have gotten to some ending state.");

	//read back operations from best:
	std::vector< Transfer > ops;

	bool is_first = true;
	while (best) {
		auto f = best_source.find(*best);
		assert(f != best_source.end());
		if (f->second.source == nullptr) break;
		State const& state = *f->second.source;
		Action const& action = f->second.action;

		int32_t top_l = -(state.l + 1);
		int32_t top_r = int32_t(top.size()) - 1 - (state.r - int32_t(bottom.size()));

		assert(state.l <= state.r);
		assert(top_l <= top_r);

		if (action.type == Action::MoveLeft) {
			if (state.l >= 0) {
				ops.emplace_back(BedNeedle(bottom_bed, bottom[state.l].needle), BedNeedle(to_bottom_bed, action.needle));
			}
			else {
				assert(state.l < 0);
				assert(top_l >= 0 && uint32_t(top_l) < top.size());
				ops.emplace_back(BedNeedle(top_bed, top[top_l].needle), BedNeedle(to_bottom_bed, action.needle));
			}
		}
		else if (action.type == Action::MoveRight) {
			if (state.r < int32_t(bottom.size())) {
				ops.emplace_back(BedNeedle(bottom_bed, bottom[state.r].needle), BedNeedle(to_bottom_bed, action.needle));
			}
			else {
				assert(state.r >= int32_t(bottom.size()));
				assert(top_r >= 0 && uint32_t(top_r) < top.size());
				ops.emplace_back(BedNeedle(top_bed, top[top_r].needle), BedNeedle(to_bottom_bed, action.needle));
			}
		}
		else if (action.type == Action::Finish) {
			assert(is_first);
		}
		else {
			assert(0 && "Invalid action type.");
		}
		if (action.type != Action::Finish) {
			assert(!ops.empty() && ops.back().why == "");
			ops.back().why = state.to_string() + "; " + action.to_string();
		}

		best = f->second.source;
		is_first = false;
	}
	std::cout.flush(); //DEBUG

	std::reverse(ops.begin(), ops.end());

	/*
		std::cout << "  Final plan:\n"; //DEBUG
		for (auto const &op : ops) {
			std::cout << "    " << op.to_string() << '\n';
		}
		std::cout.flush(); //DEBUG
	*/

	/*
		std::cout << "Before Expand:\n"; //DEBUG
		draw_beds(top_bed, top, bottom_bed, bottom); //DEBUG
	*/
	run_transfers(constraints,
		top_bed, top,
		bottom_bed, bottom,
		ops,
		to_top_bed, &to_top,
		to_bottom_bed, &to_bottom);
	/*
		std::cout << "After Expand:\n"; //DEBUG
		draw_beds(to_top_bed, to_top, to_bottom_bed, to_bottom); //DEBUG
	*/

	plan.insert(plan.end(), ops.begin(), ops.end());


}



void best_shift(
	Constraints const& constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal >* to_top_,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal >* to_bottom_,
	std::vector< Transfer >* plan_
) {
	//Shift will change both beds:
	assert(top_bed != to_top_bed);
	assert(bottom_bed != to_bottom_bed);

	//make sure all output arrays exist:
	assert(to_top_);
	auto& to_top = *to_top_;
	assert(to_bottom_);
	auto& to_bottom = *to_bottom_;
	assert(plan_);
	auto& plan = *plan_;

	//Clear output:
	to_top.clear();
	to_bottom.clear();
	//NOTE: don't clear plan, only append.

	//------------

	int32_t best_ofs = 0;
	uint32_t best_penalty = std::numeric_limits< uint32_t >::max();

	//Find best offset to shift to:

	auto do_ofs = [&](int32_t ofs) {
		//would shifting to this offset move outside limits?
		if (!top.empty()) {
			if (top[0].needle + ofs < constraints.min_free) return;
			if (top.back().needle + ofs > constraints.max_free) return;
		}
		if (!bottom.empty()) {
			if (bottom[0].needle + ofs < constraints.min_free) return;
			if (bottom.back().needle + ofs > constraints.max_free) return;
		}
		uint32_t penalty = 0;
		for (auto const& nrg : top) {
			penalty += nrg.after_offset_and_roll(ofs, 0).penalty(constraints.min_free, constraints.max_free);
		}
		for (auto const& nrg : bottom) {
			penalty += nrg.after_offset_and_roll(ofs, 0).penalty(constraints.min_free, constraints.max_free);
		}
		if (penalty < best_penalty) {
			best_penalty = penalty;
			best_ofs = ofs;
		}
		//std::cout << "  offset " << ofs << " has penalty " << penalty << std::endl; //DEBUG
	};

	for (int32_t ofs = 0; ofs < int32_t(constraints.max_racking); ++ofs) {
		do_ofs(ofs);
		if (ofs != 0) do_ofs(-ofs);
	}
	//std::cout << "Best offset: " << best_ofs << std::endl; //DEBUG

	assert(best_penalty < std::numeric_limits< uint32_t >::max());

	//shift by that offset:

	std::vector< Transfer > ops;
	for (auto const& nrg : top) {
		ops.emplace_back(BedNeedle(top_bed, nrg.needle), BedNeedle(to_top_bed, nrg.needle + best_ofs));
	}
	for (auto const& nrg : bottom) {
		ops.emplace_back(BedNeedle(bottom_bed, nrg.needle), BedNeedle(to_bottom_bed, nrg.needle + best_ofs));
	}

	/*
		std::cout << "Before Shift:\n"; //DEBUG
		draw_beds(top_bed, top, bottom_bed, bottom); //DEBUG
	*/

	run_transfers(constraints,
		top_bed, top,
		bottom_bed, bottom,
		ops,
		to_top_bed, &to_top,
		to_bottom_bed, &to_bottom
	);

	/*
		std::cout << "After Shift:\n"; //DEBUG
		draw_beds(to_top_bed, to_top, to_bottom_bed, to_bottom); //DEBUG
	*/

	plan.insert(plan.end(), ops.begin(), ops.end());


}



void draw_beds(
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom
) {
	//Should make something like:
	//         -10       -9        -8          -7     <-- needle #
	// bs:    [9+2>-2-[9+2]-------1------[            <-- blob per nrg, line per edge w/ slack
	// f :     ...


	int32_t min_needle = std::numeric_limits< int32_t >::max();
	int32_t max_needle = std::numeric_limits< int32_t >::min();
	for (auto const& nrg : top) {
		min_needle = std::min(min_needle, nrg.needle);
		max_needle = std::max(max_needle, nrg.needle);
	}
	for (auto const& nrg : bottom) {
		min_needle = std::min(min_needle, nrg.needle);
		max_needle = std::max(max_needle, nrg.needle);
	}

	auto make_label = [](NeedleRollGoal const& nrg) -> std::string {
		std::string ret = "";
		if (nrg.left_slack != SlackForNoYarn) ret += std::to_string(nrg.left_slack);
		ret += (nrg.can_stack_left ? "<" : "[");
		ret += std::to_string(nrg.goal);
		if (nrg.roll > 0) ret += "+" + std::to_string(nrg.roll);
		if (nrg.roll < 0) ret += std::to_string(nrg.roll);
		ret += (nrg.can_stack_right ? ">" : "]");
		if (nrg.right_slack != SlackForNoYarn) ret += std::to_string(nrg.right_slack);
		return ret;
	};

	std::vector< std::string > needle_labels(max_needle - min_needle + 1, "");
	std::vector< std::string > top_labels(max_needle - min_needle + 1, "");
	std::vector< std::string > bottom_labels(max_needle - min_needle + 1, "");

	for (int32_t n = min_needle; n <= max_needle; ++n) {
		needle_labels[n - min_needle] = std::to_string(n);
	}

	for (auto const& s : top) {
		assert(top_labels[s.needle - min_needle] == "");
		top_labels[s.needle - min_needle] = make_label(s);
	}
	for (auto const& s : bottom) {
		assert(bottom_labels[s.needle - min_needle] == "");
		bottom_labels[s.needle - min_needle] = make_label(s);
	}

	int32_t column_width = 1;
	for (auto const& l : needle_labels) column_width = std::max< int32_t >(column_width, l.size());
	for (auto const& l : top_labels) column_width = std::max< int32_t >(column_width, l.size());
	for (auto const& l : bottom_labels) column_width = std::max< int32_t >(column_width, l.size());
	column_width += 1;

	std::ostringstream needle_str, top_str, bottom_str;

	{ //row labels:
		auto make_row_header = [](BedNeedle::Bed bed) {
			if (bed == BedNeedle::Front) return "f :";
			else if (bed == BedNeedle::FrontSliders) return "fs:";
			else if (bed == BedNeedle::BackSliders) return "bs:";
			else if (bed == BedNeedle::Back) return "b :";
			else {
				assert(0 && "Doesn't work.");
				return "!"; //never executed
			}
		};
		needle_str << "   ";
		top_str << make_row_header(top_bed);
		bottom_str << make_row_header(bottom_bed);
	}

	auto pad = [](std::string str, char p, uint32_t width) -> std::string {
		while (str.size() < width) {
			str += p;
			if (str.size() < width) str = p + str;
		}
		return str;
	};

	for (auto& l : needle_labels) l = pad(l, ' ', column_width);
	for (auto& l : top_labels) l = pad(l, ' ', column_width);
	for (auto& l : bottom_labels) l = pad(l, ' ', column_width);

	/*
		//Someday, nice edge labels
		for (uint32_t i = 0; i + 1 < top.size(); ++i) {
			if (top[i].needle == top[i+1].needle) continue;
			std::string label = std::to_string(top[i].right_slack);
			uint32_t space = 0;
			for (uint32_t n = top[i].needle; n <= top[i+1].needle; ++n) {
				if (n == top[i].needle) {
				} else if (n == top[i+1].needle) {
				} else {
					assert(top_labels[i]
				}
			}
			assert(top_labels[top[i].needle - min_needle] != "");
			uint32_t space =
			if (top[i].
		}
	*/

	for (int32_t n = min_needle; n <= max_needle; ++n) {
		needle_str << needle_labels[n - min_needle];
		top_str << top_labels[n - min_needle];
		bottom_str << bottom_labels[n - min_needle];
	}

	std::cout << needle_str.str() << '\n';
	std::cout << top_str.str() << '\n';
	std::cout << bottom_str.str() << '\n';
	std::cout.flush();
}

void minimize_winding(std::vector< int32_t >* winding_) {
	assert(winding_);
	auto& winding = *winding_;

	if (winding.empty()) return;
	//can be minimized by subtracting one of the multiples of two that bracket the median:
	std::vector< int32_t > temp = winding;

	//DEBUG: dump numbers
	//for (auto t : temp) std::cout << ' ' << t;
	//std::cout << '\n';

	std::sort(temp.begin(), temp.end());
	int32_t twice_median;
	if (temp.size() == 1) twice_median = 2 * temp[0];
	else twice_median = temp[temp.size() / 2] + temp[(temp.size() + 1) / 2];

	int32_t close_multiple_a = (twice_median / 4) * 2;
	int32_t close_multiple_b = ((twice_median + (twice_median < 0 ? -4 : 4)) / 4) * 2;


	int32_t sum_a = 0;
	int32_t sum_b = 0;
	for (auto& w : winding) {
		sum_a += std::abs(w - close_multiple_a);
		sum_b += std::abs(w - close_multiple_b);
	}
	int32_t close_multiple = (sum_a <= sum_b ? close_multiple_a : close_multiple_b);

	//std::cout << "close multiples: " << close_multiple_a << " (sum:" << sum_a << "), " << close_multiple_b << " (sum: " << sum_b << ")" << std::endl;

	int32_t DEBUG_before = 0;
	int32_t DEBUG_after = 0;
	int32_t DEBUG_after_minus = 0;
	int32_t DEBUG_after_plus = 0;

	for (auto& w : winding) {
		DEBUG_before += std::abs(w);
		w -= close_multiple;
		DEBUG_after += std::abs(w);
		DEBUG_after_minus += std::abs(w - 2);
		DEBUG_after_plus += std::abs(w + 2);
	}
	assert(DEBUG_after <= DEBUG_before);
	assert(DEBUG_after <= DEBUG_after_minus);
	assert(DEBUG_after <= DEBUG_after_plus);
}

//!!NOTE: transfers don't know if they are rolling or not. This might be an issue!

void run_transfers(
	Constraints const& constraints,
	BedNeedle::Bed top_bed, std::vector< NeedleRollGoal > const& top,
	BedNeedle::Bed bottom_bed, std::vector< NeedleRollGoal > const& bottom,
	std::vector< Transfer > const& plan,
	BedNeedle::Bed to_top_bed, std::vector< NeedleRollGoal >* to_top_,
	BedNeedle::Bed to_bottom_bed, std::vector< NeedleRollGoal >* to_bottom_
) {
	//make sure all output arrays exist:
	assert(to_top_);
	auto& to_top = *to_top_;
	assert(to_bottom_);
	auto& to_bottom = *to_bottom_;

	//clear output arrays:
	to_top.clear();
	to_bottom.clear();

	//no stitches -> nothing to do!
	if (top.size() + bottom.size() == 0) {
		assert(plan.empty());
		return;
	}


	//transform stitches into ccw list:
	//NOTE: this is ccw when viewed with 'top' bed above 'bottom' bed; this can actually be cw in practice, but what matters is that this computation is consistent with how stitches are extracted later.
	std::vector< std::pair< BedNeedle::Bed, NeedleRollGoal > > ccw;
	ccw.reserve(bottom.size() + top.size());
	for (auto si = bottom.begin(); si != bottom.end(); ++si) {
		ccw.emplace_back(std::make_pair(bottom_bed, *si));
	}
	for (auto si = top.rbegin(); si != top.rend(); ++si) {
		ccw.emplace_back(std::make_pair(top_bed, *si));
	}
	assert(ccw.size() == bottom.size() + top.size());

	/*//DEBUG:
	std::cout << "before:";
	for (auto const &bs : ccw) {
		std::cout << ' ' << BedNeedle(bs.first, bs.second.needle).to_string();
	}
	std::cout << "\n";*/

	//everything must be on the 'from' beds:
	for (auto const& bs : ccw) {
		assert(bs.first == bottom_bed || bs.first == top_bed);
	}

	//run transfers on the ccw list:
	for (auto const& t : plan) {
		assert(t.from.bed == top_bed || t.from.bed == bottom_bed);
		assert(t.to.bed == to_top_bed || t.to.bed == to_bottom_bed);

		assert(constraints.min_free <= t.to.needle && t.to.needle <= constraints.max_free); //make sure the transfer goes to a valid needle

		//POTENTIAL OPTIMIZATION: use a hash table to avoid visiting every stitch here:
		uint32_t source = -1U;
		uint32_t target = -1U;
		for (auto& bs : ccw) {
			if (bs.first == t.to.bed && bs.second.needle == t.to.needle) {
				assert(target == -1U);
				target = &bs - &ccw[0];
			}
			if (bs.first == t.from.bed && bs.second.needle == t.from.needle) {
				assert(source == -1U);
				source = &bs - &ccw[0];
			}
		}
		assert(source != -1U);

		auto& bs = ccw[source];
		{ //transform source:
			bs.first = t.to.bed;
			bool from_is_top = (t.from.bed == top_bed);
			bool to_is_top = (t.to.bed == to_top_bed);

			//NOTE: will recompute roll, so only parity matters:
			int32_t roll = (from_is_top == to_is_top ? 0 : 1);
			bs.second = bs.second.after_offset_and_roll(t.to.needle - t.from.needle, roll);
			assert(bs.second.needle == t.to.needle);
		}

		if (target != -1U) {
			auto& tbs = ccw[target];
			assert(tbs.first == bs.first && tbs.second.needle == bs.second.needle);

			//can only stack things with the same eventual goal:
			assert(bs.second.has_same_real_goal_as(tbs.second));

			//handle stacking (merge slack / can_stack):
			if ((source + 1) % ccw.size() == target) {
				if (bs.first == to_top_bed) {
					//target is left of source, and source is stacking atop
					assert(bs.second.can_stack_left);

					bs.second.can_stack_left = tbs.second.can_stack_left;
					bs.second.left_slack = tbs.second.left_slack;
					tbs.second.can_stack_right = bs.second.can_stack_right;
					tbs.second.right_slack = bs.second.right_slack;
				}
				else {
					assert(bs.first == to_bottom_bed);
					//target is right of source, and source is stacking under
					assert(tbs.second.can_stack_left);

					tbs.second.can_stack_left = bs.second.can_stack_left;
					tbs.second.left_slack = bs.second.left_slack;
					bs.second.can_stack_right = tbs.second.can_stack_right;
					bs.second.right_slack = tbs.second.right_slack;
				}
			}
			else {
				assert((target + 1) % ccw.size() == source);
				if (bs.first == to_top_bed) {
					//target is right of source, and source is stacking atop
					assert(bs.second.can_stack_right);

					tbs.second.can_stack_left = bs.second.can_stack_left;
					tbs.second.left_slack = bs.second.left_slack;
					bs.second.can_stack_right = tbs.second.can_stack_right;
					bs.second.right_slack = tbs.second.right_slack;

				}
				else {
					assert(bs.first == to_bottom_bed);
					//target is left of source, and source is stacking under
					assert(tbs.second.can_stack_right);

					bs.second.can_stack_left = tbs.second.can_stack_left;
					bs.second.left_slack = tbs.second.left_slack;
					tbs.second.can_stack_right = bs.second.can_stack_right;
					tbs.second.right_slack = bs.second.right_slack;
				}
			}

			//TODO: eventually track how many stitches are on this location with 'weight':
			//bs.weight = tbs.weight = bs.weight + tbs.weight;

			//delete source or target stitch:
			ccw.erase(ccw.begin() + std::max(source, target));
		}
	}

	/*//DEBUG:
	std::cout << "after:";
	for (auto const &bs : ccw) {
		std::cout << ' ' << BedNeedle(bs.first, bs.second.needle).to_string();
	}
	std::cout << "\n";*/

	//must have placed everything on the 'to' beds:
	for (auto const& bs : ccw) {
		assert(bs.first == to_bottom_bed || bs.first == to_top_bed);
	}

	{ //update roll:
		auto swaps = [&to_bottom_bed, &to_top_bed](BedNeedle const& p, BedNeedle const& n) -> int8_t {
			if (p.bed == n.bed) {
				if (p.bed == to_bottom_bed) {
					return (p.needle <= n.needle ? 0 : 2);
				}
				else {
					assert(p.bed == to_top_bed);
					return (p.needle >= n.needle ? 0 : 2);
				}
			}
			else {
				return 1;
			}
		};

		auto from_bn = [](std::pair< BedNeedle::Bed, NeedleRollGoal > const& bs) -> BedNeedle {
			return BedNeedle(bs.first, bs.second.needle);
		};
		auto to_bn = [&to_bottom_bed, &to_top_bed](std::pair< BedNeedle::Bed, NeedleRollGoal > const& bs) -> BedNeedle {
			if ((bs.first == to_bottom_bed) == (bs.second.roll % 2 == 0)) {
				//on bottom bed and not rolling, or on top bed and rolling
				return BedNeedle(to_bottom_bed, bs.second.goal);
			}
			else {
				//on top bed and not rolling, or on bottom bed and rolling
				return BedNeedle(to_top_bed, bs.second.goal);
			}
		};

		std::vector< int32_t > winding;
		winding.reserve(ccw.size());
		winding.emplace_back(from_bn(ccw[0]).bed == to_bn(ccw[0]).bed ? 0 : 1);
		for (uint32_t i = 1; i < ccw.size(); ++i) {
			winding.emplace_back(
				winding.back()
				- swaps(from_bn(ccw[i - 1]), from_bn(ccw[i]))
				+ swaps(to_bn(ccw[i - 1]), to_bn(ccw[i]))
			);
		}
		//make sure winding is consistent:
		assert((winding.back()
			- swaps(from_bn(ccw.back()), from_bn(ccw[0]))
			+ swaps(to_bn(ccw.back()), to_bn(ccw[0]))
			- winding[0]) % 2 == 0);

		minimize_winding(&winding);

		assert(winding.size() == ccw.size());
		for (uint32_t i = 0; i < ccw.size(); ++i) {
			//make sure winding minimization hasn't changed meaning of roll:
			int32_t new_roll = (ccw[i].first == to_bottom_bed ? winding[i] : -winding[i]);
			assert((ccw[i].second.roll - new_roll) % 2 == 0);
			ccw[i].second.roll = new_roll;
		}
	}


	//transform ccw list back into beds:
	{
		uint32_t first_stitch = 0; //bottom-most, ccw-most stitch
		for (uint32_t i = 1; i < ccw.size(); ++i) {
			if (ccw[i].first == to_bottom_bed) {
				if (ccw[first_stitch].first == to_top_bed || ccw[first_stitch].second.needle > ccw[i].second.needle) {
					first_stitch = i;
				}
			}
			else {
				assert(ccw[i].first == to_top_bed);
				if (ccw[first_stitch].first == to_top_bed && ccw[first_stitch].second.needle < ccw[i].second.needle) {
					first_stitch = i;
				}
			}
		}
		std::rotate(ccw.begin(), ccw.begin() + first_stitch, ccw.end());

		//ccw now looks like either bbbbbbtttttt [or just ttttt]
		auto ci = ccw.begin();

		to_bottom.clear();
		while (ci != ccw.end() && ci->first == to_bottom_bed) {
			to_bottom.emplace_back(ci->second);
			++ci;
		}

		to_top.clear();
		while (ci != ccw.end() && ci->first == to_top_bed) {
			to_top.emplace_back(ci->second);
			++ci;
		}

		assert(ci == ccw.end());

		std::reverse(to_top.begin(), to_top.end());
	}

	//PARANOIA: did all the slacks get maintained properly?
	bool good = true;

	for (uint32_t i = 0; i + 1 < to_top.size(); ++i) {
		good = good && (to_top[i].right_slack == to_top[i + 1].left_slack);
	}
	for (uint32_t i = 0; i + 1 < to_bottom.size(); ++i) {
		good = good && (to_bottom[i].right_slack == to_bottom[i + 1].left_slack);
	}
	if (!to_top.empty() && !to_bottom.empty()) {
		good = good && (to_top[0].left_slack == to_bottom[0].left_slack);
		good = good && (to_top.back().right_slack == to_bottom.back().right_slack);
	}
	else if (!to_top.empty() && to_bottom.empty()) {
		good = good && (to_top[0].left_slack == to_top.back().right_slack);
	}
	else if (to_top.empty() && !to_bottom.empty()) {
		good = good && (to_bottom[0].left_slack == to_bottom.back().right_slack);
	}

	if (!good) {
		std::cout << "!!!!! bad slacks after run_transfers !!!!!\n";
		draw_beds(to_top_bed, to_top, to_bottom_bed, to_bottom);
		exit(1);
	}

}
