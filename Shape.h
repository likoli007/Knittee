#pragma once

#include <vector>
#include <cstdint>
#include <cassert>
#include <limits>

/*
*	Helper class copied over from autoknit, used to represent the shape of a cycle of stitches on the bed.
*/


//Shape represents the shape of a cycle of stitches on the bed.
// It is a map from stitch indices to needles, such that each stitch
//is on a needle adjacent to its index-wise neighbors.

//For a given count, shapes are controlled by a type and offset:
//Valid Shapes:
//   . . o o o o . . .
//   . . o o o o . . .

//   . . . o o . . . .
//   . . o o o o . . .

//   . . . o o o . . .
//   . . o o o . . . .

//   . . o o o . . . .
//   . . . o o o . . .

//   . . o o o o . . .
//   . . . o o . . . .

//   . . o o o o . . .
//   . . o o o . . . .

//   . . o o o . . . .
//   . . o o o o . . .

//   . . . o o o . . .
//   . . o o o o . . .

//   . . o o o o . . .
//   . . . o o o . . .

typedef uint32_t PackedShape;
struct Shape {
	uint32_t roll = 0; //cycle starts at front left at roll zero; roll is ccw steps.

	//shape will be evenly balanced between the beds, with some nibbles taken out of both sides:
	//NOTE: will never have both front + back nibble on a side.
	//NOTE: #nibbles + #stitches must be even
	enum Nibble : uint8_t {
		BackLeft = (1 << 0),
		BackRight = (1 << 1),
		FrontLeft = (1 << 2),
		FrontRight = (1 << 3),
	};
	uint8_t nibbles = 0;

	//------------------

	Shape(uint32_t roll_, uint8_t nibbles_) : roll(roll_), nibbles(nibbles_) {
	}

	void size_to_range(uint32_t size, int32_t* front_min_, int32_t* front_max_, int32_t* back_min_, int32_t* back_max_, int32_t left = 0) const {

		uint32_t width = size
			+ ((nibbles & BackLeft) ? 1 : 0)
			+ ((nibbles & FrontLeft) ? 1 : 0)
			+ ((nibbles & BackRight) ? 1 : 0)
			+ ((nibbles & FrontRight) ? 1 : 0);
		assert(width % 2 == 0);
		width /= 2;

		uint32_t on_front = width
			- ((nibbles & FrontLeft) ? 1 : 0)
			- ((nibbles & FrontRight) ? 1 : 0);
		uint32_t on_back = width
			- ((nibbles & BackLeft) ? 1 : 0)
			- ((nibbles & BackRight) ? 1 : 0);
		assert(on_front + on_back == size);

		int32_t front_left, back_left;

		if (on_back && on_front) {
			front_left = ((nibbles & FrontLeft) ? 1 : 0);
			back_left = ((nibbles & BackLeft) ? 1 : 0);
		}
		else {
			front_left = back_left = 0;
		}


		if (on_front) {
			if (front_min_) *front_min_ = front_left + left;
			if (front_max_) *front_max_ = front_left + int32_t(on_front) - 1 + left;
		}
		else {
			if (front_min_) *front_min_ = std::numeric_limits< int32_t >::max();
			if (front_max_) *front_max_ = std::numeric_limits< int32_t >::min();
		}
		if (on_back) {
			if (back_min_) *back_min_ = back_left + left;
			if (back_max_) *back_max_ = back_left + int32_t(on_back) - 1 + left;
		}
		else {
			if (back_min_) *back_min_ = std::numeric_limits< int32_t >::max();
			if (back_max_) *back_max_ = std::numeric_limits< int32_t >::min();
		}

		//TODO: add some PARANOIA checks to make sure this lines up with append_to_beds...
	}
	//------------------
	void size_index_to_bed_needle(uint32_t size, uint32_t index, char* bed_, int32_t* needle_) const {
		assert(index < size);

		assert(bed_);
		auto& bed = *bed_;

		assert(needle_);
		auto& needle = *needle_;

		//figure out how many items go on the back and how many on the front:
		uint32_t width = size
			+ ((nibbles & BackLeft) ? 1 : 0)
			+ ((nibbles & FrontLeft) ? 1 : 0)
			+ ((nibbles & BackRight) ? 1 : 0)
			+ ((nibbles & FrontRight) ? 1 : 0);
		assert(width % 2 == 0);
		width /= 2;

		uint32_t on_front = width
			- ((nibbles & FrontLeft) ? 1 : 0)
			- ((nibbles & FrontRight) ? 1 : 0);
		uint32_t on_back = width
			- ((nibbles & BackLeft) ? 1 : 0)
			- ((nibbles & BackRight) ? 1 : 0);
		assert(on_front + on_back == size);

		index = (index + roll) % size;

		int32_t front_left, back_left;

		if (on_back && on_front) {
			front_left = ((nibbles & FrontLeft) ? 1 : 0);
			back_left = ((nibbles & BackLeft) ? 1 : 0);
		}
		else {
			front_left = back_left = 0;
		}

		if (index < on_front) {
			bed = 'f';
			needle = front_left + index;
		}
		else {
			assert(index < size);
			bed = 'b';
			index -= on_front;
			assert(index < on_back);
			needle = back_left + on_back - 1 - index;
		}
	}

	//------------------
	//append_to_beds will position the current shape as-leftward-as-possible:
	template< typename C, typename T >
	void append_to_beds(C const& data, T const& gap, std::vector< T >* front_, std::vector< T >* back_) const {
		assert(front_);
		auto& front = *front_;
		assert(back_);
		auto& back = *back_;

		//figure out how many items go on the back and how many on the front:
		uint32_t width = data.size()
			+ ((nibbles & BackLeft) ? 1 : 0)
			+ ((nibbles & FrontLeft) ? 1 : 0)
			+ ((nibbles & BackRight) ? 1 : 0)
			+ ((nibbles & FrontRight) ? 1 : 0);
		assert(width % 2 == 0);
		width /= 2;

		uint32_t on_front = width
			- ((nibbles & FrontLeft) ? 1 : 0)
			- ((nibbles & FrontRight) ? 1 : 0);
		uint32_t on_back = width
			- ((nibbles & BackLeft) ? 1 : 0)
			- ((nibbles & BackRight) ? 1 : 0);
		assert(on_front + on_back == data.size());

		if (on_back && on_front) {
			//add gaps as needed to get beds evened up (taking into account nibbles):
			uint32_t front_add = ((nibbles & BackLeft) ? 1 : 0);
			uint32_t back_add = ((nibbles & FrontLeft) ? 1 : 0);
			while (front.size() + front_add < back.size() + back_add) front.emplace_back(gap);
			while (back.size() + back_add < front.size() + front_add) back.emplace_back(gap);
			if (nibbles & BackLeft) {
				assert(!(nibbles & FrontLeft));
				assert(front.size() + 1 == back.size());
			}
			else if (nibbles & FrontLeft) {
				assert(front.size() == back.size() + 1);
			}
			else {
				assert(front.size() == back.size());
			}
		}

		//starting element in 'data' from roll; cw roll would make this much easier, but so be it:
		uint32_t start = (-int32_t(roll) % int32_t(data.size())) + int32_t(data.size());

		for (uint32_t i = 0; i < on_front; ++i) {
			front.emplace_back(data[(start + i) % data.size()]);
		}
		for (uint32_t i = on_back - 1; i < on_back; --i) {
			back.emplace_back(data[(start + on_front + i) % data.size()]);
		}

		if (on_front && on_back) {
			//make sure beds are as (un-)even as expected:
			assert(front.size() + ((nibbles & FrontRight) ? 1 : 0) == back.size() + ((nibbles & BackRight) ? 1 : 0));
		}
	}

	//------------------

	static uint32_t count_shapes_for(uint32_t count) {
		if (count == 1) {
			return 2;
		}
		else if (count % 2 == 0) {
			return 5 * count;
		}
		else {
			return 4 * count;
		}
	}

	static std::vector< Shape > make_shapes_for(uint32_t count) {
		std::vector< Shape > ret;
		if (count == 1) {
			ret.reserve(2);
			ret.emplace_back(0, BackLeft);
			ret.emplace_back(0, FrontLeft);
		}
		else if (count % 2 == 0) {
			ret.reserve(5 * count);
			for (uint32_t i = 0; i < count; ++i) {
				ret.emplace_back(i, 0);
				ret.emplace_back(i, BackLeft | BackRight);
				ret.emplace_back(i, BackLeft | FrontRight);
				ret.emplace_back(i, FrontLeft | BackRight);
				ret.emplace_back(i, FrontLeft | FrontRight);
			}
			assert(ret.size() == count * 5);
		}
		else {
			ret.reserve(4 * count);
			for (uint32_t i = 0; i < count; ++i) {
				ret.emplace_back(i, BackLeft);
				ret.emplace_back(i, FrontLeft);
				ret.emplace_back(i, BackRight);
				ret.emplace_back(i, FrontRight);
			}
			assert(ret.size() == count * 4);
		}
		return ret;
	}

	uint32_t index_for(uint32_t count) const {
		if (count % 2 == 0) {
			uint32_t ret = roll * 5;
			if (nibbles == 0) ret += 0;
			else if (nibbles == (BackLeft | BackRight)) ret += 1;
			else if (nibbles == (BackLeft | FrontRight)) ret += 2;
			else if (nibbles == (FrontLeft | BackRight)) ret += 3;
			else if (nibbles == (FrontLeft | FrontRight)) ret += 4;
			else assert(0 && "invalid even nibbles config");
			return ret;
		}
		else {
			uint32_t ret = roll * 4;
			if (nibbles == BackLeft) ret += 0;
			else if (nibbles == FrontLeft) ret += 1;
			else if (nibbles == BackRight) ret += 2;
			else if (nibbles == FrontRight) ret += 3;
			else assert(0 && "invalid odd nibbles config");
			return ret;
		}
	}

	//------------------

	PackedShape pack() const {
		return uint32_t(nibbles) | (roll << 4);
	}
	static Shape unpack(PackedShape ps) {
		return Shape(ps >> 4, ps & 0x0f);
	}
};
