#pragma once
#include <cstdint>
#include <glm/glm.hpp>

struct TracedStitch {
	uint32_t yarn = -1U; //yarn ID (why is this on a yarn_in? I guess the schedule.cpp code will tell me someday.
	//ins and outs are in construction order (OLD was: CW direction):
	uint32_t ins[2] = { -1U, -1U };
	uint32_t outs[2] = { -1U, -1U };
	enum Type : char {
		None = '\0',
		Start = 's',
		End = 'e',
		Tuck = 't',
		Miss = 'm',
		Knit = 'k',
		//I am going to add these because the scheduling re-write cares about them, though I'm not sure if they are a good idea to have in a general sense:
		Increase = 'i',
		Decrease = 'd',
	} type = None;
	enum Dir : char {
		CW = 'c', Clockwise = CW,
		AC = 'a', Anticlockwise = AC, CCW = AC, Counterclocwise = AC,
	} dir = AC;

	//useful for debugging and visualization:
	uint32_t vertex = -1U; //vertex of rowcolgraph where created
	glm::vec3 at = glm::vec3(std::numeric_limits< float >::quiet_NaN());
};
