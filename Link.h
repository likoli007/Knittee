#pragma once
#include <cstdint>

struct Link {
	uint32_t from_chain, from_stitch;
	uint32_t to_chain, to_stitch;
};
