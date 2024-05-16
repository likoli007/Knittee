#pragma once
#include <cstdint>


/*
*	Helper struct used to represent a link between two stitches.
*	Each link has a source stitch and a destination stitch.
*   Widely used during the linkChains() step of the algorithm.
*/
struct Link {
	uint32_t from_chain, from_stitch;
	uint32_t to_chain, to_stitch;
};
