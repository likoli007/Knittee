#pragma once


/*
*	Helper class representing a stitch on a KnitGraph (RowColGraph), used mainly to implement valid linking of stitches 
*/

struct Stitch {
	float t; //position along chain [0,1)
	unsigned int vertex; //used to track stitches when building RowColGraph
	enum Flag : char {
		FlagDiscard = -1,
		FlagLinkOne = 1,
		FlagLinkAny = 2,
	};
	Flag flag; //what sort of links are allowed to this stitch (or if this stitch is just marked for removal)
	Stitch(float t_, Flag flag_, unsigned int vertex_ = -1U) : t(t_), vertex(vertex_), flag(flag_) { }
};

