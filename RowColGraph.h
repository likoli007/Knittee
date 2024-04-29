#pragma once
#include <vector>
#include <cassert>
#include <cstdint>
#include <QOpenGLExtraFunctions>
#include <GL/gl.h>
#include "EmbeddedVertex.h"

struct RowColGraph {
	struct Vertex {
		EmbeddedVertex at;
		uint32_t row_in = -1U;
		uint32_t row_out = -1U;
		uint32_t col_in[2] = { -1U, -1U };
		uint32_t col_out[2] = { -1U, -1U };
		void add_col_in(uint32_t i) {
			if (col_in[0] == -1U) col_in[0] = i;
			else if (col_in[1] == -1U) col_in[1] = i;
			else assert(col_in[0] == -1U || col_in[1] == -1U); //no room!
		}
		void add_col_out(uint32_t i) {
			if (col_out[0] == -1U) col_out[0] = i;
			else if (col_out[1] == -1U) col_out[1] = i;
			else assert(col_out[0] == -1U || col_out[1] == -1U); //no room!
		}
	};
	std::vector< Vertex > vertices;
	void clear() {
		vertices.clear();
	}
};