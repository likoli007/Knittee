#pragma once
#include <QVector3D>
#include <glm/glm.hpp>
#include <QOpenGLExtraFunctions>
#include <GL/gl.h>

struct EmbeddedVertex {
	glm::uvec3 simplex;
	QVector3D weights;
	EmbeddedVertex() = default;
	EmbeddedVertex(glm::uvec3 const& simplex_, QVector3D const& weights_) : simplex(simplex_), weights(weights_) { }

	bool operator==(EmbeddedVertex const& o) const {
		return simplex == o.simplex && weights == o.weights;
	}
	bool operator!=(EmbeddedVertex const& o) const {
		return !(*this == o);
	}

	static EmbeddedVertex on_vertex(GLuint a) {
		return EmbeddedVertex(glm::uvec3(a, -1U, -1U), QVector3D(1.0f, 0.0f, 0.0f));
	}
	static EmbeddedVertex on_edge(GLuint a, GLuint b, float mix) {
		if (a > b) {
			std::swap(a, b);
			mix = 1.0f - mix;
		}
		return EmbeddedVertex(glm::uvec3(a, b, -1U), QVector3D(1.0f - mix, mix, 0.0f));
	}

	static EmbeddedVertex canonicalize(glm::uvec3 s, QVector3D w) {
		if (s.x > s.y) {
			std::swap(s.x, s.y);
			float temp = w.x();
			w.setX(w.y());
			w.setY(temp);
		}
		if (s.y > s.z) {
			std::swap(s.y, s.z);
			float temp = w.y();
			w.setY(w.z());
			w.setZ(temp);
		}
		// why still?
		if (s.x > s.y) {
			std::swap(s.x, s.y);

			float temp = w.x();
			w.setX(w.y());
			w.setY(temp);
		}

		return EmbeddedVertex(s, w);
	}

	static EmbeddedVertex mix(EmbeddedVertex const& a, EmbeddedVertex const& b, float m) {
		glm::uvec3 common = common_simplex(a.simplex, b.simplex);
		QVector3D result = a.weights_on(common) * (1.0f - m) + b.weights_on(common) * m;
		return EmbeddedVertex(
			common,
			result
		);
	}

	template< typename T >
	T interpolate(std::vector< T > const& values) const {
		T ret = values[simplex.x] * weights.x();
		if (simplex.y != -1U) ret += values[simplex.y] * weights.y();
		if (simplex.z != -1U) ret += values[simplex.z] * weights.z();
		return ret;
	}

	QVector3D weights_on(glm::uvec3 simplex2) const {
		QVector3D ret(0, 0, 0);
		uint32_t o = 0;
		for (uint32_t i = 0; i < 3; ++i) {
			if (simplex[i] == -1U) break;
			while (simplex2[o] < simplex[i]) {
				++o;
				//assert(o < 3);
			}
			//assert(simplex2[o] == simplex[i]);
			ret[o] = weights[i];
		}
		return ret;
	}

	static glm::uvec3 common_simplex(const glm::uvec3& a, const glm::uvec3& b) {
		glm::ivec3 ret;
		uint32_t ia = 0;
		uint32_t ib = 0;
		for (uint32_t o = 0; o < 3; ++o) {
			if (a[ia] == b[ib]) {
				ret[o] = a[ia];
				++ia; ++ib;
			}
			else if (a[ia] < b[ib]) {
				ret[o] = a[ia];
				++ia;
			}
			else {
				assert(a[ia] > b[ib]);
				ret[o] = b[ib];
				++ib;
			}
		}
		//assert(ia == 3 || a[ia] == -1U);
		//assert(ib == 3 || b[ib] == -1U);
		return ret;
	}
};