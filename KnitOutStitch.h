#pragma once
#include <vector>
#include <string>
#include <glm/glm.hpp>
#include <cstdint>
#include <QString>
#include <QFile>
#include <QDataStream>
#include <QDebug>
/*
*	Originally called 'Stitch' in autoknit, however since there are many 'Stitch' classes in the project,
		(including one directly named 'Stitch'), renamed it to 'KnitOutStitch'
*/

struct KnitOutStitch
{
	//which yarn the stitch is being made with:
	uint32_t yarn = 0;
	//type of stitch (determines how many of in/out are used):
	enum : char {
		//0-in, 1-out:
		Start = 's',
		//1-in, 0-out:
		End = 'e',
		//1-in, 1-out:
		Tuck = 't',
		Miss = 'm',
		Knit = 'k',
		//TODO: Purl = 'p',
		//1-in, 2-out:
		Increase = 'i',
		//2-in, 1-out:
		Decrease = 'd',
	};
	char type = Knit;
	//direction of stitch (relative to current tube):
	enum : char {
		CW = 'c', Clockwise = CW,
		AC = 'a', Anticlockwise = AC, CCW = AC, Counterclockwise = AC,
	};
	char direction = CW;
	//ins and outs are in construction order:
	uint32_t in[2] = { -1U, -1U };
	uint32_t out[2] = { -1U, -1U };
	//DEBUG info:
	glm::vec3 at;

	//helpers:
	uint32_t find_in(uint32_t s) const {
		if (in[0] == s) return 0;
		else if (in[1] == s) return 1;
		else return -1U;
	}
	uint32_t find_out(uint32_t s) const {
		if (out[0] == s) return 0;
		else if (out[1] == s) return 1;
		else return -1U;
	}
	bool check_type() const {
		if (type == Start) {
			return (in[0] == -1U && in[1] == -1U && out[0] != -1U && out[1] == -1U);
		}
		else if (type == End) {
			return (in[0] != -1U && in[1] == -1U && out[0] == -1U && out[1] == -1U);
		}
		else if (type == Tuck || type == Miss || type == Knit) {
			return (in[0] != -1U && in[1] == -1U && out[0] != -1U && out[1] == -1U);
		}
		else if (type == Increase) {
			return (in[0] != -1U && in[1] == -1U && out[0] != -1U && out[1] != -1U);
		}
		else if (type == Decrease) {
			return (in[0] != -1U && in[1] != -1U && out[0] != -1U && out[1] == -1U);
		}
		else {
			return false;
		}
	}
	bool check_direction() const {
		return (direction == Clockwise || direction == Anticlockwise);
	}

};

inline QString to_string(KnitOutStitch const& s) {
	QString ret;
	ret += "(";
	ret += std::to_string(int32_t(s.in[0]));
	ret += ",";
	ret += std::to_string(int32_t(s.in[1]));
	ret += ")-";
	ret += s.type;
	ret += s.direction;
	ret += "->(";
	ret += std::to_string(int32_t(s.out[0]));
	ret += ",";
	ret += std::to_string(int32_t(s.out[1]));
	ret += ")";
	return ret;
};

bool load_stitches(QString const& filename, std::vector< KnitOutStitch >* into);
//void save_stitches(std::string const& filename, std::vector< KnitOutStitch > const& from);


