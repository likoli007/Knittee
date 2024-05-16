#include "KnitOutStitch.h"



bool load_stitches(QString const& filename, std::vector< KnitOutStitch >* into_) {
	assert(into_);
	auto& into = *into_;
	into.clear();

	QFile file(filename);

	if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
		qDebug() << "Failed to open file.";
		return false;
	}

	//std::string line;
	QTextStream in(&file);
	while (!in.atEnd()) {
		QString line = in.readLine();
		QStringList parts = line.split(" ");
		int32_t in[2];
		int32_t out[2];

		if (parts.size() != 10) {
			qDebug() << "ERROR: Failed to read stitch.";
			return false;
		}

		KnitOutStitch temp;

		temp.yarn = parts[0].toInt();
		temp.type = parts[1].at(0).toLatin1();
		temp.direction = parts[2].at(0).toLatin1();
		in[0] = parts[3].toInt();
		in[1] = parts[4].toInt();
		out[0] = parts[5].toInt();
		out[1] = parts[6].toInt();
		temp.at.x = parts[7].toFloat();
		temp.at.y = parts[8].toFloat();
		temp.at.z = parts[9].toFloat();

		temp.in[0] = in[0];
		temp.in[1] = in[1];
		temp.out[0] = out[0];
		temp.out[1] = out[1];

		qDebug() << "Read stitch: " << temp.yarn << " " << temp.type << " " << temp.direction << " " << temp.in[0] << " " << temp.in[1] << " " << temp.out[0] << " " << temp.out[1] << " " << temp.at.x << " " << temp.at.y << " " << temp.at.z;

		if (!temp.check_type()) {
			qDebug() << "ERROR: Stitch does not have proper in/out for type.";
			qDebug() << "  line: '" << line << "'";
			return false;
		}


		into.emplace_back(temp);
	}
	for (auto const& s : into) {
		uint32_t idx = &s - &into[0];
		auto check_in = [&](uint32_t in_idx) -> bool {
			if (in_idx == -1U) {
				return true;
			}
			else {
				if (in_idx >= idx) return false;
				if (into[in_idx].find_out(idx) == -1U) return false;
				return true;
			}
		};
		if (!check_in(s.in[0]) || !check_in(s.in[1])) {
			qDebug() << "Stitch does not have proper 'in' array.";
			return false;
		}
		auto check_out = [&](uint32_t out_idx) -> bool {
			if (out_idx == -1U) {
				return true;
			}
			else {
				if (out_idx >= into.size()) return false;
				if (out_idx <= idx) return false;
				if (into[out_idx].find_in(idx) == -1U) return false;
				return true;
			}
		};
		if (!check_out(s.out[0]) || !check_out(s.out[1])) {
			qDebug() << "Stitch does not have proper 'out' array.";
			return false;
		}
	}
	return true;
}
