#include "LaceKnitter.h"




LaceKnitter::LaceKnitter(QObject* parent) : QObject(parent) {
	qDebug() << "loaded LaceKnitter";
}

void LaceKnitter::createSheet(int w, int h, int racking)
{
	sheet.clear();

	sheetWidth = w;
	sheetHeight = h;
	maxOffset = racking;
	for (int i = 0; i < h; i++)
	{
		std::vector<FlatPoint> row;
		for (int j = 0; j < w; j++)
		{
			FlatPoint fp;
			row.push_back(fp);
		}
		sheet.push_back(row);
	}
}

void LaceKnitter::setDimensions(int w, int h, int racking)
{
	sheetWidth = w;
	sheetHeight = h;
	maxOffset = racking;
}

void LaceKnitter::moveLoop(QPair<int, int> from, QPair<int, int> to)
{

	if (abs(to.second - from.second) > maxOffset) {
		qDebug() << "cannot rack to this value!";
		return;
	}

	sheet[from.first][from.second].offset = to.second - from.second;

	emit sheetChanged(&sheet);
}

void LaceKnitter::rackingChanged(int r) {
	maxOffset = r;
}

void LaceKnitter::widthChanged(int increase, int left) {


	if (sheetWidth + increase < 2) {
		qDebug() << "cannot have a width of 0!";
		return;
	}

	sheetWidth += increase;

	if (increase == 1) {
		if (left == 1) {
			for (int i = 0; i < sheet.size(); i++) {
				FlatPoint fp;
				sheet[i].insert(sheet[i].begin(), fp);
			}
		}
		else {
			for (int i = 0; i < sheet.size(); i++) {
				FlatPoint fp;
				sheet[i].push_back(fp);
			}
		}
	}
	else {
		if (left == 1) {
			for (int i = 0; i < sheet.size(); i++) {
				sheet[i].erase(sheet[i].begin());
			}
		}
		else {
			for (int i = 0; i < sheet.size(); i++) {
				sheet[i].pop_back();
			}
		}
	}
	emit(sheetChanged(&sheet));
}

void LaceKnitter::heightChanged(int increase, int top) {

	//qDebug() << "arrrivaaa" << increase;
	if (sheetHeight + increase < 2) {
		qDebug() << "cannot have a height of 0!";
		return;
	}

	sheetHeight += increase;

	if (increase == 1) {
		if (top == 1) {
			std::vector<FlatPoint> row;
			for (int i = 0; i < sheetWidth; i++) {
				FlatPoint fp;
				row.push_back(fp);
			}
			sheet.insert(sheet.begin(), row);
		}
		else {
			std::vector<FlatPoint> row;
			for (int i = 0; i < sheetWidth; i++) {
				FlatPoint fp;
				row.push_back(fp);
			}
			sheet.push_back(row);
		}
	}
	else {
		if (top == 1) {
			sheet.erase(sheet.begin());
		}
		else {
			sheet.pop_back();
		}
	}
	emit sheetChanged(&sheet);
}

void LaceKnitter::loadFromFile(QString filename)
{
	sheet.clear();

	qDebug() << "loading from file";

	filename += "/sheet.dat";
	QFile file(filename);
	if (!file.open(QIODevice::ReadOnly))
	{
		qDebug() << "could not open file";
		return;
	}

	QDataStream in(&file);

	while (!in.atEnd())
	{

		for (int i = 0; i < sheetHeight; i++)
		{
			std::vector<FlatPoint> row;
			for (int j = 0; j < sheetWidth; j++) {
				FlatPoint fp;
				in >> fp;
				row.push_back(fp);
			}
			sheet.push_back(row);
		}

	}

	file.close();
	qDebug() << "loaded from file";
	emit sheetChanged(&sheet);
}

int LaceKnitter::getWidth()
{
	return sheetWidth;
}
int LaceKnitter::getHeight()
{
	return sheetHeight;
}
int LaceKnitter::getRacking()
{
	return maxOffset;
}

void LaceKnitter::saveToFile(QString filename)
{
	qDebug() << "saving to file";

	filename += "/sheet.dat";
	QFile file(filename);

	if (!file.open(QIODevice::WriteOnly))
	{
		qDebug() << "could not open file";
		return;
	}

	QDataStream out(&file);

	for (int i = 0; i < sheet.size(); i++)
	{
		for (int j = 0; j < sheet[i].size(); j++)
		{
			out << sheet[i][j];
		}
	}

	file.close();
	//qDebug() << "saved to file";
}

void LaceKnitter::generateKnitout(int algo) {
	k.clear();
	QString knitout = ";;Carriage settings\n";
	chosenAlgorithm = algo;
	//iterate through sheet backwards


	QString dir = "-";
	QString carrier = "5";


	for (int i = 0; i <= sheet.size() - 1; i++) {
		offsets.clear();
		for (int j = 0; j < sheet[i].size(); j++) {
			offsets.push_back(sheet[i][j].offset);
		}
		qDebug() << offsets;

		if (dir == "-") {
			for (int j = sheet[i].size(); j >= 1; j--) {
				k.append("knit - f" + QString::number(j) + " " + carrier);
			}
			dir = "+";
		}
		else {
			for (int j = 1; j <= sheet[i].size(); j++) {
				k.append("knit + f" + QString::number(j) + " " + carrier);
			}
			dir = "-";
		}

		k.append("; xfer pass [line " + QString::number(i) + "]");
		if (chosenAlgorithm == 1) {
			schoolbus(offsets, -maxOffset, maxOffset);
		}
		else if (chosenAlgorithm == 0) {
			cse(maxOffset);
		}
		else {
			qDebug() << "algorithm not implemented";
		}

	}
	k.push_back("outhook " + carrier);
	for (QString line : k) {
		qDebug() << line;
	}

}

void LaceKnitter::xfer(QString& fromBed, int fromIndex, QString& toBed, int toIndex) {

	int fromNeedle = fromIndex + 1;
	int toNeedle = toIndex + 1;
	int r = (fromBed[0] == 'b' ? toIndex - fromIndex : fromIndex - toIndex);
	if (racking != r) {
		k.push_back("rack " + QString::number(r));
		racking = r;
	}
	k.push_back("xfer " + fromBed + QString::number(fromNeedle) + " " + toBed + QString::number(toNeedle));
	if (racking != 0) {
		k.push_back("rack 0");
		racking = 0;
	}
}


void LaceKnitter::schoolbus(const std::vector<int>& offsets, int minRacking, int maxRacking) {

	QString f = "f";
	QString b = "b";

	for (size_t i = 0; i < offsets.size(); i++) {
		xfer(f, i, b, i);
	}

	int currentOffset = minRacking;
	while (currentOffset <= maxRacking) {
		for (size_t i = 0; i < offsets.size(); i++) {
			if (offsets[i] == currentOffset) {
				xfer(b, i, f, i + currentOffset);
			}
		}
		currentOffset++;
	}


}


constexpr int top = 0;

constexpr int parent(int i) {
	return ((i + 1) >> 1) - 1;
}

constexpr int left(int i) {
	return (i << 1) + 1;
}

constexpr int right(int i) {
	return (i + 1) << 1;
}

int penalty(const State& state) {
	int p = 0;
	for (int i = 0; i < state.offsets.size(); i++) {
		p += abs(state.offsets[i]);
	}
	return p;
}


struct CompareState {
	bool operator()(const State& s1, const State& s2) {
		if (penalty(s1) == penalty(s2)) {
			if (s1.path.size() == s2.path.size()) {
				return s1.r - s1.l > s2.r - s2.l;
			}
			return s1.path.size() < s2.path.size();
		}
		return penalty(s1) < penalty(s2);
	}
};


class PriorityQueue {
private:
	std::vector<State> heap;
	CompareState comparator;

	int parent(int i) {
		return ((i + 1) >> 1) - 1;
	}

	int left(int i) {
		return (i << 1) + 1;
	}

	int right(int i) {
		return (i + 1) << 1;
	}

	void swap(int i, int j) {
		State temp = heap[i];
		heap[i] = heap[j];
		heap[j] = temp;
	}

	void siftUp(int node) {
		while (node > 0 && comparator(heap[node], heap[parent(node)])) {
			swap(node, parent(node));
			node = parent(node);
		}
	}

	void siftDown(int node) {
		while ((left(node) < heap.size() && comparator(heap[left(node)], heap[node])) ||
			(right(node) < heap.size() && comparator(heap[right(node)], heap[node]))) {
			int maxChild = (right(node) < heap.size() && comparator(heap[right(node)], heap[left(node)])) ? right(node) : left(node);
			swap(node, maxChild);
			node = maxChild;
		}
	}

public:
	int size() {
		return heap.size();
	}

	bool isEmpty() {
		return size() == 0;
	}

	State peek() {
		return heap[0];
	}

	void push(State value) {
		heap.push_back(value);
		siftUp(size() - 1);
	}

	State pop() {
		State poppedValue = peek();
		if (size() > 1) {
			swap(0, size() - 1);
		}
		heap.pop_back();
		siftDown(0);
		return poppedValue;
	}
};


State new_state(State state = State()) {
	State s;
	if (state.current.size() != 0) {
		s.current = state.current;
		s.offsets = state.offsets;
		s.path = state.path;
		s.l = state.l;
		s.r = state.r;
		s.prev = state.prev;
		s.Do = state.Do;
		s.rack = state.rack;
		s.chain = state.chain;
		s.chain.push_back(state);
	}
	return s;
}

bool LaceKnitter::state_respects_slack(const State& state) {
	std::vector<char> beds;
	for (int i = 0; i < n_stitches; i++) {
		if (i < state.l) {
			beds.push_back('f');
		}
		else if (i >= state.l && i <= state.r) {
			beds.push_back('b');
		}
		else if (i > state.r) {
			beds.push_back('f');
		}
		else {
			assert(false && "Invalid bed state");
		}
	}

	for (int i = 1; i < n_stitches; i++) {
		int slack = std::max(1, std::abs((i + offsets[i]) - (i - 1 + offsets[i - 1])));
		int sep = std::abs(state.current[i] - state.current[i - 1]);
		if (beds[i] == beds[i - 1] && sep > slack) {
			//std::cout << "stretching too much on the same bed" << std::endl;
			return false;
		}
		if (beds[i] != beds[i - 1]) {
			int back = (beds[i] == 'b' ? state.current[i] : state.current[i - 1]);
			int front = (beds[i] == 'f' ? state.current[i] : state.current[i - 1]);
			int stretch = std::abs(back + state.rack - front);
			if (stretch > slack) {
				//std::cout << "stretching too much " << stretch << " between beds at current rack " << state.rack << " slack " << slack << std::endl;
				return false;
			}
		}
	}

	return true;
}



void LaceKnitter::generate_transfers(std::vector<LaceTransfer>& path, bool log) {
	// keeping track of multiple transfers, similar to driver
	// std::cout << path << std::endl;
	std::map<QString, std::vector<int>> innerneedles;
	qDebug() << "generating";
	n_stitches = offsets.size(); // Assuming n_stitches is known or passed as an argument
	log = false;
	for (int i = 0; i < n_stitches; ++i) {
		innerneedles["f" + QString::number(i)] = { i };
	}

	for (auto& e : path) {
		if (e.from_bed != 'f' && e.from_bed != 'b') {
			assert(false && "Invalid from_bed value");
		}
		if (e.to_bed != 'f' && e.to_bed != 'b') {
			assert(false && "Invalid to_bed value");
		}

		if (false) {
			qDebug() << "\txfer " << e.from_bed << e.from << " " << e.to_bed << e.to;
		}
		else {
			qDebug() << "\txfer " << e.from_bed << e.from << " " << e.to_bed << e.to;
			auto& from = innerneedles[e.from_bed + QString::number(e.from)];
			if (from.empty()) continue; // empty needle
			xfer(e.from_bed, e.from, e.to_bed, e.to); // Assuming xfer function exists
			if (needles.find(e.to_bed + QString::number(e.to)) == needles.end()) {
				QString line = e.to_bed + QString::number(e.to);
				needles.value(line) = {};
			}

			QString line = e.to_bed + QString::number(e.to);

			QVector<int> to = needles.value(line);

			while (!from.empty()) { to.push_back(from.back()); from.pop_back(); }
			needles.insert(line, to);
		}
	}
}





QString stringify(LaceTransfer& lace) {
	QJsonObject jsonObject;

	jsonObject["from_bed"] = lace.from_bed;
	jsonObject["from"] = lace.from;
	jsonObject["to_bed"] = lace.to_bed;
	jsonObject["to"] = lace.to;

	QJsonDocument jsonDoc(jsonObject);

	// Convert JSON document to a string
	QString jsonString = jsonDoc.toJson(QJsonDocument::Compact);

	return jsonString;
}

QString stringify(State& state) {
	QJsonObject jsonObject;

	QJsonArray jsonArray;
	for (int value : state.current) {
		jsonArray.append(value);
	}

	jsonObject["current"] = jsonArray;

	QJsonArray jsonArray2;
	for (int value : state.offsets) {
		jsonArray2.append(value);
	}

	jsonObject["offsets"] = jsonArray2;

	jsonObject["l"] = state.l;
	jsonObject["r"] = state.r;
	jsonObject["rack"] = state.rack;
	jsonObject["Do"] = state.Do;

	QJsonArray jsonArray3;
	for (int value : state.prev) {
		jsonArray3.append(value);
	}

	jsonObject["prev"] = jsonArray3;

	QJsonArray jsonArray4;
	for (LaceTransfer value : state.path) {
		jsonArray4.append(stringify(value));
	}

	jsonObject["path"] = jsonArray4;

	if (state.chain.size() == 0) {
		jsonObject["chain"] = QJsonValue::Null;
	}
	else {
		QJsonArray jsonArray5;
		for (State value : state.chain) {
			jsonArray5.append(stringify(value));
		}
		jsonObject["chain"] = jsonArray5;
	}

	QJsonDocument jsonDoc(jsonObject);

	// Convert JSON document to a string
	QString jsonString = jsonDoc.toJson(QJsonDocument::Compact);

	return jsonString;

}

void LaceKnitter::visit_state(State state, std::set<QString>& visited) {
	State s = new_state(state);
	s.path.clear();
	s.chain.clear();
	//s.chain.push_back(state);
	//qDebug() << s.path.size() << stringify(s) ;
	visited.insert(stringify(s)); // Assuming JSON.stringify(s) converts the state to a string
}

bool LaceKnitter::has_state(State state, std::set<QString>& visited) {
	State s = new_state(state);
	s.path.clear();
	s.chain.clear();

	return visited.count(stringify(s)) > 0;
}

std::vector<State> LaceKnitter::priority_order(std::vector<State> states) {
	std::sort(states.begin(), states.end(), [this](State& a, State& b) {
		int keyA = penalty(a);
		int keyB = penalty(b);

		if (keyA < keyB) return true;
		if (keyA > keyB) return false;

		// If penalty is identical, pick the one with fewer transfers in its path
		if (a.path.size() < b.path.size()) return true;
		if (a.path.size() > b.path.size()) return false;

		// Break ties by something else
		if (a.l - a.r > b.l - b.r) return true;
		if (a.l - a.r < b.l - b.r) return false;

		return false; // Default case
		});
	return states;
}

bool LaceKnitter::okay_to_move_index_by_offset(const State& state, int idx, int ofs) {
	assert(idx >= 0 && idx < n_stitches && "idx is not valid");
	assert(state.r >= 0 && state.r < n_stitches && "r is not valid");
	assert(state.l >= 0 && state.l < n_stitches && "l is not valid");

	// not causing slack problems
	std::vector<char> beds(n_stitches);
	for (int i = 0; i < n_stitches; i++) {
		if (i < state.l) {
			beds[i] = 'f';
		}
		else if (i >= state.l && i <= state.r) {
			beds[i] = 'b';
		}
		else if (i > state.r) {
			beds[i] = 'f';
		}
		else {
			assert(false && "?");
		}
	}
	beds[idx] = 'f'; // idx will move to the front

	if (idx > 0 && beds[idx] != beds[idx - 1]) {
		int back = state.current[idx - 1];
		int front = state.current[idx] + ofs;
		int stretch = std::abs(back + ofs - front);
		int slack = std::max(1, std::abs(idx + offsets[idx] - (idx - 1 + offsets[idx - 1])));
		if (stretch > slack) return false;
	}
	if (idx + 1 < n_stitches && beds[idx] != beds[idx + 1]) {
		int back = state.current[idx + 1];
		int front = state.current[idx] + ofs;
		int stretch = std::abs(back + ofs - front);
		int slack = std::max(1, std::abs(idx + offsets[idx] - (idx + 1 + offsets[idx + 1])));
		if (stretch > slack) return false;
	}

	// not stacking over a stitch that wants to go elsewhere, not tangling
	int c = state.current[idx] + ofs;
	int o = state.offsets[idx] - ofs;
	for (int i = 0; i < state.l; i++) {
		if (state.current[i] == c && (state.offsets[i]) != (o)) {
			assert(i != idx && "shouldn't happen, right?");
			return false;
		}
		if (state.current[i] > c) {
			return false;
		}
	}
	for (int i = state.r + 1; i < n_stitches; i++) {
		if (state.current[i] == c && (state.offsets[i]) != (o)) {
			assert(i != idx && "shouldn't happen, right?");
			return false;
		}
		if (state.current[i] < c) {
			return false;
		}
	}

	return true;
}

void LaceKnitter::cse(int maxRacking) {

	bool verbose = false;
	bool strict_order = false; // !options[0];

	const int Expanding = 0;
	const int StretchToBack = 1;
	const int ExpandToFront = 2;
	const std::vector<QString> action = { "Expanding......", "Stretch-to-back", "Expand-to-front" };

	n_stitches = offsets.size();

	qDebug() << maxRacking;
	std::vector<int> slack_forward(n_stitches, maxRacking);
	std::vector<int> slack_backward(n_stitches, maxRacking);
	std::vector<int> current(n_stitches);
	std::set<QString> visited;

	auto accum_add = [](int accumulator, int currentValue) { return abs(accumulator) + abs(currentValue); };


	std::vector<int> target;
	for (int i = 0; i < n_stitches; i++) {
		current.push_back(i);
		target.push_back(i + offsets[i]);
	}

	for (int i = 0; i < n_stitches - 1; i++) {
		slack_forward[i] = std::max(1, std::abs(target[i] - target[i + 1]));
	}
	for (int i = 1; i < n_stitches; i++) {
		slack_backward[i] = std::max(1, std::abs(target[i - 1] - target[i]));
	}




	auto comparator = [this](const State& a, const State& b) {
		int penaltyA = penalty(a);
		int penaltyB = penalty(b);

		if (penaltyA == penaltyB) {
			if (a.path.size() == b.path.size()) {
				return a.r - a.l > b.r - b.l;
			}
			else {
				return a.path.size() < b.path.size();
			}
		}
		else {
			return penaltyA < penaltyB;
		}
	};

	PriorityQueue PQ;


	std::vector<LaceTransfer> path;
	std::vector<State> chain;
	State first_state = { current, current, offsets, StretchToBack, path, 0, n_stitches - 1, chain ,  0 };


	PQ.push(first_state);

	qDebug() << "original:" << current << offsets;

	int last_penalty = INT_MAX;

	while (!PQ.isEmpty()) {
		State s = PQ.pop();


		qDebug() << PQ.size() << s.Do << s.offsets;
		//qDebug() << s.Do;

		visit_state(s, visited);


		if (s.Do == StretchToBack) {
			assert(penalty(s) < last_penalty);
			last_penalty = penalty(s);
		}


		//assert(state_respects_slack(s));

		if (penalty(s) == 0) {

			if (s.Do == Expanding) {
				for (int i = s.l; i <= s.r; i++) {
					LaceTransfer t = { "b", s.current[i], "f", s.current[i] };
					s.path.push_back(t);
				}
			}
			generate_transfers(s.path);

			return;

		}

		if (s.Do == StretchToBack) {

			//qDebug() << "stretch";
			//console.log('\t\tStretch state to back');
			// the easy case, only translate.
			// as long as within min-max range of needles, always safe
			// TODO add a min-max needle range
			for (int r = -maxRacking; r <= maxRacking; r++) {
				State ss = new_state(s);
				for (int i = 0; i < n_stitches; i++) {
					//assert(ss.current[i] !=  undefined);
					//assert(ss.current[i] - r !=  undefined);
					ss.path.push_back({ "f", ss.current[i], "b", ss.current[i] - r });

					ss.current[i] -= r;
					ss.offsets[i] += r;
				}
				ss.Do = Expanding;
				ss.prev = ss.current;
				ss.l = 0;
				ss.r = n_stitches - 1;
				ss.rack = 0; // at the end of this operation, rack can be reset
				if (state_respects_slack(ss) /*&& penalty(ss) <= last_penalty*/ && !has_state(ss, visited)) {
					//States.push(ss);
					qDebug() << "stretch respects";
					PQ.push(ss);
					//print out state info
					//qDebug() << ss->current << ss->offsets << ss->l << ss->r << ss->rack;
				}
			}
		}//collapse
		else if (s.Do == Expanding) {
			//qDebug() << "expand";
			assert(s.l <= s.r);
			int min_l = -maxRacking;
			int min_r = -maxRacking;
			int max_l = maxRacking;
			int max_r = maxRacking;
			// Find the next legal state(s) to move s.l 
			{

				int n = s.current[s.l];
				//console.log('\t\t ** Moving l in range ', min_l, max_l);
				// go through each offset
				for (auto o = min_l; o <= max_l; o++) {

					//console.log('attempting to move by', o);
					{

						if (true) {
							auto next = new_state(s);

							next.rack = o;
							next.current[next.l] += o;
							next.offsets[next.l] -= o;
							qDebug() << "\tnext" << next.offsets;
						}

						qDebug() << "checking okay to move for" << s.l << o;
						if (okay_to_move_index_by_offset(s, s.l, o)) {
							// launch a new version 
							auto next = new_state(s);

							next.rack = o;
							//qDebug() << "\tnext" << next.offsets;
							if (!state_respects_slack(next)) continue;

							next.path.push_back({ "b", next.current[next.l], "f",next.current[next.l] + o });
							next.current[next.l] += o;
							next.offsets[next.l] -= o;
							//	console.log('(L)Found offset ' + o.toString() + ' that works.');
							//	console.log('\t l:', next.l,'current', next.current[next.l],'offset', next.offsets[next.l]);
							qDebug() << "\t\trespects";


							//qDebug() << next.path.size();
							if (next.l < next.r) {
								next.l++;
								next.Do = Expanding;

								if (!has_state(next, visited) /*&& penalty(next) <= last_penalty*/ && state_respects_slack(next)) {
									PQ.push(next);
									//qDebug() << "aoru";
									//States.push(next);
									//console.log('\t adding next l ('+next.l.toString()+')');
								}
							}
							else if (next.l == next.r) {
								next.Do = ExpandToFront;
								next.l = 0;
								next.r = n_stitches - 1;
								if (!has_state(next, visited) /*&& penalty(next) <= last_penalty*/ && state_respects_slack(next)) {
									PQ.push(next);
									//qDebug() << "mauru";
									//States.push(next);
									//console.log('\t adding next(l) for expansion');
								}

							}
						}


					}
				}
			}

		}


		else if (s.Do == ExpandToFront) { // if collapsed state expand to front
			//qDebug() << "expand to front";
			//console.log('\t\tExpand state to front', s.prev, s.current, s.offsets);
			auto ss = new_state(s);

			// expanding has already done all the work, so this is trivial 
			// assuming prev was not clobbered by anything (it shouldn't be)
			// but also need to do this in order

			ss.Do = StretchToBack;
			// possibly figure out if fewer passes are possible 
			ss.rack = 0; // at the end of this operation rack can be reset
			//assert(state_respects_slack(ss), 'this should be slack friendly, right?');

			//qDebug() << penalty(ss) << last_penalty;
			if (penalty(ss) < last_penalty && !has_state(ss, visited)) {
				qDebug() << "wowee";
				PQ.push(ss);
				//States.push(ss);
			}

		}


	}

	qDebug() << "ayayayayayay";

}