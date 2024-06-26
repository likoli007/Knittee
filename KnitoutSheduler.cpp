﻿#include "KnitoutSheduler.h"



KnitoutScheduler::KnitoutScheduler(QObject* parent) : QObject(parent)
{
	maxRacking = 3;
}

void KnitoutScheduler::setFilePath(QString path)
{
	filePath = path;
}

void KnitoutScheduler::schedule()
{
	QString head1 = QString(";!knitout-2");
	QString head2 = QString(";;Carriers: 1 2 3 4 5 6 7 8 9 10");


	out(head1);
	out(head2);



	std::vector< KnitOutStitch > stitches;
	std::vector< std::pair<BedNeedle, BedNeedle>> stitch_locations, input_locations;
	if (!load_stitches(filePath, &stitches)) {
		qDebug() << "ERROR: failed to load stitches";
		return;
	}
	emit helpBoxCommunication(  "Read " +QString::number(stitches.size()) + " stitches");
	stitch_locations.assign(stitches.size(), std::make_pair(BedNeedle(BedNeedle::FrontSliders, -1), BedNeedle(BedNeedle::FrontSliders, -1)));
	input_locations.assign(stitches.size(), std::make_pair(BedNeedle(BedNeedle::FrontSliders, -1), BedNeedle(BedNeedle::FrontSliders, -1)));

	typedef uint32_t StorageIdx;
	struct Storage : std::deque< Loop > {
		//Shape shape;
	};

	struct ScheduleOption {
		ScheduleCost cost;
		std::vector< uint32_t > in_order; //indexes into step.in
		std::vector< PackedShape > in_shapes; //in_shapes[i] is shape for step.in[i]
		PackedShape inter_shape;
		std::vector< uint32_t > out_order; //indexes into step.out
		std::vector< PackedShape > out_shapes; //out_shapes[i] is shape for step.out[i]
	};

	struct Step {
		uint32_t begin = 0, end = 0; //stitch range constructed in this step
		std::vector< StorageIdx > in, out;
		Storage inter; //configuration of loops after input arrangement but before performing stitches
		std::vector< uint32_t > inter_to_out; //how stitches in inter are transformed to get to out

		//Step reinterprets storages[in[*]] as inter + storages[out[1...]]
		// stretches and rolls inter,
		// and performs knitting on inter to produce storages[out[0]]

		std::vector< ScheduleOption > options;
	};
	
	std::vector< Storage > storages;
	std::vector< Step > steps;

	storages.reserve(2 * stitches.size()); //just how much storage can there actually be?
	steps.reserve(stitches.size());

	//get the out loop for a given 'in':
	auto source_loop = [&stitches](uint32_t si, uint32_t ii) {
		assert(stitches[si].in[ii] < stitches.size());
		Loop ret(stitches[si].in[ii], 0);
		if (stitches[ret.stitch].out[ret.idx] != si) ret.idx = 1;
		assert(stitches[ret.stitch].out[ret.idx] == si);
		return ret;
	};

	{ //build steps:
		//current yarn position w.r.t. active loops:
		struct YarnInfo {
			Loop loop = INVALID_LOOP;
			char direction = '\0';
		};
		std::map< uint32_t, YarnInfo > active_yarns;

		std::map< Loop, StorageIdx > active_loops;

		struct Bridge {
			Loop a, b;
			uint32_t distance = 0;
		};
		std::vector< Bridge > bridges;

		uint32_t next_begin = 0;
		while (next_begin < stitches.size()) {
			//std::cout << "steps[" << steps.size() << "]:\n"; //DEBUG

			steps.emplace_back();
			Step& step = steps.back();

			{ //figure out begin/end range for step:
				step.begin = next_begin;
				step.end = step.begin + 1;
				//steps take as many stitches as they can, given some constraints:
				uint32_t increases = 0, decreases = 0;
				uint32_t last_shaping = step.begin;
				if (stitches[step.begin].type == KnitOutStitch::Increase) {
					++increases;
				}
				if (stitches[step.begin].type == KnitOutStitch::Decrease) {
					++decreases;
				}
				auto basic_type = [](KnitOutStitch const& s) {
					if (s.type == KnitOutStitch::Start) return '+';
					else if (s.type == KnitOutStitch::End) return '-';
					else return '|';
				};
				while (step.end < stitches.size()) {
					//same basic type:
					if (basic_type(stitches[step.end]) != basic_type(stitches[step.begin])) break;
					//same direction:
					if (stitches[step.end].direction != stitches[step.begin].direction) break;
					//same yarn:
					if (stitches[step.end].yarn != stitches[step.begin].yarn) break;
					//same course: (not sure this is really needed, but it does mean that every loop ends up in a Storage at some point)
					if (stitches[step.end].in[0] != -1U && stitches[step.end].in[0] >= step.begin) break;
					if (stitches[step.end].in[1] != -1U && stitches[step.end].in[1] >= step.begin) break;
					//not too many decreases/increases in a segment (...and, if there might be, cut things back a bit)
					if (stitches[step.end].type == KnitOutStitch::Increase) {
						++increases;
						if (increases >= 4) {
							step.end = (step.end - last_shaping + 1) / 2 + last_shaping;
							break;
						}
						last_shaping = step.end;
					}
					if (stitches[step.end].type == KnitOutStitch::Decrease) {
						++decreases;
						if (decreases >= 4) {
							step.end = (step.end - last_shaping + 1) / 2 + last_shaping;
							break;
						}
						last_shaping = step.end;
					}
					++step.end;
				}
				next_begin = step.end;
			}

			std::vector< Loop > in_chain;
			std::vector< Loop > out_chain;

			//stitches take in_chain -> out_chain:
			for (uint32_t si = step.begin; si != step.end; ++si) {
				for (uint32_t i = 0; i < 2; ++i) {
					if (stitches[si].in[i] == -1U) continue;
					Loop l = source_loop(si, i);
					in_chain.emplace_back(l);
				}
				for (uint32_t o = 0; o < 2; ++o) {
					if (stitches[si].out[o] == -1U) continue;
					out_chain.emplace_back(si, o);
				}
			}

			if (stitches[step.begin].direction == KnitOutStitch::CW) {
				std::reverse(in_chain.begin(), in_chain.end());
				std::reverse(out_chain.begin(), out_chain.end());
			}

			YarnInfo& yarn = active_yarns[stitches[step.begin].yarn];

			{ //fill in step's 'in' array based on storages holding in_chain:
				std::set< StorageIdx > used;
				for (auto const& l : in_chain) {
					auto f = active_loops.find(l);
					assert(f != active_loops.end() && "stitches should depend only on active stuff");
					used.insert(f->second);
				}
				if (yarn.loop != INVALID_LOOP) {
					auto f = active_loops.find(yarn.loop);
					assert(f != active_loops.end() && "active yarn should reference active stuff");
					used.insert(f->second);
				}

				step.in.assign(used.begin(), used.end());
			}

			{ //do some loop surgery to figure out 'out' storages:
				std::list< Storage > outs;
				for (auto storage_idx : step.in) {
					outs.emplace_back(storages[storage_idx]);
				}

				//flip and re-jigger yarn storages:
				auto link_ccw = [&outs](Loop const& a, Loop const& b) {
					qDebug() << "Linking " << a.to_string() << " -> " << b.to_string(); //DEBUG
					std::list< Storage >::iterator sa = outs.end();
					uint32_t la = -1U;
					std::list< Storage >::iterator sb = outs.end();
					uint32_t lb = -1U;

					for (auto oi = outs.begin(); oi != outs.end(); ++oi) {
						for (uint32_t li = 0; li < oi->size(); ++li) {
							if ((*oi)[li] == a) {
								assert(sa == outs.end());
								sa = oi;
								assert(la == -1U);
								la = li;
							}
							if ((*oi)[li] == b) {
								assert(sb == outs.end());
								sb = oi;
								assert(lb == -1U);
								lb = li;
							}
						}
					}
					assert(sa != outs.end() && la < sa->size());
					assert(sb != outs.end() && lb < sb->size());
					if (sa == sb) {
						if ((la + 1) % sa->size() == lb) {
							//already ccw! great.
						}
						else {
							//must split loop.
							Storage non_ab;
							Storage ab;
							if (la < lb) {
								ab.insert(ab.end(), sa->begin(), sa->begin() + la + 1);
								non_ab.insert(non_ab.end(), sa->begin() + la + 1, sa->begin() + lb);
								ab.insert(ab.end(), sa->begin() + lb, sa->end());
							}
							else {
								assert(lb < la);
								non_ab.insert(non_ab.end(), sa->begin(), sa->begin() + lb);
								ab.insert(ab.end(), sa->begin() + lb, sa->begin() + la + 1);
								non_ab.insert(non_ab.end(), sa->begin() + la + 1, sa->end());
							}
							//PARANOIA:
							assert(non_ab.size() + ab.size() == sa->size());
							la = -1U;
							lb = -1U;
							for (uint32_t li = 0; li < ab.size(); ++li) {
								if (ab[li] == a) {
									assert(la == -1U);
									la = li;
								}
								if (ab[li] == b) {
									assert(lb == -1U);
									lb = li;
								}
							}
							assert(la != -1U && lb != -1U);
							assert((la + 1) % ab.size() == lb);
							for (auto const& l : non_ab) {
								assert(l != a && l != b);
							}
							//end PARANOIA
							*sa = std::move(ab);
							outs.emplace_back(non_ab);
						}
					}
					else {
						//must merge loops
						std::rotate(sa->begin(), sa->begin() + (la + 1), sa->end());
						assert(sa->back() == a);
						std::rotate(sb->begin(), sb->begin() + lb, sb->end());
						assert(sb->front() == b);
						//TODO: add_bridge(sa->front(), sa->back())
						//TODO: add_bridge(sb->front(), sb->back())
						sa->insert(sa->end(), sb->begin(), sb->end());
						outs.erase(sb);
					}

				};

				if (yarn.loop != INVALID_LOOP) {
					if (yarn.direction != stitches[step.begin].direction) {
						//reversing direction should happen on the same stitch:
						assert(!in_chain.empty() && yarn.loop == (stitches[step.begin].direction == KnitOutStitch::CCW ? in_chain[0] : in_chain.back()));
					}
					else {
						assert(!in_chain.empty() && "Don't handle linking yarn into starts, though that might be useful.");
						if (yarn.direction == KnitOutStitch::CCW) { //formed ccw, so yarn link is before first elt of chain:
							link_ccw(yarn.loop, in_chain[0]);
						}
						else {
							link_ccw(in_chain.back(), yarn.loop);
						}
					}
				}

				for (uint32_t i = 0; i + 1 < in_chain.size(); ++i) {
					link_ccw(in_chain[i], in_chain[i + 1]);
				}

				//Now 'outs' should have in_chain in one ccw list.


				// --> replace with out_chain.
				if (in_chain.empty()) {
					//empty in chain -> create a new cycle~

					step.inter = Storage(); //inter is an empty cycle

					Storage result;
					result.insert(result.end(), out_chain.begin(), out_chain.end());
					outs.emplace_front(result);
				}
				else {
					//check that in_chain is indeed in order:
					std::list< Storage >::iterator s = outs.end();
					uint32_t l = -1U;

					for (auto oi = outs.begin(); oi != outs.end(); ++oi) {
						for (uint32_t li = 0; li < oi->size(); ++li) {
							if ((*oi)[li] == in_chain[0]) {
								assert(s == outs.end());
								s = oi;
								assert(l == -1U);
								l = li;
							}
						}
					}
					assert(s != outs.end());
					assert(l < s->size());
					std::rotate(s->begin(), s->begin() + l, s->end());

					//move 's' to outs[0]:
					std::swap(*outs.begin(), *s);
					s = outs.begin();

					//store the pre-transform version of outs[0] as the step's "inter":
					step.inter = *s;

					assert(s->size() >= in_chain.size());
					for (uint32_t i = 0; i < in_chain.size(); ++i) {
						assert((*s)[i] == in_chain[i]);
					}

					s->erase(s->begin(), s->begin() + in_chain.size());
					s->insert(s->begin(), out_chain.begin(), out_chain.end());
					if (s->empty()) {
						outs.erase(s);
						//TODO: special flag for steps whose inter doesn't actually turn into outs[0] ?
					}
				}
				uint32_t base = storages.size();
				storages.insert(storages.end(), outs.begin(), outs.end());
				for (uint32_t i = base; i < storages.size(); ++i) {
					step.out.push_back(i);
				}
			}

			{ //DEBUG: dump step ins/outs
				auto storage_to_string = [&](uint32_t idx) {
					assert(idx < storages.size());
					std::string ret = std::to_string(idx) + "[";
					for (auto const& s : storages[idx]) {
						if (&s != &storages[idx][0]) ret += ' ';
						ret += s.to_string();
					}
					ret += "]";
					return ret;
				};
				//qDebug() << "step[" << (&step - &steps[0]) << "]:";
				//qDebug() << "  in:";
				for (auto i : step.in) qDebug() << ' ' << storage_to_string(i);
				//qDebug() << " int:";
				for (auto const& l : step.inter) qDebug() << ' ' << l.to_string();
				//qDebug() << " out:";
				for (auto o : step.out) qDebug() << ' ' << storage_to_string(o);
			}

			if (step.in.size() > 0 && step.out.size() > 0) {
				//make a map from step.inter -> step.out[0]
				std::vector< uint32_t > inter_to_out(step.inter.size(), -1U);

				std::map< Loop, CycleIndex > inter_loops;
				std::map< Loop, CycleIndex > out_loops;

				for (uint32_t oi = 0; oi < step.out.size(); ++oi) {
					Storage const& inter = (oi == 0 ? step.inter : storages[step.out[oi]]);
					for (uint32_t i = 0; i < inter.size(); ++i) {
						auto const& l = inter[i];
						auto ret = inter_loops.insert(std::make_pair(l, CycleIndex(oi, i)));
						assert(ret.second);
					}
					Storage const& out = storages[step.out[oi]];
					for (uint32_t i = 0; i < out.size(); ++i) {
						auto const& l = out[i];
						auto ret = out_loops.insert(std::make_pair(l, CycleIndex(oi, i)));
						assert(ret.second);
					}
				}

				//First, the loops that actually change because of stitches:
				for (uint32_t si = step.begin; si != step.end; ++si) {
					KnitOutStitch const& stitch = stitches[si];
					if (stitch.type == KnitOutStitch::Tuck || stitch.type == KnitOutStitch::Miss || stitch.type == KnitOutStitch::Knit
						|| stitch.type == KnitOutStitch::Increase) {
						auto f = inter_loops.find(source_loop(si, 0));
						assert(f != inter_loops.end());
						assert(f->second.cycle == 0 && f->second.index < inter_to_out.size());

						auto t = out_loops.find(Loop(si, 0));
						assert(t != out_loops.end());
						assert(t->second.cycle == 0 && t->second.index < storages[step.out[0]].size());

						assert(inter_to_out[f->second.index] == -1U);
						inter_to_out[f->second.index] = t->second.index;

					}
					else if (stitch.type == KnitOutStitch::Decrease) {
						auto t = out_loops.find(Loop(si, 0));
						assert(t != out_loops.end());
						assert(t->second.cycle == 0 && t->second.index < storages[step.out[0]].size());
						for (uint32_t i = 0; i < 2; ++i) {
							auto f = inter_loops.find(source_loop(si, i));
							assert(f != inter_loops.end());
							assert(f->second.cycle == 0 && f->second.index < inter_to_out.size());

							assert(inter_to_out[f->second.index] == -1U);
							inter_to_out[f->second.index] = t->second.index;
						}
					}
					else {
						//qDebug() << "What is '" << stitch.type << "'?"; //DEBUG
						assert(0 && "Unsupported stitch type.");
					}
				}

				//add things in the cycle that aren't touched by stitches:
				for (uint32_t ii = 0; ii < step.inter.size(); ++ii) {
					if (inter_to_out[ii] != -1U) continue;
					auto t = out_loops.find(step.inter[ii]);
					assert(t != out_loops.end());
					assert(t->second.cycle == 0 && t->second.index < storages[step.out[0]].size());
					inter_to_out[ii] = t->second.index;
				}

				step.inter_to_out = inter_to_out;
			}

			{
				if (step.end >= stitches.size() || stitches[step.end].yarn != stitches[step.end - 1].yarn || stitches[step.end - 1].type == KnitOutStitch::End) {
					auto f = active_yarns.find(stitches[step.begin].yarn);
					assert(f != active_yarns.end());
					active_yarns.erase(f);
					assert(active_yarns.empty()); //<-- should only have had one yarn active at a time anyway.
				}
				else {
					yarn.loop = Loop(
						step.end - 1,
						stitches[step.end - 1].out[1] != -1U ? 1 : 0
					);
					assert(stitches[step.end - 1].out[yarn.loop.idx] != -1U);
					yarn.direction = stitches[step.end - 1].direction;
				}
			}

			{ //update active_loops:
				for (StorageIdx i : step.in) {
					for (auto const& l : storages[i]) {
						auto f = active_loops.find(l);
						assert(f != active_loops.end());
						active_loops.erase(f);
					}
				}
				for (StorageIdx i : step.out) {
					for (auto const& l : storages[i]) {
						assert(!active_loops.count(l));
						active_loops[l] = i;
					}
				}
			}

		} //while (more stitches)

	} //end of build steps


	//Figure out possible shapes for storages near *interesting* steps:
	//qDebug() << "Figuring out shapes for interesting steps:"; //DEBUG
	for (auto& step : steps) {
		//interesting steps have more than one out/in:
		if (step.in.size() <= 1 && step.out.size() <= 1) continue;

		//qDebug() << "steps[" << (&step - &steps[0]) << "]:"; //DEBUG

		{ //DEBUG
			for (StorageIdx s : step.in) {
				qDebug() << "  uses storages[" << s << "]:";
				for (auto const& l : storages[s]) qDebug() << " " << l.to_string();
			}
			for (StorageIdx s : step.out) {
				qDebug() << "  makes storages[" << s << "]:";
				for (auto const& l : storages[s]) qDebug() << " " << l.to_string();
			}
		}

		//Construction model:
		// input cycles come in some order and shape
		// (input cycles are perhaps reshaped) <-- "late reshape"
		// middle (intermediate) cycle is created
		// intermediate cycle is reshaped
		// knitting happens


		//get local copies of storages:
		std::vector< Storage > ins, outs, intermediates;
		for (auto si : step.in) {
			ins.emplace_back(storages[si]);
		}
		for (auto si : step.out) {
			outs.emplace_back(storages[si]);
		}
		intermediates.emplace_back(step.inter);
		for (uint32_t i = 1; i < step.out.size(); ++i) {
			intermediates.emplace_back(storages[step.out[i]]);
		}

		//The contents of 'inter + outs[1...]' and 'ins' should be the same. Construct a mapping:
		std::vector< std::vector< CycleIndex > > in_to_inter;
		{
			std::map< Loop, CycleIndex > loops;
			for (auto const& inter : intermediates) {
				for (uint32_t i = 0; i < inter.size(); ++i) {
					auto const& l = inter[i];
					auto ret = loops.insert(std::make_pair(l, CycleIndex(&inter - &intermediates[0], i)));
					assert(ret.second);
				}
			}
			uint32_t used = 0;
			in_to_inter.reserve(ins.size());
			for (auto const& in : ins) {
				in_to_inter.emplace_back();
				in_to_inter.back().reserve(in.size());
				for (auto const& l : in) {
					auto f = loops.find(l);
					assert(f != loops.end());
					in_to_inter.back().emplace_back(f->second);
					++used;
				}
				assert(in_to_inter.back().size() == in.size());
			}
			assert(in_to_inter.size() == ins.size());
			assert(used == loops.size());
		}

		//intermediates[0] (== step.inter) gets knit on and changes contents before it is output:
		assert(step.inter_to_out.size() == step.inter.size());

		//TODO: ~cache of transfer costs~


		//enumerate input orders and shapes:
		//  -- for each order, test if the intermediate storages are all cycles
		//     then figure out cost to reshape for knitting (for all possible reshapes)
		//     add [cost, input order, input shapes, output shapes, output order] to possible shapes list
		std::vector< uint32_t > remaining;
		remaining.reserve(ins.size());
		for (uint32_t i = 0; i < ins.size(); ++i) {
			remaining.emplace_back(i);
		}
		std::vector< uint32_t > in_order;
		std::vector< PackedShape > in_shapes(ins.size(), -1U);
		std::vector< CycleIndex > front;
		std::vector< CycleIndex > back;

		//Given a current order, test to see if it leaves intermediates in a valid configuration:
		auto test_order = [&]() -> const char* {
			struct Info {
				uint32_t back_min = -1U;
				uint32_t back_max = 0;
				uint32_t front_min = -1U;
				uint32_t front_max = 0;
				uint32_t zero_index = -1U;
				bool zero_on_back = false;
			};
			std::vector< Info > inter_info(intermediates.size());
			for (uint32_t i = 0; i < back.size(); ++i) {
				auto const& ci = back[i];
				if (ci.cycle == -1U) continue; //gap
				assert(ci.cycle < inter_info.size());
				auto& info = inter_info[ci.cycle];
				info.back_min = std::min(info.back_min, i);
				info.back_max = std::max(info.back_max, i);
				if (ci.index == 0) {
					assert(info.zero_index == -1U);
					info.zero_index = i;
					info.zero_on_back = true;
				}
			}
			for (uint32_t i = 0; i < front.size(); ++i) {
				auto const& ci = front[i];
				if (ci.cycle == -1U) continue; //gap
				assert(ci.cycle < inter_info.size());
				auto& info = inter_info[ci.cycle];
				info.front_min = std::min(info.front_min, i);
				info.front_max = std::max(info.front_max, i);
				if (ci.index == 0) {
					assert(info.zero_index == -1U);
					info.zero_index = i;
					info.zero_on_back = false;
				}
			}

			std::vector< PackedShape > inter_shapes;
			inter_shapes.reserve(intermediates.size());
			for (auto const& info : inter_info) {
				uint32_t cycle = &info - &inter_info[0];
				//check for gaps inside cycle:
				uint32_t count = 0;
				if (info.back_min <= info.back_max) count += (info.back_max - info.back_min + 1);
				if (info.front_min <= info.front_max) count += (info.front_max - info.front_min + 1);
				if (count > intermediates[cycle].size()) {
					//gap found; bail out
					return "gap";
				}
				assert(count == intermediates[cycle].size()); //must have seen all loops
				assert(info.zero_index != -1U); //including the zero element, of course

				Shape shape(0, 0);

				//check that things are aligned:
				if (info.back_min > info.back_max) {
					//only on front bed
					assert(info.front_min <= info.front_max);
					if (info.front_min != info.front_max) {
						return "flat (front)";
					}
					shape.nibbles = Shape::BackLeft;
					shape.roll = 0;
				}
				else if (info.front_min > info.front_max) {
					//only on back bed
					assert(info.back_min <= info.back_max);
					if (info.back_min != info.back_max) {
						return "flat (back)";
					}
					shape.nibbles = Shape::FrontLeft;
					shape.roll = 0;
				}
				else {
					//on both beds, so check alignment between beds:
					if (std::abs(int32_t(info.back_min) - int32_t(info.front_min)) > 1) {
						return "mis-aligned (left)";
					}
					if (std::abs(int32_t(info.back_max) - int32_t(info.front_max)) > 1) {
						return "mis-aligned (right)";
					}
					shape.nibbles =
						(info.back_min > info.front_min ? Shape::BackLeft : 0)
						| (info.back_min < info.front_min ? Shape::FrontLeft : 0)
						| (info.back_max < info.front_max ? Shape::BackRight : 0)
						| (info.back_max > info.front_max ? Shape::FrontRight : 0)
						;
					if (info.zero_on_back) {
						shape.roll = (info.front_max - info.front_min + 1) + (info.back_max - info.zero_index);
					}
					else {
						shape.roll = info.zero_index - info.front_min;
					}
				}

				std::vector< CycleIndex > data;
				data.reserve(intermediates[cycle].size());
				for (uint32_t i = 0; i < intermediates[cycle].size(); ++i) {
					data.emplace_back(cycle, i);
				}
				std::vector< CycleIndex > expected_front, expected_back;
				shape.append_to_beds(data, CycleIndex(-1U, -1U), &expected_front, &expected_back);

				for (uint32_t i = info.back_min; i <= info.back_max; ++i) {
					uint32_t r = i - std::min(info.back_min, info.front_min);
					assert(r < expected_back.size());
					if (back[i] != expected_back[r]) {
						return "bad shape (back)";
					}
				}
				for (uint32_t i = info.front_min; i <= info.front_max; ++i) {
					uint32_t r = i - std::min(info.back_min, info.front_min);
					assert(r < expected_front.size());
					if (front[i] != expected_front[r]) {
						return "bad shape (front)";
					}
				}

				inter_shapes.emplace_back(shape.pack());
			}

			std::vector< uint32_t > out_order;
			out_order.reserve(outs.size());
			{ //figure out order of outs:
				std::vector< bool > seen(outs.size(), false);
				auto do_out = [&seen, &out_order](uint32_t cycle) {
					if (cycle == -1U) return;
					assert(cycle < seen.size());
					if (!seen[cycle]) {
						seen[cycle] = true;
						out_order.emplace_back(cycle);
					}
				};
				for (uint32_t i = 0; i < std::max(front.size(), back.size()); ++i) {
					if (i < front.size()) do_out(front[i].cycle);
					if (i < back.size()) do_out(back[i].cycle);
				}
				//make sure we saw all out cycles:
				for (auto s : seen) {
					assert(s);
				}
				assert(out_order.size() == outs.size());
			}


			//Order passed the test, now add possible transforms with their costs:


			{ //new method: no baked-in transfers (will bake into edges instead):
				//so just penalize for the intermediate shapes (might be double-charging?):
				ScheduleCost cost;
				for (auto const& inter_shape : inter_shapes) {
					cost += ScheduleCost::shape_cost(Shape::unpack(inter_shape));
				}

				ScheduleOption option;
				option.cost = cost;

				option.in_order = in_order;
				option.in_shapes = in_shapes;

				option.inter_shape = inter_shapes[0];

				option.out_order = out_order;
				option.out_shapes = inter_shapes;

				//out0 needs to be assigned a shape, so pick the cheapest one:
				std::vector< Shape > out0_shapes = Shape::make_shapes_for(outs[0].size());
				ScheduleCost best_cost = ScheduleCost::max();
				for (auto const& out0_shape : out0_shapes) {
					ScheduleCost cost = ScheduleCost::shape_cost(out0_shape);
					cost += ScheduleCost::transfer_cost(
						intermediates[0].size(), Shape::unpack(inter_shapes[0]),
						outs[0].size(), out0_shape,
						step.inter_to_out);
					if (cost < best_cost) {
						best_cost = cost;
						option.out_shapes[0] = out0_shape.pack();
					}
				}

				step.options.emplace_back(option);
			}

			return nullptr;

		};

		std::function< void() > enumerate_orders = [&]() {
			if (remaining.empty()) {		
				test_order();
				return;
			}
			//order isn't done yet, try all shapes:
			for (uint32_t r = 0; r < remaining.size(); ++r) {
				std::swap(remaining[r], remaining.back());
				uint32_t idx = remaining.back();
				remaining.pop_back();
				in_order.emplace_back(idx);

				std::vector< Shape > shapes = Shape::make_shapes_for(ins[idx].size());
				for (auto const& shape : shapes) {
					in_shapes[idx] = shape.pack();

					uint32_t front_size_before = front.size();
					uint32_t back_size_before = back.size();

					shape.append_to_beds(in_to_inter[idx], CycleIndex(-1U, -1U), &front, &back);

					enumerate_orders();

					front.erase(front.begin() + front_size_before, front.end());
					back.erase(back.begin() + back_size_before, back.end());

					in_shapes[idx] = -1U;
				}

				assert(!in_order.empty() && in_order.back() == idx);
				in_order.pop_back();
				remaining.push_back(idx);
				std::swap(remaining[r], remaining.back());
			}
		};

		enumerate_orders();

		//qDebug() << "In total, step had " << step.options.size() << " scheduling options.";

	} //interesting steps

	//qDebug() << "(done)" ;


	//Figure out possible shapes for storages near *boring* steps:
	qDebug() << "Figuring out shapes for boring steps:"; //DEBUG
	for (auto& step : steps) {
		//qDebug() << "step = " << &step - &steps[0];
		if (!(step.in.size() <= 1 && step.out.size() <= 1)) continue; //skip exciting steps
		uint32_t inter_roll = 0; //such that: storages[step.in[0]][i] == step.inter[i + inter_roll]
		//used to fix up inter_shape relative to in_shape so the stitches are in the same places.
		//
		if (step.in.size() == 1) {
			assert(step.inter.size() == storages[step.in[0]].size());
			inter_roll = -1U;
			for (uint32_t roll = 0; roll < step.inter.size(); ++roll) {
				if (storages[step.in[0]][0] == step.inter[roll]) {
					assert(inter_roll == -1U);
					inter_roll = roll;
				}
			}
			assert(inter_roll != -1U);
			for (uint32_t i = 0; i < step.inter.size(); ++i) {
				assert(storages[step.in[0]][i] == step.inter[(i + inter_roll) % step.inter.size()]);

			}
		}
		else {
			assert(step.inter.empty());
		}

		if (step.in.size() == 0 && step.out.size() == 1) {
			//Starts a cycle -- can do so in any order.
			assert(stitches[step.begin].type == KnitOutStitch::Start);

			std::vector< Shape > out_shapes = Shape::make_shapes_for(storages[step.out[0]].size());
			ScheduleOption option;
			option.in_order.clear();
			option.in_shapes.clear();
			option.inter_shape = 0;
			option.out_order = { 0 };
			for (auto const& out_shape : out_shapes) {
				option.cost = ScheduleCost::shape_cost(out_shape);
				option.out_shapes = { out_shape.pack() };
				step.options.emplace_back(option);
			}

		}
		else if (step.in.size() == 1 && step.out.size() == 0) {
			//Ends a cycle -- can do so in any order.
			assert(stitches[step.begin].type == KnitOutStitch::End);

			std::vector< Shape > in_shapes = Shape::make_shapes_for(storages[step.in[0]].size());
			ScheduleOption option;
			option.out_order.clear();
			option.out_shapes.clear();
			option.in_order = { 0 };
			for (auto const& in_shape : in_shapes) {
				option.cost = ScheduleCost::shape_cost(in_shape);
				Shape inter_shape = in_shape;
				inter_shape.roll = (in_shape.roll + step.inter.size() - inter_roll) % step.inter.size();
				option.inter_shape = inter_shape.pack();
				option.in_shapes = { in_shape.pack() };

				step.options.emplace_back(option);
			}

		}
		else {
			assert(step.in.size() == 1 && step.out.size() == 1);

			std::vector< Shape > in_shapes = Shape::make_shapes_for(storages[step.in[0]].size());
			std::vector< Shape > out_shapes = Shape::make_shapes_for(storages[step.out[0]].size());
			assert(step.inter.size() == storages[step.in[0]].size());
			assert(step.inter.size() == step.inter_to_out.size());

			ScheduleOption option;
			option.in_order = { 0 };
			option.out_order = { 0 };
			for (auto const& in_shape : in_shapes) {
				option.in_shapes = { in_shape.pack() };
				Shape inter_shape = in_shape;
				inter_shape.roll = (in_shape.roll + step.inter.size() - inter_roll) % step.inter.size();
				option.inter_shape = inter_shape.pack();
				for (auto const& out_shape : out_shapes) {
					option.out_shapes = { out_shape.pack() };

					option.cost = ScheduleCost::shape_cost(in_shape);
					option.cost += ScheduleCost::shape_cost(out_shape);
					option.cost += ScheduleCost::transfer_cost(
						step.inter.size(), inter_shape,
						storages[step.out[0]].size(), out_shape,
						step.inter_to_out
					);
					step.options.emplace_back(option);
				}
			}

		}
		//qDebug() << "steps[" << (&step - &steps[0]) << "] had " << step.options.size() << " scheduling options.";
	}

	std::vector< int32_t > storage_positions(storages.size(), std::numeric_limits< int32_t >::max());
	std::vector< PackedShape > storage_shapes(storages.size(), -1U); //shapes of storages just after creation (should be same as the step_options of the step that made them)
	std::vector< uint32_t > step_options(steps.size(), -1U); //shapes of in and out cycles
	{ //Next up: do upward-planar embedding.

		//For every 'interesting' step, build a DAGNode
		//For every chain of 1-in / 1-out steps, build a DAGEdge

		std::vector< DAGNode > nodes;
		std::vector< DAGEdge > edges;

		std::map< uint32_t, DAGEdgeIndex > active_edges; //maps from storages to edges currently in progress

		std::vector< uint32_t > node_steps;

		//dag edges corresponding to each step's ins/outs:
		std::vector< std::vector< DAGEdgeIndex > > step_in_edges;
		std::vector< std::vector< DAGEdgeIndex > > step_out_edges;
		step_in_edges.reserve(steps.size());
		step_out_edges.reserve(steps.size());

		//qDebug() << "Building DAG "; 
		emit helpBoxCommunication("Building DAG, this may take a while...");
		for (auto const& step : steps) {
			//qDebug() << "."; 
			step_in_edges.emplace_back();
			step_out_edges.emplace_back();
			if (step.in.size() == 1 && step.out.size() == 1) {
				//1-1 step; accumulate DAGEdge cost matrix:
				auto f = active_edges.find(step.in[0]);
				assert(f != active_edges.end());
				assert(f->second < edges.size());
				DAGEdge& edge = edges[f->second];
				assert(edge.to == step.in[0]);
				active_edges.erase(f);

				std::vector< DAGCost > old_costs = std::move(edge.costs);
				uint32_t old_to_shapes = edge.to_shapes;
				assert(edge.to_shapes == Shape::count_shapes_for(storages[step.in[0]].size()));

				edge.to = step.out[0];
				edge.to_shapes = Shape::count_shapes_for(storages[step.out[0]].size());
				edge.costs.resize(edge.from_shapes * edge.to_shapes, ScheduleCost::max());

				//accumulate step costs into edge costs:
				for (auto const& option : step.options) {
					assert(option.in_order.size() == 1 && option.in_order[0] == 0);
					assert(option.in_shapes.size() == 1);
					assert(option.out_order.size() == 1 && option.out_order[0] == 0);
					assert(option.out_shapes.size() == 1);

					uint32_t in_shape = Shape::unpack(option.in_shapes[0]).index_for(storages[step.in[0]].size());
					uint32_t out_shape = Shape::unpack(option.out_shapes[0]).index_for(storages[step.out[0]].size());

					for (uint32_t f = 0; f < edge.from_shapes; ++f) {
						DAGCost& cost = edge.costs[f * edge.to_shapes + out_shape];
						DAGCost test = old_costs[f * old_to_shapes + in_shape];
						if (!(test == DAGCost::max())) {
							test += option.cost;
							cost = std::min(cost, test);
						}
					}
				}

				//step is <within> this edge:
				step_in_edges.back().emplace_back(&edge - &edges[0]);
				step_out_edges.back().emplace_back(&edge - &edges[0]);

				auto ret = active_edges.insert(std::make_pair(step.out[0], &edge - &edges[0]));
				assert(ret.second);
			}
			else {
				assert(step.in.size() != 1 || step.out.size() != 1);
				//Exciting step; build a DAGNode, consume some DAGEdges, start some DAGEdges.

				node_steps.emplace_back(&step - &steps[0]);
				nodes.emplace_back();
				DAGNode& node = nodes.back();

				std::vector< DAGEdgeIndex > in_edges;
				std::vector< DAGEdgeIndex > out_edges;

				//look up indices of in edges:
				for (uint32_t i = 0; i < step.in.size(); ++i) {
					auto f = active_edges.find(step.in[i]);
					assert(f != active_edges.end());
					DAGEdge& edge = edges[f->second];
					assert(edge.to == step.in[i]);
					edge.to = &node - &nodes[0];
					in_edges.emplace_back(f->second);
					active_edges.erase(f);
				}

				//create out edges:
				for (uint32_t i = 0; i < step.out.size(); ++i) {
					auto ret = active_edges.insert(std::make_pair(step.out[i], edges.size()));
					assert(ret.second);
					out_edges.emplace_back(edges.size());

					edges.emplace_back();
					DAGEdge& edge = edges.back();
					edge.from = &node - &nodes[0];
					edge.from_shapes = Shape::count_shapes_for(storages[step.out[i]].size());
					edge.to = step.out[i];
					edge.to_shapes = Shape::count_shapes_for(storages[step.out[i]].size());

					edge.costs.assign(edge.from_shapes * edge.to_shapes, ScheduleCost::max());
					assert(edge.from_shapes == edge.to_shapes);
					if (step.out.size() > 1 || step.in.size() > 1) {
						std::vector< Shape > shapes = Shape::make_shapes_for(storages[step.out[i]].size());
						assert(shapes.size() == edge.from_shapes);
						assert(shapes.size() == edge.to_shapes);
						std::vector< uint32_t > map(storages[step.out[i]].size());
						for (auto& m : map) {
							m = &m - &map[0];
						}
						//interesting steps get a ("free") xfer pass on their output. Maybe all steps should? Hmm.
						for (uint32_t f = 0; f < edge.from_shapes; ++f) {
							for (uint32_t t = 0; t < edge.to_shapes; ++t) {
								edge.costs[f * edge.to_shapes + t] = ScheduleCost::transfer_cost(
									storages[step.out[i]].size(), shapes[f],
									storages[step.out[i]].size(), shapes[t],
									map);
							}
						}
					}
					else {
						//NOTE: this ends up costing a bit of extra time, because multiple-source-shape edges always need a square cost matrix, but we don't actually have any costs to put into it, so we just make up a matrix that does not permit rotation:
						for (uint32_t s = 0; s < edge.from_shapes; ++s) {
							edge.costs[s * edge.to_shapes + s] = ScheduleCost::zero();
						}
					}
				}

				step_in_edges.back() = in_edges;
				step_out_edges.back() = out_edges;

				//PARANOIA: no duplicate in_edges or out_edges, right?
				assert(std::set< uint32_t >(in_edges.begin(), in_edges.end()).size() == in_edges.size());
				assert(std::set< uint32_t >(out_edges.begin(), out_edges.end()).size() == out_edges.size());

				node.options.reserve(step.options.size());
				for (auto const& option : step.options) {
					node.options.emplace_back();
					DAGOption& dag_option = node.options.back();
					dag_option.cost = option.cost;
					dag_option.in_order.reserve(step.in.size());
					dag_option.in_shapes.reserve(step.in.size());
					dag_option.out_order.reserve(step.out.size());
					dag_option.out_shapes.reserve(step.out.size());
					for (uint32_t i = 0; i < step.in.size(); ++i) {
						uint32_t in = option.in_order[i];
						assert(in < in_edges.size());
						dag_option.in_order.emplace_back(in_edges[in]);
						assert(in < step.in.size());
						dag_option.in_shapes.emplace_back(Shape::unpack(option.in_shapes[in]).index_for(storages[step.in[in]].size()));
						assert(dag_option.in_shapes.back() < edges[in_edges[in]].to_shapes); //DEBUG
					}
					assert(dag_option.in_order.size() == step.in.size());
					assert(dag_option.in_shapes.size() == step.in.size());
					for (uint32_t o = 0; o < step.out.size(); ++o) {
						uint32_t out = option.out_order[o];
						assert(out < out_edges.size());
						dag_option.out_order.emplace_back(out_edges[out]);
						assert(out < step.out.size());
						assert(out < option.out_shapes.size());
						dag_option.out_shapes.emplace_back(Shape::unpack(option.out_shapes[out]).index_for(storages[step.out[out]].size())); //<-- WHY IS THIS FAILING?
						assert(dag_option.out_shapes.back() < edges[out_edges[out]].from_shapes); //DEBUG
					}
					assert(dag_option.out_order.size() == step.out.size());
					assert(dag_option.out_shapes.size() == step.out.size());
					//PARANOIA: no duplicates in order, right?
					assert(std::set< uint32_t >(dag_option.in_order.begin(), dag_option.in_order.end()).size() == dag_option.in_order.size());
					assert(std::set< uint32_t >(dag_option.out_order.begin(), dag_option.out_order.end()).size() == dag_option.out_order.size());
				}
				assert(node.options.size() == step.options.size());
			}
		} //for(steps)

		
		//qDebug() << " done.";

		assert(node_steps.size() == nodes.size());

		
		emit helpBoxCommunication("DAG built! Have "+QString::number(edges.size())+" edges and "+QString::number(nodes.size()) + " nodes.\n" +
			"Embedding the graph, this may take a while...");
		std::vector< uint32_t > node_options;
		std::vector< int32_t > node_positions, edge_positions;

		if (!embed_DAG(nodes, edges, &node_options, &node_positions, &edge_positions)) {
			emit helpBoxCommunication( "ERROR: failed to find an upward-planar embedding.");
			return;
		}

		assert(node_options.size() == nodes.size());
		assert(node_positions.size() == nodes.size());
		assert(edge_positions.size() == edges.size());

		emit helpBoxCommunication("Going through scheduling steps, this make take a while...");
		//Go through steps and assign selected option:
		assert(node_options.size() == nodes.size());
		for (uint32_t n = 0; n < nodes.size(); ++n) {
			assert(node_steps[n] < steps.size());
			//WARNING COST THING HERE
			//qDebug() << "Step " << node_steps[n] << " gets option " << node_options[n] << " of " << nodes[n].options.size() << " for cost ";
			assert(node_options[n] < steps[node_steps[n]].options.size());
			assert(nodes[n].options.size() == steps[node_steps[n]].options.size());
			step_options[node_steps[n]] = node_options[n];
		}
		for (auto const& e : edges) {
			//qDebug() << "Edge " << (&e - &edges[0]) << " from " << e.from << " to " << e.to << " ends up with cost "; 
			uint32_t from_shape = -1U;
			auto const& from_option = nodes[e.from].options[node_options[e.from]];
			for (uint32_t out = 0; out < from_option.out_order.size(); ++out) {
				if (from_option.out_order[out] == (&e - &edges[0])) {
					assert(from_shape == -1U);
					from_shape = from_option.out_shapes[out];
				}
			}
			uint32_t to_shape = -1U;
			auto const& to_option = nodes[e.to].options[node_options[e.to]];
			for (uint32_t in = 0; in < to_option.in_order.size(); ++in) {
				if (to_option.in_order[in] == (&e - &edges[0])) {
					assert(to_shape == -1U);
					to_shape = to_option.in_shapes[in];
				}
			}
			assert(from_shape != -1U);
			assert(to_shape != -1U);

			assert(from_shape < e.from_shapes);
			assert(to_shape < e.to_shapes);
			//qDebug() << e.costs[from_shape * e.to_shapes + to_shape] << std::endl;

		}

		//qDebug() << "assigning storage positions";
		//assign storage positions:
		for (auto const& step : steps) {
			uint32_t si = &step - &steps[0];
			if (step_in_edges[si].size() == 1 && step_out_edges[si].size() == 1) {
				assert(step_options[si] == -1U); //wasn't part of the dag nodes.
				assert(step_in_edges[si][0] == step_out_edges[si][0]);
			}
			assert(step_in_edges[si].size() == step.in.size());
			for (uint32_t in = 0; in < step.in.size(); ++in) {
				assert(step_in_edges[si][in] < edge_positions.size());
				assert(storage_positions[step.in[in]] == edge_positions[step_in_edges[si][in]]);
			}
			assert(step_out_edges[si].size() == step.out.size());
			for (uint32_t out = 0; out < step.out.size(); ++out) {
				assert(step_out_edges[si][out] < edge_positions.size());
				assert(storage_positions[step.out[out]] == std::numeric_limits< int32_t >::max());
				storage_positions[step.out[out]] = edge_positions[step_out_edges[si][out]];
			}
		}

		//build up cost matrices to figure out all unassigned options:
		struct ShapeCost {
			ScheduleCost cost = ScheduleCost::max();
			uint32_t option = -1U; //option that cost comes via
		};
		std::vector< std::vector< ShapeCost > > storage_costs;
		storage_costs.reserve(storages.size());
		for (auto const& storage : storages) {
			storage_costs.emplace_back(Shape::count_shapes_for(storage.size()));
		}
		std::vector< uint32_t > storage_steps(storages.size(), -1U);
		for (auto const& step : steps) {
			for (auto o : step.out) {
				assert(o < storage_steps.size());
				assert(storage_steps[o] == -1U);
				storage_steps[o] = &step - &steps[0];
			}
		}

		//qDebug() << "iterating thru steps";
		for (auto const& step : steps) {
			uint32_t si = &step - &steps[0];
			if (step_options[si] != -1U) {
				//already-assigned option:
				auto const& option = step.options[step_options[si]];
				//chase back storage_costs chains from ins:
				for (uint32_t in = 0; in < step.in.size(); ++in) {
					uint32_t storage = step.in[in];
					uint32_t idx = Shape::unpack(option.in_shapes[in]).index_for(storages[step.in[in]].size());
					while (true) {
						uint32_t si = storage_steps[storage];
						assert(si < steps.size());
						if (step_options[si] != -1U) break; //finished
						assert(steps[si].in.size() == 1 && steps[si].out.size() == 1);
						assert(steps[si].out[0] == storage);

						assert(!(storage_costs[storage][idx].cost == ScheduleCost::max()));
						uint32_t opt = storage_costs[storage][idx].option;
						assert(opt < steps[si].options.size());

						uint32_t in_idx = Shape::unpack(steps[si].options[opt].in_shapes[0]).index_for(storages[steps[si].in[0]].size());
						uint32_t out_idx = Shape::unpack(steps[si].options[opt].out_shapes[0]).index_for(storages[steps[si].out[0]].size());
						assert(out_idx == idx);

						step_options[si] = opt;

						idx = in_idx;
						storage = steps[si].in[0];
					}
				}
				//start storage_costs chains from outs:
				for (uint32_t out = 0; out < step.out.size(); ++out) {
					uint32_t idx = Shape::unpack(option.out_shapes[out]).index_for(storages[step.out[out]].size());
					assert(idx < storage_costs[step.out[out]].size());
					storage_costs[step.out[out]][idx].cost = ScheduleCost::zero();
					storage_costs[step.out[out]][idx].option = step_options[si];
				}
			}
			else {
				//was part of an edge, need to propagate:
				assert(step.in.size() == 1 && step.out.size() == 1);

				for (auto const& option : step.options) {
					assert(option.in_order.size() == 1);
					assert(option.in_shapes.size() == 1);
					assert(option.out_order.size() == 1);
					assert(option.out_shapes.size() == 1);

					uint32_t in_shape = Shape::unpack(option.in_shapes[0]).index_for(storages[step.in[0]].size());
					uint32_t out_shape = Shape::unpack(option.out_shapes[0]).index_for(storages[step.out[0]].size());

					assert(in_shape < storage_costs[step.in[0]].size());
					ScheduleCost test = storage_costs[step.in[0]][in_shape].cost;
					if (!(test == ScheduleCost::max())) {
						test += option.cost;
						assert(out_shape < storage_costs[step.out[0]].size());
						if (test < storage_costs[step.out[0]][out_shape].cost) {
							storage_costs[step.out[0]][out_shape].cost = test;
							storage_costs[step.out[0]][out_shape].option = &option - &step.options[0];
						}
					}
				}
			}
		}


		//record shapes:
		//qDebug() << "record shapes";
		uint32_t pre_xfers = 0;
		for (uint32_t si = 0; si < steps.size(); ++si) {
			//qDebug() << "Step " << si << " option " << step_options[si] << " says "; //DEBUG
			assert(step_options[si] < steps[si].options.size());
			auto const& option = steps[si].options[step_options[si]];
			for (uint32_t in = 0; in < steps[si].in.size(); ++in) {
				//qDebug() << " i" << steps[si].in[in] << "=" << option.in_shapes[in]; //DEBUG
				assert(storage_shapes[steps[si].in[in]] != -1U);
				if (storage_shapes[steps[si].in[in]] != option.in_shapes[in]) {
					++pre_xfers;
					//qDebug() << "*"; //DEBUG
				}
			}
			for (uint32_t out = 0; out < steps[si].out.size(); ++out) {
				//qDebug() << " o" << steps[si].out[out] << "=" << option.out_shapes[out]; //DEBUG
				assert(storage_shapes[steps[si].out[out]] == -1U);
				storage_shapes[steps[si].out[out]] = option.out_shapes[out];
			}
			//endline here
			//PARANOIA: check order vs storage positions:
			for (uint32_t i = 1; i < option.in_order.size(); ++i) {
				assert(option.in_order[i - 1] < steps[si].in.size());
				assert(option.in_order[i] < steps[si].in.size());
				assert(
					storage_positions[steps[si].in[option.in_order[i - 1]]]
					< storage_positions[steps[si].in[option.in_order[i]]]
				);
			}
			for (uint32_t i = 1; i < steps[si].options[step_options[si]].out_order.size(); ++i) {
				assert(option.out_order[i - 1] < steps[si].out.size());
				assert(option.out_order[i] < steps[si].out.size());
				assert(
					storage_positions[steps[si].out[option.out_order[i - 1]]]
					< storage_positions[steps[si].out[option.out_order[i]]]
				);
			}
		}
		if (pre_xfers > 0) {
			qDebug() << "NOTE: will need to reshape " << pre_xfers << " inputs before steps.";
		}

		//TODO: PARANOIA: check that step option orders match recorded storage positions.
		
		auto summarize = [](Shape const& shape, uint32_t size, uint32_t mark) -> std::string {
			//shape as braille dots
			// e.g.  ⠲⠶⠷⠶
			std::vector< char > data(size, '.');
			if (mark < size) data[mark] = '*';
			std::vector< char > front, back;
			shape.append_to_beds(data, ' ', &front, &back);

			if (front.size() % 2 == 1) front.emplace_back(' ');
			if (back.size() % 2 == 1) back.emplace_back(' ');

			std::vector< uint8_t > dots(std::max(front.size(), back.size()) / 2, 0);
			for (uint32_t i = 0; i + 1 < front.size(); i += 2) {
				if (front[i] == '.') dots[i / 2] |= 0x4;
				if (front[i + 1] == '.') dots[i / 2] |= 0x20;
				if (front[i] == '*') dots[i / 2] |= 0x4 | 0x40;
				if (front[i + 1] == '*') dots[i / 2] |= 0x20 | 0x80;
			}
			for (uint32_t i = 0; i + 1 < back.size(); i += 2) {
				if (back[i] == '.') dots[i / 2] |= 0x2;
				if (back[i + 1] == '.') dots[i / 2] |= 0x10;
				if (back[i] == '*') dots[i / 2] |= 0x2 | 0x1;
				if (back[i + 1] == '*') dots[i / 2] |= 0x10 | 0x8;
			}

			std::string ret;
			for (auto d : dots) {
				uint16_t c = 0x2800 | d;
				ret += char(0xe0 | ((c >> 12) & 0x1f));
				ret += char(0x80 | ((c >> 6) & 0x3f));
				ret += char(0x80 | ((c) & 0x3f));
			}
			return ret;
		};

		//qDebug() << "Storage positions";
		//DEBUG: dump selections:
		for (auto const& step : steps) {
			uint32_t si = &step - &steps[0];
			qDebug() << "Step " << si << " gets option " << step_options[si] << " of " << steps[si].options.size();
			auto const& option = steps[si].options[step_options[si]];
			qDebug() << "  takes";
			for (auto i : option.in_order) {
				uint32_t in = step.in[i];
				uint32_t mark = 0;
				for (uint32_t s = 1; s < storages[in].size(); ++s) {
					if (storages[in][mark] < storages[in][s]) mark = s;
				}
				qDebug() << " s" << in << " in shape " << storage_shapes[in] << " " << summarize(Shape::unpack(storage_shapes[in]), storages[in].size(), mark) << " in lane " << storage_positions[in] << ",";
			}
			qDebug() << "  makes";
			for (auto o : option.out_order) {
				uint32_t out = step.out[o];
				uint32_t mark = 0;
				for (uint32_t s = 1; s < storages[out].size(); ++s) {
					if (storages[out][mark] < storages[out][s]) mark = s;
				}
				qDebug() << " s" << out << " in shape " << storage_shapes[out] << " " << summarize(Shape::unpack(storage_shapes[out]), storages[out].size(), mark) << " in lane " << storage_positions[out] << ",";
			}
		}
	}

	struct StorageLeft {
		StorageLeft(uint32_t storage_, int32_t left_) : storage(storage_), left(left_) { }
		uint32_t storage;
		int32_t left;
	};
	std::vector< std::vector< StorageLeft > > step_storages; //left-edge positions for active storages just after the step
	std::vector< uint32_t > storage_widths; //widths of each storage in their just-after-step shapes.

	{ //Determine left-edge position for each active storage just after each step:

		//this will be handy; max width (over front/back bed) of each storage:
		storage_widths.reserve(storages.size());
		for (auto const& storage : storages) {
			Shape shape = Shape::unpack(storage_shapes[&storage - &storages[0]]);
			std::vector< char > stitches(storage.size(), '.');
			std::vector< char > front, back;
			shape.append_to_beds(stitches, ' ', &front, &back);
			storage_widths.emplace_back(std::max(front.size(), back.size()));
		}

		//Now, for each stitch, record:
		std::vector< std::vector< uint32_t > > active; //indices of active (sorted left-to-right)
		std::vector< std::vector< int32_t > > gaps; //spaces between active
		active.reserve(steps.size());
		gaps.reserve(steps.size());
		{ //determine what is active and start with a non-intersecting layout:
			std::set< uint32_t > current;
			for (auto const& step : steps) {
				for (auto in : step.in) {
					auto f = current.find(in);
					assert(f != current.end());
					current.erase(f);
				}
				for (auto out : step.out) {
					auto ret = current.insert(out);
					assert(ret.second);
				}
				std::vector< uint32_t > act(current.begin(), current.end());
				std::sort(act.begin(), act.end(), [&](uint32_t a, uint32_t b) {
					return storage_positions[a] < storage_positions[b];
					});
				gaps.emplace_back(std::max(int32_t(act.size()) - 1, 0), 0);
				active.emplace_back(std::move(act));
			}
			assert(current.empty());
		}

		//TODO: build penalties between active storages at each step by looking at stitch ancestors

		//TODO: some sort of actual optimization or something!
		for (uint32_t si = 0; si < active.size(); ++si) {
			gaps[si].assign(gaps[si].size(), 1); //just a 'lil space between everything
		}

		std::vector< std::vector< int32_t > > lefts; //left edge of each active storage
		lefts.reserve(steps.size());
		for (uint32_t si = 0; si < steps.size(); ++si) {
			std::vector< int32_t > left;
			if (!active[si].empty()) {
				left.emplace_back(0);
				for (uint32_t i = 0; i < gaps[si].size(); ++i) {
					left.emplace_back(left.back() + storage_widths[active[si][i]] + gaps[si][i]);
				}
			}
			//TODO: shift storage left positions based on penalties
			lefts.emplace_back(std::move(left));
		}

		step_storages.reserve(steps.size());
		for (uint32_t si = 0; si < steps.size(); ++si) {
			assert(lefts[si].size() == active[si].size());
			step_storages.emplace_back();
			auto& step_storage = step_storages.back();
			step_storage.reserve(active[si].size());
			for (uint32_t i = 0; i < active[si].size(); ++i) {
				step_storage.emplace_back(active[si][i], lefts[si][i]);
			}
		}

		{ //DEBUG: show where the steps sit on actual needles
			int32_t min_needle = std::numeric_limits< int32_t >::max();
			int32_t max_needle = std::numeric_limits< int32_t >::min();
			for (uint32_t si = 0; si < steps.size(); ++si) {
				if (active[si].empty()) continue;
				min_needle = std::min(min_needle, lefts[si][0]);
				max_needle = std::max(max_needle, lefts[si].back() + int32_t(storage_widths[active[si].back()]) - 1);
			}
			for (uint32_t si = 0; si < steps.size(); ++si) {
				std::string num = std::to_string(si);
				while (num.size() < 4) num = ' ' + num;
				//qDebug() << "after " << num << ": ";

				std::vector< bool > on_back(max_needle - min_needle + 1, false);
				std::vector< bool > on_front(max_needle - min_needle + 1, false);

				for (uint32_t a = 0; a < active[si].size(); ++a) {
					Shape shape = Shape::unpack(storage_shapes[active[si][a]]);
					std::vector< bool > data(storages[active[si][a]].size(), true);
					std::vector< bool > back, front;
					shape.append_to_beds(data, false, &front, &back);
					for (uint32_t i = 0; i < back.size(); ++i) {
						if (back[i]) {
							assert(on_back[i + (lefts[si][a] - min_needle)] == false);
							on_back[i + (lefts[si][a] - min_needle)] = true;
						}
					}
					for (uint32_t i = 0; i < front.size(); ++i) {
						if (front[i]) {
							assert(on_front[i + (lefts[si][a] - min_needle)] == false);
							on_front[i + (lefts[si][a] - min_needle)] = true;
						}
					}
				}

				if (on_front.size() % 2 == 1) on_front.emplace_back(false);
				if (on_back.size() % 2 == 1) on_back.emplace_back(false);

				std::vector< uint8_t > dots(std::max(on_front.size(), on_back.size()) / 2, 0);
				for (uint32_t i = 0; i + 1 < on_front.size(); i += 2) {
					if (on_front[i]) dots[i / 2] |= 0x4;
					if (on_front[i + 1]) dots[i / 2] |= 0x20;
					//if (front[i]   == '*') dots[i/2] |= 0x4 | 0x40;
					//if (front[i+1] == '*') dots[i/2] |= 0x20 | 0x80;
				}
				for (uint32_t i = 0; i + 1 < on_back.size(); i += 2) {
					if (on_back[i]) dots[i / 2] |= 0x2;
					if (on_back[i + 1]) dots[i / 2] |= 0x10;
					//if (back[i]   == '*') dots[i/2] |= 0x2 | 0x1;
					//if (back[i+1] == '*') dots[i/2] |= 0x10 | 0x8;
				}

				for (auto d : dots) {
					uint16_t c = 0x2800 | d;
					qDebug() << char(0xe0 | ((c >> 12) & 0x1f));
					qDebug() << char(0x80 | ((c >> 6) & 0x3f));
					qDebug() << char(0x80 | ((c) & 0x3f));
				}

				qDebug() << " |";
				for (auto l : lefts[si]) {
					qDebug() << " " << l;
				}
				qDebug();
			}
		}
	}

	//--- Dump knitting instructions steps ----
	std::vector< std::string > instructions;
	auto add_instr = [&instructions](std::string const& instr) {
		instructions.emplace_back(instr);
		//std::cout << instr << std::endl; //DEBUG
	};

	//Not necessary since i rewrote the javascript
	//add_instr("const autoknit = require('autoknit');");
	//add_instr("let h = new autoknit.Helpers;");
	emit helpBoxCommunication("Instructions scheduled! now generating knitout instructions...");
	//helper:
	auto typeset_bed_needle = [](char bed, int32_t needle) -> std::string {
		std::string ret;
		//ret += '\'';
		static_assert(BedNeedle::Front == 'f' && BedNeedle::Back == 'b', "BedNeedle and our ad-hoc bed characters agree.");
		if (bed == 'f' || bed == 'b') ret += bed;
		else if (bed == BedNeedle::FrontSliders) ret += "fs";
		else if (bed == BedNeedle::BackSliders) ret += "bs";
		else assert(0 && "unhandled bed name");
		ret += std::to_string(needle);
		//ret += '\'';
		return ret;
	};

	//keep track of every created loop:
	//std::unordered_multimap< BedNeedle, Loop > bn_to_loop; <-- someday, update this incrementally...
	std::map< Loop, BedNeedle > loop_to_bn;

	auto make_start = [&loop_to_bn](char bed0, int32_t needle0, Loop const& out0) {
		BedNeedle bn0(bed0 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle0);
		for (auto const& lbn : loop_to_bn) {
			assert(lbn.second != bn0); //needle should be vacant (NOTE: inefficient check!)
		}
		auto ret = loop_to_bn.insert(std::make_pair(out0, bn0));
		assert(ret.second); //loop shouldn't be in map already!
	};

	auto make_stitch = [&loop_to_bn](char bed0, int32_t needle0, Loop const& in0, Loop const& out0) {
		BedNeedle bn0(bed0 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle0);

		{ //remove input loop:
			auto f = loop_to_bn.find(in0);
			assert(f != loop_to_bn.end());
			assert(f->second == bn0); //input to stitch should be on this needle
			loop_to_bn.erase(f);
		}

		for (auto const& lbn : loop_to_bn) {
			assert(lbn.second != bn0); //needle should be vacant now that input is removed (NOTE: inefficient check!)
		}

		auto ret = loop_to_bn.insert(std::make_pair(out0, bn0));
		assert(ret.second); //output loop shouldn't be in map already!
	};

	auto make_increase = [&loop_to_bn](char bed0, int32_t needle0, char bed1, int32_t needle1, Loop const& in0, Loop const& out0, Loop const& out1) {
		BedNeedle bn0(bed0 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle0);
		BedNeedle bn1(bed1 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle1);

		{ //remove input loop:
			auto f = loop_to_bn.find(in0);
			assert(f != loop_to_bn.end());
			assert(f->second == bn0); //input to stitch should be on this needle
			loop_to_bn.erase(f);
		}

		for (auto const& lbn : loop_to_bn) {
			assert(lbn.second != bn0); //needle should be vacant now that input is removed (NOTE: inefficient check!)
			assert(lbn.second != bn1); //other needle should also be vacant now
		}

		{
			auto ret = loop_to_bn.insert(std::make_pair(out0, bn0));
			assert(ret.second); //output loop shouldn't be in map already!
		}

		{
			auto ret = loop_to_bn.insert(std::make_pair(out1, bn1));
			assert(ret.second); //output loop shouldn't be in map already!
		}
	};

	auto make_decrease = [&loop_to_bn](char bed0, int32_t needle0, Loop const& in0, Loop const& in1, Loop const& out0) {
		BedNeedle bn0(bed0 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle0);

		{ //remove input loop 0:
			auto f = loop_to_bn.find(in0);
			assert(f != loop_to_bn.end());
			assert(f->second == bn0); //input to stitch should be on this needle
			loop_to_bn.erase(f);
		}
		{ //remove input loop 1:
			auto f = loop_to_bn.find(in1);
			assert(f != loop_to_bn.end());
			assert(f->second == bn0); //input to stitch should be on this needle
			loop_to_bn.erase(f);
		}

		for (auto const& lbn : loop_to_bn) {
			assert(lbn.second != bn0); //needle should be vacant now that input is removed (NOTE: inefficient check!)
		}

		auto ret = loop_to_bn.insert(std::make_pair(out0, bn0));
		assert(ret.second); //output loop shouldn't be in map already!
	};

	auto make_end = [&loop_to_bn](char bed0, int32_t needle0, Loop const& in0) {
		BedNeedle bn0(bed0 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle0);

		{ //remove input loop:
			auto f = loop_to_bn.find(in0);
			assert(f != loop_to_bn.end());
			assert(f->second == bn0); //input to stitch should be on this needle
			loop_to_bn.erase(f);
		}

		for (auto const& lbn : loop_to_bn) {
			assert(lbn.second != bn0); //needle should be vacant now that input is removed (NOTE: inefficient check!)
		}
	};

	auto make_xfers = [&loop_to_bn](std::vector< BedNeedle > const& from, std::vector< BedNeedle > const& to, Storage const& loops) {
		assert(from.size() == loops.size());
		assert(to.size() == loops.size());
		for (uint32_t li = 0; li < loops.size(); ++li) {
			//std::cout << "xfer " << loops[li].to_string() << " at " << typeset_bed_needle(from[li].bed, from[li].needle) << " to " << typeset_bed_needle(to[li].bed, to[li].needle) << std::endl; //DEBUG

			auto f = loop_to_bn.find(loops[li]);
			assert(f != loop_to_bn.end());
			/*if (f->second != from[li]) {
				std::cout << "  is at " << typeset_bed_needle(f->second.bed, f->second.needle) << " instead?!?" << std::endl;
			}*/
			assert(f->second == from[li]);
			loop_to_bn.erase(f);
		}
		//TODO: check that range used during xfer planning is actually empty
		for (uint32_t li = 0; li < loops.size(); ++li) {
			auto ret = loop_to_bn.insert(std::make_pair(loops[li], to[li]));
			assert(ret.second);
		}
	};


	//helper: flop to one bed (returns bed flopped to, if no preference given)
	auto stash = [&add_instr, &loop_to_bn, &typeset_bed_needle, &make_xfers, this](Storage const& storage, char const bed = '\0', int32_t offset = 0) -> char { //n.b. nonzero offset only makes sense when already stashed on *other* bed (maybe?)

		if (offset != 0) assert(bed != '\0');

		char target;
		{
			char stashed_to = '\0';
			for (auto const& loop : storage) {
				auto f = loop_to_bn.find(loop);
				assert(f != loop_to_bn.end());
				if (f->second.bed == BedNeedle::FrontSliders) {
					if (stashed_to == '\0') stashed_to = 'f';
					assert(stashed_to == 'f');
				}
				if (f->second.bed == BedNeedle::BackSliders) {
					if (stashed_to == '\0') stashed_to = 'b';
					assert(stashed_to == 'b');
				}
			}
			if (bed == '\0') {
				assert(offset == 0);
				target = (stashed_to != '\0' ? stashed_to : 'b');
			}
			else {
				if (offset != 0) {
					assert(stashed_to != '\0' && stashed_to != bed);
				}
				target = bed;
			}
		}
		assert(target == 'f' || target == 'b');

		std::string from_str = "";
		std::string to_str = "";
		std::vector< BedNeedle > from_vec, to_vec;
		Storage loops;

		QStringList fromList;
		QStringList toList;

		for (auto const& loop : storage) {
			auto f = loop_to_bn.find(loop);
			assert(f != loop_to_bn.end());

			BedNeedle from = f->second;
			BedNeedle to = f->second;

			if (target == 'b' && from.bed == BedNeedle::FrontSliders) to.bed = BedNeedle::Back;
			if (target == 'b' && from.bed == BedNeedle::Front) to.bed = BedNeedle::BackSliders;
			if (target == 'f' && from.bed == BedNeedle::BackSliders) to.bed = BedNeedle::Front;
			if (target == 'f' && from.bed == BedNeedle::Back) to.bed = BedNeedle::FrontSliders;

			to.needle += offset;

			if (from != to) {
				from_vec.emplace_back(from);
				to_vec.emplace_back(to);
				loops.emplace_back(loop);

				if (from_str != "") from_str += ", ";
				from_str += typeset_bed_needle(from.bed, from.needle);
				fromList.append(QString::fromStdString(typeset_bed_needle(from.bed, from.needle)));
				if (to_str != "") to_str += ", ";
				to_str += typeset_bed_needle(to.bed, to.needle);
				toList.append(QString::fromStdString(typeset_bed_needle(to.bed, to.needle)));
			}
		}
		assert(from_vec.size() == to_vec.size());

		if (!from_vec.empty()) {
			std::string instr =
				"stash([" + from_str + "],\n"
				"        [" + to_str + "]);";

			add_instr(instr);
			this->kStash(fromList, toList);

			make_xfers(from_vec, to_vec, loops);
		}
		return target;
	};

	//helper: flop from one bed
	auto unstash = [this, &add_instr, &loop_to_bn, &typeset_bed_needle, &make_xfers](Storage const& storage) {
		std::string from_str = "";
		std::string to_str = "";
		std::vector< BedNeedle > from, to;
		Storage loops;

		QStringList fromList;
		QStringList toList;

		for (auto const& loop : storage) {
			auto f = loop_to_bn.find(loop);
			assert(f != loop_to_bn.end());
			if (f->second.bed == BedNeedle::FrontSliders || f->second.bed == BedNeedle::BackSliders) {
				from.emplace_back(f->second);
				to.emplace_back((f->second.bed == BedNeedle::FrontSliders ? BedNeedle::Back : BedNeedle::Front), f->second.needle);
				loops.emplace_back(loop);
				if (from_str != "") from_str += ", ";
				from_str += typeset_bed_needle(from.back().bed, from.back().needle);
				fromList.append(QString::fromStdString(typeset_bed_needle(from.back().bed, from.back().needle)));
				if (to_str != "") to_str += ", ";
				to_str += typeset_bed_needle(to.back().bed, to.back().needle);
				toList.append(QString::fromStdString(typeset_bed_needle(to.back().bed, to.back().needle)));
			}
		}
		assert(from.size() == to.size());

		if (!from.empty()) {
			std::string instr =
				"unstash([" + from_str + "],\n"
				"          [" + to_str + "]);";

			add_instr(instr);
			this->kUnstash(fromList, toList);
			make_xfers(from, to, loops);
		}
	};
	auto unstash_rec = [&unstash](Storage const& storage) {
		//TODO: should unstash any other storages that overlap this one
		unstash(storage);
	};

	int max_racking = maxRacking; //TODO: all of these auto funcs will be their own func 
	//helper:
	auto reshape = [&max_racking, &add_instr, &stash, &unstash, &typeset_bed_needle, &make_xfers](std::unordered_map< Storage const*, std::pair< int32_t, Shape > > from_layout,
		std::unordered_map< Storage const*, std::pair< int32_t, Shape > > const& to_layout) {
			//move things between different shapes:
			assert(from_layout.size() == to_layout.size()); //layouts should have the same elements

			//sort storages left-to-right:
			std::vector< Storage const* > active;
			for (auto const& sls : from_layout) {
				active.emplace_back(sls.first);
			}
			std::sort(active.begin(), active.end(), [&](Storage const* a, Storage const* b) {
				auto fa = from_layout.find(a);
				assert(fa != from_layout.end());
				auto fb = from_layout.find(b);
				assert(fb != from_layout.end());
				return fa->second.first < fb->second.first;
				});

			//check for anything that needs rolling:
			for (auto storage : active) {
				Shape* from_shape;
				int32_t* from_left;
				{
					auto f = from_layout.find(storage);
					assert(f != from_layout.end());
					from_left = &f->second.first;
					from_shape = &f->second.second;
				}
				Shape const* to_shape;
				int32_t const* to_left;
				{
					auto f = to_layout.find(storage);
					assert(f != to_layout.end());
					to_left = &f->second.first;
					to_shape = &f->second.second;
				}
				if (from_shape->pack() != to_shape->pack()) {
					//must roll the shape.
					{
						uint32_t storage_idx = -1U;
						for (uint32_t i = 0; i < active.size(); ++i) {
							if (active[i] == storage) {
								assert(storage_idx == -1U);
								storage_idx = i;
							}
						}
						assert(storage_idx != -1U);

						int32_t front_min, front_max, back_min, back_max;
						from_shape->size_to_range(storage->size(), &front_min, &front_max, &back_min, &back_max, *from_left);

						//TODO: this would be improved a whole lot by having a stash function that knows about the order of storages, then it could do this recursively, which would be much cleaner!

						{ //stash stuff to the left:
							char prev_bed = '\0';
							int32_t prev_back_min = back_min;
							int32_t prev_front_min = front_min;
							for (uint32_t i = storage_idx - 1; i < active.size(); --i) {
								Storage const* other = active[i];
								assert(other != storage);
								Shape const* other_shape;
								int32_t const* other_left;
								{
									auto f = from_layout.find(other);
									assert(f != from_layout.end());
									other_left = &f->second.first;
									other_shape = &f->second.second;
								}
								int32_t other_front_min, other_front_max, other_back_min, other_back_max;
								other_shape->size_to_range(other->size(), &other_front_min, &other_front_max, &other_back_min, &other_back_max, *other_left);
								assert(other_front_max < prev_front_min);
								assert(other_back_max < prev_back_min);

								char stash_bed;
								if (other_front_max == prev_back_min) {
									// o p p
									// o o p
									stash_bed = (prev_bed == 'b' ? '\0' : 'f');
								}
								else if (other_back_max == prev_front_min) {
									// o o p
									// o p p
									stash_bed = (prev_bed == 'f' ? '\0' : 'b');
								}
								else {
									stash_bed = '\0';
								}

								prev_bed = stash(*other, stash_bed);
								prev_front_min = other_front_min;
								prev_back_min = other_back_min;
							}
						}

						{ //stash stuff to the right:
							char prev_bed = '\0';
							int32_t prev_back_max = back_max;
							int32_t prev_front_max = front_max;
							for (uint32_t i = storage_idx - 1; i < active.size(); --i) {
								Storage const* other = active[i];
								assert(other != storage);
								Shape const* other_shape;
								int32_t const* other_left;
								{
									auto f = from_layout.find(other);
									assert(f != from_layout.end());
									other_left = &f->second.first;
									other_shape = &f->second.second;
								}
								int32_t other_front_min, other_front_max, other_back_min, other_back_max;
								other_shape->size_to_range(other->size(), &other_front_min, &other_front_max, &other_back_min, &other_back_max, *other_left);
								assert(prev_front_max < other_front_min);
								assert(prev_back_max < other_back_min);

								char stash_bed;
								if (prev_front_max == other_back_min) {
									// p o o
									// p p o
									stash_bed = (prev_bed == 'b' ? '\0' : 'f');
								}
								else if (prev_back_max == other_front_min) {
									// p p o
									// p o o
									stash_bed = (prev_bed == 'f' ? '\0' : 'b');
								}
								else {
									stash_bed = '\0';
								}

								prev_bed = stash(*other, stash_bed);
								prev_front_max = other_front_max;
								prev_back_max = other_back_max;
							}
						}
					}

					unstash(*storage);

					std::string from_str = "";
					std::string to_str = "";
					std::vector< BedNeedle > from;
					std::vector< BedNeedle > to;

					for (uint32_t i = 0; i < storage->size(); ++i) {
						char from_bed;
						int32_t from_needle;
						from_shape->size_index_to_bed_needle(storage->size(), i, &from_bed, &from_needle);
						from_needle += *from_left;

						from.emplace_back(from_bed == 'f' ? BedNeedle::Front : BedNeedle::Back, from_needle);

						if (i != 0) from_str += ", ";
						from_str += typeset_bed_needle(from_bed, from_needle);

						char to_bed;
						int32_t to_needle;
						to_shape->size_index_to_bed_needle(storage->size(), i, &to_bed, &to_needle);
						to_needle += *to_left;

						to.emplace_back(to_bed == 'f' ? BedNeedle::Front : BedNeedle::Back, to_needle);

						if (i != 0) to_str += ", ";
						to_str += typeset_bed_needle(to_bed, to_needle);
					}

					add_instr(
						"roll_cycle([" + from_str + "]\n"
						"             [" + to_str + "]);"
					);

					make_xfers(from, to, *storage);

					*from_shape = *to_shape;
				}
			}

			bool need_shift = false;
			for (auto storage : active) {
				Shape* from_shape;
				int32_t* from_left;
				{
					auto f = from_layout.find(storage);
					assert(f != from_layout.end());
					from_left = &f->second.first;
					from_shape = &f->second.second;
				}
				Shape const* to_shape;
				int32_t const* to_left;
				{
					auto f = to_layout.find(storage);
					assert(f != to_layout.end());
					to_left = &f->second.first;
					to_shape = &f->second.second;
				}
				assert(from_shape->pack() == to_shape->pack());

				if (*from_left != *to_left) {
					//there is a need to shift!
					need_shift = true;
				}
			}

			char from_bed = 'b';
			char to_bed = 'f';

			if (need_shift) {
				for (auto storage : active) {
					stash(*storage, from_bed);
				}
			}

			while (need_shift) {
				need_shift = false;

				//TODO: need to do this in the proper order to avoid loops getting snagged at irregular edges:
				for (auto storage : active) {
					Shape* from_shape;
					int32_t* from_left;
					{
						auto f = from_layout.find(storage);
						assert(f != from_layout.end());
						from_left = &f->second.first;
						from_shape = &f->second.second;
					}
					Shape const* to_shape;
					int32_t const* to_left;
					{
						auto f = to_layout.find(storage);
						assert(f != to_layout.end());
						to_left = &f->second.first;
						to_shape = &f->second.second;
					}
					assert(from_shape->pack() == to_shape->pack());

					int32_t offset = *to_left - *from_left;
					//racking limit (attempt!)

					if (offset < -(int32_t)max_racking) offset = -(int32_t)max_racking;
					if (offset > (int32_t)max_racking) offset = (int32_t)max_racking;

					stash(*storage, to_bed, offset);

					*from_left += offset;

					if (*from_left != *to_left) {
						//there is (still) a need to shift!
						need_shift = true;
					}
				}

				std::swap(from_bed, to_bed);
			}

	};

	std::unordered_map< Storage const*, std::pair< int32_t, Shape > > storage_layouts; //<-- is this really needed?
	for (uint32_t stepi = 0; stepi < steps.size(); ++stepi) {

		auto check_storage_layout = [&loop_to_bn](Storage const& storage, int32_t left, Shape const& shape) {
			bool front_stashed = false;
			bool front_unstashed = false;
			bool back_stashed = false;
			bool back_unstashed = false;

			for (uint32_t li = 0; li < storage.size(); ++li) {
				Loop const& loop = storage[li];
				char bed;
				int32_t needle;
				shape.size_index_to_bed_needle(storage.size(), li, &bed, &needle);
				needle += left;

				//std::cout << "Expecting " << loop.to_string() << " at " << typeset_bed_needle(bed, needle) << std::endl; //DEBUG

				auto fl = loop_to_bn.find(loop);
				assert(fl != loop_to_bn.end());
				assert(fl->second.needle == needle);
				if (bed == 'f') {
					assert(fl->second.bed == BedNeedle::Front || fl->second.bed == BedNeedle::BackSliders);
					if (fl->second.bed == BedNeedle::Front) front_unstashed = true;
					if (fl->second.bed == BedNeedle::BackSliders) front_stashed = true;
				}
				else {
					assert(bed == 'b');
					assert(fl->second.bed == BedNeedle::Back || fl->second.bed == BedNeedle::FrontSliders);
					if (fl->second.bed == BedNeedle::Back) back_unstashed = true;
					if (fl->second.bed == BedNeedle::FrontSliders) back_stashed = true;
				}
			}
			assert(!(front_stashed && front_unstashed));
			assert(!(back_stashed && back_unstashed));
			assert(!(front_stashed && back_stashed));
		};
		auto check_storage_layouts = [&step_storages, &storages, &storage_shapes, &storage_layouts, &check_storage_layout](uint32_t stepi) {
			//Check that locations in loop-tracking arrays are as expected:
			for (auto const& sl : step_storages[stepi]) {
				Storage const& storage = storages[sl.storage];
				auto f = storage_layouts.find(&storage);
				assert(f != storage_layouts.end());

				Shape shape = f->second.second;
				int32_t left = f->second.first;

				assert(shape.pack() == storage_shapes[sl.storage]);
				assert(left == sl.left);

				check_storage_layout(storage, left, shape);
			}
		};

		//qDebug() << "Step[" << stepi << "]:"; //DEBUG
		add_instr("//steps[" + std::to_string(stepi) + "]:");
		auto const& step = steps[stepi];
		assert(step_options[stepi] < step.options.size());
		auto const& option = step.options[step_options[stepi]];

		//partition storages into stuff to the left of the step and stuff to the right (useful if we need to shove):
		std::vector< uint32_t > left_of_step;
		std::vector< uint32_t > right_of_step;
		{
			//use 'storage_positions' which are abstract indices that yield a total order:
			int32_t step_min = std::numeric_limits< int32_t >::max();
			int32_t step_max = std::numeric_limits< int32_t >::min();
			for (auto si : step.in) {
				step_min = std::min(step_min, storage_positions[si]);
				step_max = std::max(step_max, storage_positions[si]);
			}
			for (auto si : step.out) {
				step_min = std::min(step_min, storage_positions[si]);
				step_max = std::max(step_max, storage_positions[si]);
			}
			assert(step_min <= step_max);

			std::set< uint32_t > in_step(step.in.begin(), step.in.end());
			for (auto const& sls : storage_layouts) {
				uint32_t idx = sls.first - &storages[0];
				assert(idx < storages.size()); //at this point, no inters should be in storages!
				if (in_step.count(idx)) continue;
				int32_t pos = storage_positions[idx];
				if (pos < step_min) {
					left_of_step.emplace_back(idx);
				}
				else if (pos > step_max) {
					right_of_step.emplace_back(idx);
				}
				else {
					assert(pos < step_min || pos > step_max); //can't have another storage sitting in the middle of the step.
				}
			}
		}
		//DEBUG:
		/*
		qDebug() << "  layout: ";
		for (auto si : left_of_step) std::cout << si << ' ';
		std::cout << '[';
		for (auto i : option.in_order) std::cout << (i == option.in_order[0] ? "" : " ") << step.in[i];
		std::cout << " -> ";
		for (auto o : option.out_order) std::cout << (o == option.out_order[0] ? "" : " ") << step.out[o];
		std::cout << ']';
		for (auto si : right_of_step) std::cout << ' ' << si;
		std::cout << std::endl;
		*/

		//Compute the location of each loop in the step inputs:
		int32_t step_left = 0;
		std::vector< Loop > step_back;
		std::vector< Loop > step_front;

		//helper: will be useful for looking up shape locations in composite shape later:
		auto left_from_step = [&step_left, &step_back, &step_front](Shape const& shape, Storage const& loops) {
			//recomputing this is wasteful but keeps the other code cleaner:
			std::map< Loop, int32_t > step_needle;
			for (uint32_t i = 0; i < step_back.size(); ++i) {
				if (step_back[i] != INVALID_LOOP) {
					auto ret = step_needle.insert(std::make_pair(step_back[i], step_left + i));
					assert(ret.second);
				}
			}
			for (uint32_t i = 0; i < step_front.size(); ++i) {
				if (step_front[i] != INVALID_LOOP) {
					auto ret = step_needle.insert(std::make_pair(step_front[i], step_left + i));
					assert(ret.second);
				}
			}

			//use step_needle as lookup structure to figure out left:
			std::vector< Loop > front, back;
			shape.append_to_beds(loops, INVALID_LOOP, &front, &back);
			int32_t left = std::numeric_limits< int32_t >::max();
			for (uint32_t i = 0; i < back.size(); ++i) {
				if (back[i] == INVALID_LOOP) continue;
				auto l = step_needle.find(back[i]);
				assert(l != step_needle.end());
				int32_t pos = l->second - int32_t(i);
				if (left == std::numeric_limits< int32_t >::max()) left = pos;
				assert(left == pos);
			}
			for (uint32_t i = 0; i < front.size(); ++i) {
				if (front[i] == INVALID_LOOP) continue;
				auto l = step_needle.find(front[i]);
				assert(l != step_needle.end());
				int32_t pos = l->second - int32_t(i);
				if (left == std::numeric_limits< int32_t >::max()) left = pos;
				assert(left == pos);
			}
			//	return step_left + left;
			//  step needle already included step left, don't double count
			return left;
		};

		{ //reshape to get ready for step -- shift everything to the correct offset and (maybe) also transform some things:
			auto old_layouts = storage_layouts;

			//move everything to the in shape: (most things shouldn't need it)
			assert(option.in_shapes.size() == step.in.size());
			for (uint32_t i = 0; i < step.in.size(); ++i) {
				auto f = storage_layouts.find(&storages[step.in[i]]);
				assert(f != storage_layouts.end());
				if (f->second.second.pack() != option.in_shapes[i]) {
					std::cout << "  NOTE: late reshape for storage " << step.in[i] << std::endl; //DEBUG
					f->second.second = Shape::unpack(option.in_shapes[i]); //<-- in some few cases, option's selected layout will be different from storage layout.
				}
			}

			//compute the composite shape:
			for (uint32_t i : option.in_order) {
				auto f = storage_layouts.find(&storages[step.in[i]]);
				assert(f != storage_layouts.end());
				f->second.second.append_to_beds(storages[step.in[i]], INVALID_LOOP, &step_front, &step_back);
			}

			{ //DEBUG:
				std::string front_string, back_string;
				typeset_beds< Loop >(step_front, step_back, [](Loop const& l) {
					return l.to_string();
					}, " ", &front_string, &back_string);
				std::cout << "  in: " << back_string << std::endl;
				std::cout << "      " << front_string << std::endl;
			}


			//figure out the left edge:
			//TODO: some sort of minimization to shift to the middle of the open space.
			if (!step.in.empty()) {
				auto f = storage_layouts.find(&storages[step.in[option.in_order[0]]]);
				assert(f != storage_layouts.end());
				step_left = f->second.first;
			}

			//reposition ins based on the chosen left edge:
			for (auto in : step.in) {
				auto f = storage_layouts.find(&storages[in]);
				assert(f != storage_layouts.end());
				int32_t new_left = left_from_step(f->second.second, storages[in]);
				assert(new_left != std::numeric_limits< int32_t >::max());
				f->second.first = new_left;
			}

			reshape(old_layouts, storage_layouts);
		}


		{ //reinterpret in[0..] as inter + outs[1...]:
			std::vector< Loop > inter_back;
			std::vector< Loop > inter_front;

			if (step.out.empty()) {
				//special case: no outputs, so inter is *just* inter:
				Shape shape = Shape::unpack(option.inter_shape);
				shape.append_to_beds(step.inter, INVALID_LOOP, &inter_front, &inter_back);
			}
			else {
				for (auto o : option.out_order) {
					if (o == 0) {
						if (!step.inter.empty()) {
							Shape shape = Shape::unpack(option.inter_shape);
							shape.append_to_beds(step.inter, INVALID_LOOP, &inter_front, &inter_back);
						}
					}
					else {
						Shape shape = Shape::unpack(option.out_shapes[o]);
						shape.append_to_beds(storages[step.out[o]], INVALID_LOOP, &inter_front, &inter_back);
					}
				}
			}

			{ //DEBUG:
				std::string front_string, back_string;
				typeset_beds< Loop >(inter_front, inter_back, [](Loop const& l) {
					return l.to_string();
					}, " ", &front_string, &back_string);
				std::cout << " int: " << back_string << std::endl;
				std::cout << "      " << front_string << std::endl;
			}

			assert(inter_back == step_back);
			assert(inter_front == step_front);

			//NOTE: if *not* unstashed, reinterpretation could lead to a mixed-stash cycle.
			for (auto in : step.in) {
				unstash(storages[in]);
			}

			//surgery on layouts to reflect the reinterpretation:
			for (uint32_t i = 0; i < step.in.size(); ++i) {
				auto f = storage_layouts.find(&storages[step.in[i]]);
				assert(f != storage_layouts.end());
				storage_layouts.erase(f);
			}

			if (!step.inter.empty()) {
				int32_t left = left_from_step(Shape::unpack(option.inter_shape), step.inter);
				auto ret = storage_layouts.insert(std::make_pair(&step.inter, std::make_pair(left, Shape::unpack(option.inter_shape))));
				assert(ret.second);
			}
			for (uint32_t o = 1; o < step.out.size(); ++o) {
				int32_t left = left_from_step(Shape::unpack(option.out_shapes[o]), storages[step.out[o]]);
				auto ret = storage_layouts.insert(std::make_pair(&storages[step.out[o]], std::make_pair(left, Shape::unpack(option.out_shapes[o]))));
				assert(ret.second);
			}

		}

		//actually get around to doing the step:
		if (stitches[step.begin].type == KnitOutStitch::Start) {
			//cast-ons should be all cast-on:
			for (uint32_t s = step.begin; s != step.end; ++s) {
				assert(stitches[s].type == KnitOutStitch::Start);
			}

			assert(step.inter.empty());
			assert(step.out.size() == 1);
			assert(left_of_step.size() + right_of_step.size() == storage_layouts.size()); //don't have layout for out[0] stored yet.

			{ //move everything else left/right to post-step locations:
				auto old_layouts = storage_layouts;
				for (uint32_t s : left_of_step) {
					auto f = storage_layouts.find(&storages[s]);
					assert(f != storage_layouts.end());
					bool found = false;
					for (auto const& sl : step_storages[stepi]) {
						if (sl.storage == s) {
							assert(!found);
							found = true;
							f->second.first = sl.left;
						}
					}
				}
				for (uint32_t s : right_of_step) {
					auto f = storage_layouts.find(&storages[s]);
					assert(f != storage_layouts.end());
					bool found = false;
					for (auto const& sl : step_storages[stepi]) {
						if (sl.storage == s) {
							assert(!found);
							found = true;
							f->second.first = sl.left;
						}
					}
				}
				reshape(old_layouts, storage_layouts);
			}

			//actually knit a cast-on:
			{ //PARANOIA: does shape-to-needle stuff actually work?
				assert(!storages[step.out[0]].empty());
				uint32_t size = storages[step.out[0]].size();
				std::vector< Loop > out_back;
				std::vector< Loop > out_front;
				Shape shape = Shape::unpack(option.out_shapes[0]);
				shape.append_to_beds(storages[step.out[0]], INVALID_LOOP, &out_front, &out_back);
				for (uint32_t i = 0; i < size; ++i) {
					char bed;
					int32_t needle;
					shape.size_index_to_bed_needle(size, i, &bed, &needle);
					if (bed == 'f') {
						assert(uint32_t(needle) < out_front.size());
						assert(out_front[needle] == storages[step.out[0]][i]);
					}
					else if (bed == 'b') {
						assert(uint32_t(needle) < out_back.size());
						assert(out_back[needle] == storages[step.out[0]][i]);
					}
					else {
						assert(bed == 'f' || bed == 'b');
					}
				}
			}

			{
				Shape shape = Shape::unpack(option.out_shapes[0]);
				int32_t left = std::numeric_limits< int32_t >::max();
				for (auto const& sl : step_storages[stepi]) {
					if (sl.storage == step.out[0]) {
						assert(left == std::numeric_limits< int32_t >::max());
						left = sl.left;
					}
				}
				assert(left != std::numeric_limits< int32_t >::max());

				std::map< Loop, uint32_t > loop_inds;
				for (uint32_t o = 0; o < storages[step.out[0]].size(); ++o) {
					auto ret = loop_inds.insert(std::make_pair(storages[step.out[0]][o], o));
					assert(ret.second);
				}

				//TODO: need to unstash left/right if there are ragged edge overlaps because cast-on would otherwise end up on the wrong side of adjacent tubes.
				
				QString dir = stitches[step.begin].direction == KnitOutStitch::Clockwise ? "clockwise" : "anticlockwise";
				QStringList bns;

				std::string instr = "start_tube(";
				instr += (stitches[step.begin].direction == KnitOutStitch::Clockwise ? "'clockwise'" : "'anticlockwise'");
				instr += ", [";
				for (uint32_t s = step.begin; s != step.end; ++s) {
					auto const& st = stitches[s];
					assert(st.type == KnitOutStitch::Start);
					assert(st.out[0] != -1U && st.out[1] == -1U);

					Loop out0(s, 0);
					auto f0 = loop_inds.find(out0);
					assert(f0 != loop_inds.end());

					{//PARANOIA: should be exactly one loop here:
						uint32_t on_needle = 0;
						for (auto const& li : loop_inds) {
							if (li.second == f0->second) ++on_needle;
						}
						assert(on_needle == 1);
					}

					char bed;
					int32_t needle;
					shape.size_index_to_bed_needle(storages[step.out[0]].size(), f0->second, &bed, &needle);
					needle += left;

					if (s != step.begin) instr += ", ";
					instr += typeset_bed_needle(bed, needle);

					bns.append(QString::fromStdString(typeset_bed_needle(bed, needle)));

					make_start(bed, needle, Loop(s, 0));

					stitch_locations[&st - &stitches[0]].first = BedNeedle(bed == 'f' ? BedNeedle::Front : BedNeedle::Back, needle);
					input_locations[&st - &stitches[0]].first = BedNeedle(bed == 'f' ? BedNeedle::Front : BedNeedle::Back, needle);

					;
				}
				instr += "]);";

				add_instr(instr);
				startTube(dir, bns);

				//add new tube to layouts:
				auto ret = storage_layouts.insert(std::make_pair(&storages[step.out[0]], std::make_pair(left, shape)));
				assert(ret.second);
			}

		}
		else if (stitches[step.begin].type == KnitOutStitch::End) {
			assert(step.out.empty());
			assert(step.in.size() == 1);

			Shape inter_shape = Shape::unpack(option.inter_shape);
			int32_t inter_left = left_from_step(inter_shape, step.inter);

			std::map< Loop, uint32_t > loop_inds;
			for (uint32_t i = 0; i < step.inter.size(); ++i) {
				auto ret = loop_inds.insert(std::make_pair(step.inter[i], i));
				assert(ret.second);
			}

			QString dir = (stitches[step.begin].direction == KnitOutStitch::Clockwise ? "clockwise" : "anticlockwise");
			QStringList bns;

			std::string instr = "end_tube(";
			instr += (stitches[step.begin].direction == KnitOutStitch::Clockwise ? "'clockwise'" : "'anticlockwise'");
			instr += ", [";
			for (uint32_t s = step.begin; s != step.end; ++s) {
				auto const& st = stitches[s];
				assert(st.type == KnitOutStitch::End);
				assert(st.in[0] != -1U && st.in[1] == -1U);

				Loop in0(st.in[0], stitches[st.in[0]].find_out(s));
				assert(in0.idx == 0 || in0.idx == 1);
				auto f0 = loop_inds.find(in0);
				assert(f0 != loop_inds.end());

				{//PARANOIA: should be exactly one loop here:
					uint32_t on_needle = 0;
					for (auto const& li : loop_inds) {
						if (li.second == f0->second) ++on_needle;
					}
					assert(on_needle == 1);
				}

				char bed;
				int32_t needle;
				inter_shape.size_index_to_bed_needle(step.inter.size(), f0->second, &bed, &needle);
				needle += inter_left;

				if (s != step.begin) instr += ", ";
				instr += typeset_bed_needle(bed, needle);
				bns.append(QString::fromStdString(typeset_bed_needle(bed, needle)));

				assert(loop_to_bn.find(in0) != loop_to_bn.end());
				input_locations[&st - &stitches[0]].first = loop_to_bn.find(in0)->second;
				make_end(bed, needle, in0);
				stitch_locations[&st - &stitches[0]].first = BedNeedle(bed == 'f' ? BedNeedle::Front : BedNeedle::Back, needle);
			}
			instr += "]);";

			add_instr(instr);
			endTube(dir, bns);

			{ //no longer need to track inter layout; will replace with out[0]:
				auto f = storage_layouts.find(&step.inter);
				assert(f != storage_layouts.end());
				storage_layouts.erase(f);
			}

		}
		else {
			assert(step.out.size() >= 1);
			assert(step.in.size() >= 1);

			Shape inter_shape = Shape::unpack(option.inter_shape);
			int32_t inter_left = left_from_step(inter_shape, step.inter);

			Shape out_shape = Shape::unpack(option.out_shapes[0]);
			int32_t out_left = std::numeric_limits< int32_t >::max();
			for (auto const& sl : step_storages[stepi]) {
				if (sl.storage == step.out[0]) {
					assert(out_left == std::numeric_limits< int32_t >::max());
					out_left = sl.left;
				}
			}
			assert(out_left != std::numeric_limits< int32_t >::max());

			//transform from inter_shape / inter_left to out_shape / out_left:
			std::string from_str, to_str;
			std::vector< BedNeedle > from, to;

			QStringList fromList;
			QStringList toList;

			bool need_xfer = false;
			for (uint32_t i = 0; i < step.inter.size(); ++i) {
				char from_bed;
				int32_t from_needle;
				inter_shape.size_index_to_bed_needle(step.inter.size(), i, &from_bed, &from_needle);
				from_needle += inter_left;

				uint32_t to_i = step.inter_to_out[i];
				assert(to_i < storages[step.out[0]].size());
				char to_bed;
				int32_t to_needle;
				out_shape.size_index_to_bed_needle(storages[step.out[0]].size(), to_i, &to_bed, &to_needle);
				to_needle += out_left;

				if (!from_str.empty()) from_str += ", ";
				from_str += typeset_bed_needle(from_bed, from_needle);
				fromList.append(QString::fromStdString(typeset_bed_needle(from_bed, from_needle)));

				if (!to_str.empty()) to_str += ", ";
				to_str += typeset_bed_needle(to_bed, to_needle);
				toList.append(QString::fromStdString(typeset_bed_needle(to_bed, to_needle)));

				from.emplace_back(from_bed == 'f' ? BedNeedle::Front : BedNeedle::Back, from_needle);

				to.emplace_back(to_bed == 'f' ? BedNeedle::Front : BedNeedle::Back, to_needle);

				if (from_bed != to_bed || from_needle != to_needle) need_xfer = true;
			}
			if (need_xfer) {
				//Figure out which needles are used during this xform:
				int32_t used_min = std::numeric_limits< int32_t >::max();
				int32_t used_max = std::numeric_limits< int32_t >::min();
				{
					int32_t front_min, front_max, back_min, back_max;
					inter_shape.size_to_range(step.inter.size(), &front_min, &front_max, &back_min, &back_max, inter_left);
					used_min = std::min(used_min, std::min(front_min, back_min));
					used_max = std::max(used_max, std::max(front_max, back_max));
				}
				{
					int32_t front_min, front_max, back_min, back_max;
					out_shape.size_to_range(storages[step.out[0]].size(), &front_min, &front_max, &back_min, &back_max, out_left);
					used_min = std::min(used_min, std::min(front_min, back_min));
					used_max = std::max(used_max, std::max(front_max, back_max));
				}
				//Figure out which other needles are occupied:
				int32_t left_max = std::numeric_limits< int32_t >::min();
				int32_t right_min = std::numeric_limits< int32_t >::max();
				bool on_right = false;
				for (auto const& sl : step_storages[stepi]) {
					if (sl.storage == step.out[0]) {
						assert(!on_right);
						on_right = true;
					}
					else {
						{ //make sure storage is the shape I expect:
							auto f = storage_layouts.find(&storages[sl.storage]);
							assert(f != storage_layouts.end());
							//assert(f->second.first == sl.left); //not sure this will always be the case
							assert(f->second.second.pack() == storage_shapes[sl.storage]);
						}
						int32_t front_min, front_max, back_min, back_max;
						Shape::unpack(storage_shapes[sl.storage]).size_to_range(storages[sl.storage].size(), &front_min, &front_max, &back_min, &back_max, sl.left);
						if (!on_right) {
							left_max = std::max(left_max, std::max(front_max, back_max));
						}
						else {
							right_min = std::min(right_min, std::min(front_min, back_min));
						}
					}
				}

				int32_t shift_left = 0;
				if (left_max >= used_min) {
					assert(left_max == used_min);
					//need to move left at least one for some space:
					shift_left = -1;
				}
				int32_t shift_right = 0;
				if (used_max >= right_min) {
					assert(used_max == right_min);
					//need to move right at least one for some space:
					shift_right += 1;
				}
				if ((left_max >= used_min - 1) && (used_max + 1 >= right_min)) {
					//don't have a gap of one on at least one side
					if (shift_left != 0) shift_left -= 1;
					else shift_right += 1;
				}

				if (shift_left != 0 || shift_right != 0) {
					auto old_layouts = storage_layouts;

					bool on_right = false;
					for (auto const& sl : step_storages[stepi]) {
						if (sl.storage == step.out[0]) {
							assert(!on_right);
							on_right = true;
						}
						else {
							auto f = storage_layouts.find(&storages[sl.storage]);
							assert(f != storage_layouts.end());

							if (!on_right) {
								f->second.first += shift_left;
							}
							else {
								f->second.first += shift_right;
							}
						}
					}

					reshape(old_layouts, storage_layouts);
				}

				unstash_rec(step.inter);

				//now do xfers:

				Constraints constraints;
				constraints.min_free = std::max(left_max + shift_left, used_min - 20);
				constraints.max_free = std::min(right_min + shift_right, used_max + 20);
				// with half-gauging max racking 4 can cause a real racking of 9
				constraints.max_racking = max_racking;

				std::vector< Slack > slack;
				slack.reserve(from.size());
				for (uint32_t i = 0; i < from.size(); ++i) {
					Slack s = 1;
					s = std::max(s, std::abs(from[i].needle - from[i + 1 < from.size() ? i + 1 : 0].needle));
					s = std::max(s, std::abs(to[i].needle - to[i + 1 < to.size() ? i + 1 : 0].needle));
					slack.emplace_back(s);
				}

				std::string error = "";
				std::vector< Transfer > transfers;

				if (!plan_transfers(constraints, from, to, slack, &transfers, &error)) {
					throw std::runtime_error("Failed to plan cycle transfer: " + error);
				}
				std::string transfers_str;
				QList<QPair<QString, QString>> transfersList;
				for (auto const& t : transfers) {
					if (transfers_str != "") transfers_str += ", ";
					transfers_str += "[" + typeset_bed_needle(t.from.bed, t.from.needle) + "," + typeset_bed_needle(t.to.bed, t.to.needle) + "]";
					transfersList.append(QPair<QString, QString>(QString::fromStdString(typeset_bed_needle(t.from.bed, t.from.needle)), 
						QString::fromStdString(typeset_bed_needle(t.to.bed, t.to.needle))));
				}

				std::string instr =
					"xfer_cycle({minFree:" + std::to_string(constraints.min_free) + ", maxFree:" + std::to_string(constraints.max_free) + ", maxRacking:" + std::to_string(constraints.max_racking) + "},\n"
					"             [" + from_str + "],\n"
					"             [" + to_str + "],\n"
					"             [" + transfers_str + "]);";

				add_instr(instr);

				QString dummy;
				
				xferCycle(dummy, dummy, dummy, transfersList);
				make_xfers(from, to, step.inter);
			}

			//make sure inter is ready to knit:
			unstash_rec(step.inter); //this needs to recursively unstash neighbors

			{ //perform stitches:

				//record where each loop in 'inter' currently resides in the out shape:
				// (doing this because inter is in ccw order, and stitches may be in either direction)
				std::map< Loop, uint32_t > loop_inds;
				for (uint32_t i = 0; i < step.inter.size(); ++i) {
					auto ret = loop_inds.insert(std::make_pair(step.inter[i], step.inter_to_out[i]));
					assert(ret.second);
				}

				for (uint32_t s = step.begin; s != step.end; ++s) {
					auto const& st = stitches[s];
					if (st.type == KnitOutStitch::Knit || st.type == KnitOutStitch::Miss || st.type == KnitOutStitch::Tuck) {
						//should have one input:
						assert(st.in[0] != -1U && st.in[1] == -1U);
						//make sure it's in the expected place in step.inter:
						Loop in0(st.in[0], stitches[st.in[0]].find_out(s));
						assert(in0.idx == 0 || in0.idx == 1);
						auto f0 = loop_inds.find(in0);
						assert(f0 != loop_inds.end());

						{//PARANOIA: should be exactly one loop here:
							uint32_t on_needle = 0;
							for (auto const& li : loop_inds) {
								if (li.second == f0->second) ++on_needle;
							}
							assert(on_needle == 1);
						}

						char bed;
						int32_t needle;
						out_shape.size_index_to_bed_needle(storages[step.out[0]].size(), f0->second, &bed, &needle);
						needle += out_left;

						{
							std::string instr;
							if (st.type == KnitOutStitch::Knit) instr += "knit(";
							else if (st.type == KnitOutStitch::Tuck) instr += "tuck(";
							else if (st.type == KnitOutStitch::Miss) instr += "miss(";
							else assert(0 && "bad stitch type");

							instr += '\'';
							QChar stitchType;
							if (st.direction == KnitOutStitch::Clockwise) {
								instr += (bed == 'f' ? '-' : '+');
								stitchType = (bed == 'f' ? '-' : '+');
							}
							else {
								assert(st.direction == KnitOutStitch::Counterclockwise);
								instr += (bed == 'f' ? '+' : '-');
								stitchType = (bed == 'f' ? '+' : '-');
							}
							instr += '\'';
							instr += ", " + typeset_bed_needle(bed, needle) + ");";

							add_instr(instr);

							
							QString bedNeedle = QString::fromStdString(typeset_bed_needle(bed, needle));

							qDebug() << instr;
							assert(stitchType == '-' || stitchType == '+');
							QString stitchTypeStr = stitchType;
							
							if (st.type == KnitOutStitch::Knit) {
								knit(stitchTypeStr, bedNeedle);
							}
							else if (st.type == KnitOutStitch::Tuck) {
								tuck(stitchTypeStr, bedNeedle);
							}
							else if (st.type == KnitOutStitch::Miss) {
								miss(stitchTypeStr, bedNeedle);
							}


							assert(loop_to_bn.find(in0) != loop_to_bn.end());
							// set this before make_stitch clears the old loop:
							input_locations[&st - &stitches[0]].first = loop_to_bn.find(in0)->second;
							make_stitch(bed, needle, in0, Loop(s, 0));
							stitch_locations[&st - &stitches[0]].first = BedNeedle(bed == 'f' ? BedNeedle::Front : BedNeedle::Back, needle);

						}
					}
					else if (st.type == KnitOutStitch::Increase) {
						//should have one input:
						assert(st.in[0] != -1U && st.in[1] == -1U);
						//make sure it's in the expected place in step.inter:
						Loop in0(st.in[0], stitches[st.in[0]].find_out(s));
						assert(in0.idx == 0 || in0.idx == 1);
						auto f0 = loop_inds.find(in0);
						assert(f0 != loop_inds.end());

						{//PARANOIA: should be exactly one loop here:
							uint32_t on_needle = 0;
							for (auto const& li : loop_inds) {
								if (li.second == f0->second) ++on_needle;
							}
							assert(on_needle == 1);
						}

						char bed0;
						int32_t needle0;
						out_shape.size_index_to_bed_needle(storages[step.out[0]].size(), f0->second, &bed0, &needle0);
						needle0 += out_left;

						uint32_t idx1 = f0->second;
						if (st.direction == KnitOutStitch::Clockwise) {
							idx1 = (idx1 > 0 ? idx1 - 1 : storages[step.out[0]].size() - 1);
						}
						else {
							idx1 = (idx1 + 1 < storages[step.out[0]].size() ? idx1 + 1 : 0);
						}

						{//PARANOIA: shouldn't be a loop here:
							for (auto const& li : loop_inds) {
								assert(li.second != idx1);
							}
						}

						char bed1;
						int32_t needle1;
						out_shape.size_index_to_bed_needle(storages[step.out[0]].size(), idx1, &bed1, &needle1);
						needle1 += out_left;

						QChar instrType1, instrType2;
						QString instrTypeStr1, instrTypeStr2;
						QString instrArg1, instrArg2;
						{
							std::string instr = "increase(";
							instr += '\'';
							if (st.direction == KnitOutStitch::Clockwise) {
								instr += (bed0 == 'f' ? '-' : '+');
								instrType1 = bed0 == 'f' ? '-' : '+';
							}
							else {
								assert(st.direction == KnitOutStitch::Counterclockwise);
								instr += (bed0 == 'f' ? '+' : '-');
								instrType1 = bed0 == 'f' ? '+' : '-';
							}
							instrTypeStr1 = instrType1;

							instr += '\'';
							instr += ", ";
							instr += typeset_bed_needle(bed0, needle0);
							instrArg1 = QString::fromStdString(typeset_bed_needle(bed0,needle0));
							instr += ", ";

							instr += '\'';
							if (st.direction == KnitOutStitch::Clockwise) {
								instr += (bed1 == 'f' ? '-' : '+');
								instrType2 = bed1 == 'f' ? '-' : '+';
							}
							else {
								assert(st.direction == KnitOutStitch::Counterclockwise);
								instr += (bed1 == 'f' ? '+' : '-');
								instrType2 = bed1 == 'f' ? '+' : '-';
							}
							instrTypeStr2 = instrType2;
							
							instr += '\'';
							instr += ", ";
							instr += typeset_bed_needle(bed1, needle1);
							instrArg2 = QString::fromStdString(typeset_bed_needle(bed1, needle1));
							instr += ");";

							add_instr(instr);
							increase(instrTypeStr1, instrArg1, instrTypeStr2, instrArg2);

							assert(loop_to_bn.find(in0) != loop_to_bn.end());
							input_locations[&st - &stitches[0]].first = loop_to_bn.find(in0)->second;
							make_increase(bed0, needle0, bed1, needle1, in0, Loop(s, 0), Loop(s, 1));
							stitch_locations[&st - &stitches[0]].first = BedNeedle(bed0 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle0);
							stitch_locations[&st - &stitches[0]].second = BedNeedle(bed1 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle1);
						}

					}
					else if (st.type == KnitOutStitch::Decrease) {
						//should have two inputs:
						assert(st.in[0] != -1U && st.in[1] != -1U);

						//make sure inputs are in the expected place in step.inter:
						Loop in0(st.in[0], stitches[st.in[0]].find_out(s));
						assert(in0.idx == 0 || in0.idx == 1);
						auto f0 = loop_inds.find(in0);
						assert(f0 != loop_inds.end());

						Loop in1(st.in[1], stitches[st.in[1]].find_out(s));
						assert(in1.idx == 0 || in1.idx == 1);
						auto f1 = loop_inds.find(in1);
						assert(f1 != loop_inds.end());
						assert(f1->second == f0->second);

						{//PARANOIA: should be exactly two loops here:
							uint32_t on_needle = 0;
							for (auto const& li : loop_inds) {
								if (li.second == f0->second) ++on_needle;
							}
							assert(on_needle == 2);
						}

						char bed0;
						int32_t needle0;
						out_shape.size_index_to_bed_needle(storages[step.out[0]].size(), f0->second, &bed0, &needle0);
						needle0 += out_left;

						QChar instrType1, instrType2;
						QString instrTypeStr1, instrTypeStr2;
						QString instrArg1, instrArg2;

						{
							std::string instr = "decrease(";
							instr += '\'';
							if (st.direction == KnitOutStitch::Clockwise) {
								instr += (bed0 == 'f' ? '-' : '+');
								instrType1 = bed0 == 'f' ? '-' : '+';
							}
							else {
								assert(st.direction == KnitOutStitch::Counterclockwise);
								instr += (bed0 == 'f' ? '+' : '-');
								instrType1 = bed0 == 'f' ? '+' : '-';
							}
							instr += '\'';
							instr += ", ";
							instr += typeset_bed_needle(bed0, needle0);
							instrArg1 = QString::fromStdString(typeset_bed_needle(bed0, needle0));
							instr += ");";

							add_instr(instr);
							instrTypeStr1 = instrType1;
							decrease(instrTypeStr1, instrArg1);
							// awkward to use loop_to_bn.find()..,since already xferred...

							input_locations[&st - &stitches[0]].first = stitch_locations[st.in[0]].first;
							input_locations[&st - &stitches[0]].second = stitch_locations[st.in[1]].first;


							make_decrease(bed0, needle0, in0, in1, Loop(s, 0));
							stitch_locations[&st - &stitches[0]].first = BedNeedle(bed0 == 'f' ? BedNeedle::Front : BedNeedle::Back, needle0);

						}

					}
					else {
						assert(0 && "bad stitch type");
					}
				}
			}

			check_storage_layout(storages[step.out[0]], out_left, out_shape); //DEBUG


			{ //no longer need to track inter layout; will replace with out[0]:
				auto f = storage_layouts.find(&step.inter);
				assert(f != storage_layouts.end());
				storage_layouts.erase(f);
			}

			//add new tube to layouts:
			auto ret = storage_layouts.insert(std::make_pair(&storages[step.out[0]], std::make_pair(out_left, out_shape)));
			assert(ret.second);

		}



		{ //make sure everything gets to the final output position:
			auto old_layouts = storage_layouts;

			for (auto const& sl : step_storages[stepi]) {
				auto f = storage_layouts.find(&storages[sl.storage]);
				assert(f != storage_layouts.end());
				assert(f->second.second.pack() == storage_shapes[sl.storage]);
				f->second.first = sl.left;
			}

			reshape(old_layouts, storage_layouts);
		}

		//Check that locations in loop-tracking arrays are as expected:
		check_storage_layouts(stepi);
	}
	add_instr("write();");


	//qDebug() << "done!!!";
	

	scheduledInstructions.clear();
	for (auto instr : instructions) {
		scheduledInstructions.push_back(QString::fromStdString(instr));
	}


	//emit instructionsCreated(instructions);

	//original also saves schedules st file, not needed?
	emit knitoutGenerated(knitout);
	emit helpBoxCommunication("DONE - resulting 'knitout.k' was saved inside the project folder!");
	emit helpBoxCommunication("Total amount of instructions: " + QString::number(knitout.size()));
	//nerateKnitout();
	//emit knitoutGenerated(knitout);

	return;
}


PortedBedNeedle KnitoutScheduler::parseBedNeedle(QString& bedNeedle)
{

	if (bedNeedle.isEmpty() || !bedNeedle.at(bedNeedle.size() - 1).isDigit())
		qDebug() << "parseBedNeedle must be called with a string";

	QRegularExpression pattern("^([fb]s?)([-+]?\\d+)$");
	QRegularExpressionMatch match = pattern.match(bedNeedle);
	if (!match.hasMatch())
		qDebug() << "string '" << bedNeedle.toStdString() << "' does not look like a needle";

	return { match.captured(1), match.captured(2).toInt() };

}

QString KnitoutScheduler::bnToHalf(QString& bn_str) {
	PortedBedNeedle bn = parseBedNeedle(bn_str);
	if (bn.bed == "f") {
		return "f" + QString::number(2 * bn.needle);
	}
	else if (bn.bed == "fs") {
		return "f" + QString::number(2 * bn.needle + 1);
	}
	else if (bn.bed == "b") {
		return "b" + QString::number(2 * bn.needle + 1);
	}
	else if (bn.bed == "bs") {
		return "b" + QString::number(2 * bn.needle);
	}
	else {
		qDebug() << "don't know how to half-guage the needle '" + bn_str.toStdString() + "'";
	}
}

void KnitoutScheduler::bringIn(QString& Carrier) {
	qDebug() << "Bringing in Carrier " << Carrier << " before use.";
	
	QString output = "inhook " + Carrier;
	out(output);

	inYarns[Carrier] = true;
	needsReleaseHook[Carrier] = true;
	yarnUsedCount[Carrier] = 0;
}

void KnitoutScheduler::knit(QString& d, QString& bn) {
	if (!inYarns[Carrier]) {
		qDebug() << "Using carrier before inhook : " << Carrier;
		this->bringIn(Carrier);
	}
	if (needsReleaseHook[Carrier] && yarnUsedCount[Carrier] > 5) {
		QString output = "releasehook " + Carrier;
		out(output);
		needsReleaseHook[Carrier] = false;
	}
	yarnUsedCount[Carrier] += 1;

	assert(!d.isEmpty());

	QString output = "knit " + d + " " + bnToHalf(bn) + " " + Carrier;
	out(output);
}

void KnitoutScheduler::drop(QString& bn) {
	QString output = "drop " + bnToHalf(bn);
	out(output);
}

void KnitoutScheduler::tuck(QString& d, QString& bn) {
	if (!inYarns[Carrier]) {
		qDebug() << "Using carrier before inhook : " << Carrier;
		bringIn(Carrier);
	}
	if (this->needsReleaseHook[Carrier] && this->yarnUsedCount[Carrier] > 5) {
		QString output = "releasehook " + Carrier;
		out(output);
		needsReleaseHook[Carrier] = false;
	}
	yarnUsedCount[Carrier] += 1;
	
	QString output = "tuck " + d + " " + bnToHalf(bn) + " " + Carrier;
	out(output);
}

void KnitoutScheduler::miss(QString& d, QString& bn) {
	if (!inYarns[Carrier]) {
		qDebug() << "Using carrier before inhook : " << Carrier;
		bringIn(Carrier);
	}
	if (this->needsReleaseHook[Carrier] && this->yarnUsedCount[Carrier] > 5) {
		QString output = "releasehook " + Carrier;
		out(output);
		needsReleaseHook[Carrier] = false;
	}

	QString output = "miss " + d + " " + bnToHalf(bn) + " " + Carrier;
	this->out(output);
}

void KnitoutScheduler::decrease(QString& d, QString& bn) {
	qDebug() << "decresing";
	knit(d, bn);
}

void KnitoutScheduler::increase(QString& d0, QString& bn0, QString& d1, QString& bn1) {
	qDebug() << "increaseing";
	knit(d0, bn0);

	QString arg1 = d1 == "+" ? "-" : "+";
	tuck(arg1, bn1);
}

void KnitoutScheduler::xfer(QString& from, QString& to) {
	// If any releasehooks are pending, do them before xfers
	for (const auto& pair : needsReleaseHook) {
		
		if (pair.second) {
			QString output = "releasehook " + pair.first;
			out(output);
			needsReleaseHook[pair.first] = false;
		}
	}
	QString output = "xfer " + bnToHalf(from) + " " + bnToHalf(to);
	out(output);
}

void KnitoutScheduler::setRacking(QString& from_str, QString& to_str) {
	int target;
	if (from_str.isEmpty() && to_str.isEmpty()) {
		target = 0;
	}
	else {

		QString from_half = bnToHalf(from_str);
		QString to_half = bnToHalf(to_str);

		PortedBedNeedle from = parseBedNeedle(from_half);
		PortedBedNeedle to = parseBedNeedle(to_half);
		if (from.bed == "f" && to.bed == "b") {
			target = from.needle - to.needle;
		}
		else {
			assert(from.bed == "b" && to.bed == "f");
			target = to.needle - from.needle;
		}
		assert(qAbs(target) <= 8);
	}
	if (racking != target) {
		racking = target;
		QString output = "rack " + QString::number(racking);
		out(output);
	}
}

void KnitoutScheduler::xferCycle(QString& opts, QString& from, QString& to, QList<QPair<QString, QString>>& xfers) {
	for (const auto& xf : xfers) {
		QString first = xf.first;
		QString second = xf.second;
		setRacking(first, second);
		xfer(first, second);
	}
}

void KnitoutScheduler::kStash(QStringList& from, QStringList& to) {
	if (from.size() != to.size()) qDebug() << "from and to arrays should be the same length";
	if (from.isEmpty()) return;

	setRacking(from.first(), to.first());

	for (int i = 0; i < from.size(); ++i) {
		xfer(from[i], to[i]);
	}

	QString empty = "";
	QString empty2 = "";
	setRacking(empty, empty2);
}

void KnitoutScheduler::kUnstash(QStringList& from, QStringList& to) {
	if (from.size() != to.size()) qDebug() << "from and to arrays should be the same length";
	if (from.isEmpty()) return;

	setRacking(from.first(), to.first());

	for (int i = 0; i < from.size(); ++i) {
		this->xfer(from[i], to[i]);
	}

	QString empty = "";
	QString empty2 = "";
	this->setRacking(empty, empty2);
}

void KnitoutScheduler::endTube(QString& dir, QStringList& bns) {
	auto d = [this, &dir](QString& bn_str) {
		PortedBedNeedle bn = parseBedNeedle(bn_str);
		if (bn.bed.startsWith('f')) {
			if (dir == "clockwise") {
				return "-";
			}
			else {
				Q_ASSERT(dir == "anticlockwise");
				return "+";
			}
		}
		else {
			Q_ASSERT(bn.bed.startsWith('b'));
			if (dir == "clockwise") {
				return "+";
			}
			else {
				Q_ASSERT(dir == "anticlockwise");
				return "-";
			}
		}
	};


	QChar c(CastOnStitchNumber);
	QString output = "x-stitch-number " + QString(QChar(CastOnStitchNumber));
	out(output);
	qDebug() << "1";
	for (int i = 0; i < bns.size(); ++i) {
		QString l = d(bns[i]);
		if (i % 2 == 0) knit(l, bns[i]);
	}
	for (int i = 0; i < bns.size(); ++i) {
		QString l = d(bns[i]);
		if (i % 2 == 1) knit(l, bns[i]);
	}

	qDebug() << "2";
	output = "x-stitch-number " + QString(QChar(PlainStitchNumber));
	out(output);
	for (int row = 0; row < EndRows; ++row) {
		for (QString& bn : bns) {
			QString l = d(bn);
			knit(l, bn);
		}
	}

	if (this->inYarns[Carrier]) {
		output = "outhook " + Carrier;
		out(output);
		inYarns[Carrier] = false;
		yarnUsedCount[Carrier] = 0;
		qDebug() << "Taking Carrier " << Carrier << " out.";
	}
	else {
		qDebug() << "Taking out carrier that was never in : " << Carrier;
	}

	for (QString& bn : bns) {
		drop(bn);
	}
}


void KnitoutScheduler::startTube(QString& dir, QStringList& bns)
{
	QList<int> front;
	QList<int> back;
	for (QString& bn_str : bns) {
		PortedBedNeedle bn = parseBedNeedle(bn_str);
		if (bn.bed == "f") {
			front.append(bn.needle);
		}
		else if (bn.bed == "b") {
			back.append(bn.needle);
		}
		else {
			qDebug() << "start_tube should only be called with 'f' or 'b' needles.";
		}
	}

	std::sort(front.begin(), front.end());
	std::sort(back.begin(), back.end());

	if (front.isEmpty() || back.isEmpty()) {
		qDebug() << "should start a tube with at least a stitch on each bed.";
	}

	// Do a tuck pattern to anchor yarn
	int n = qMax(front.last(), back.last());
	QList<QString> toDrop;
	auto initTuck = [this, &toDrop](QString& d, QString& bn) {
		tuck(d, bn);
		toDrop.append(bn);
	};

	QString output = "x-stitch-number " + QString(QChar(YarnInStitchNumber));
	out(output);


	if (inYarns.find(Carrier) == inYarns.end()) {
		QString output = "inhook " + Carrier;
		out(output);
		inYarns[Carrier] = true;
		yarnUsedCount[Carrier] = 0;
		qDebug() << "Bringing Carrier " << Carrier << " in.";
	}
	else {
		qDebug() << "Bringing in carrier that was already in " << Carrier;
	}

	QString line1 = "f" + QString::number(n);
	QString line2 = "f" + QString::number(n - 2);
	QString line3 = "f" + QString::number(n - 4);
	QString line4 = "f" + QString::number(n - 3);
	QString line5 = "f" + QString::number(n - 1);
	QString minus = "-";
	QString plus = "+";

	initTuck(minus, line1);
	initTuck(minus, line2);
	initTuck(minus, line3);
	initTuck(plus, line4);
	initTuck(plus, line5);


	if (inYarns.find(Carrier) != inYarns.end()) {
		QString line = "releasehook " + Carrier;
		out(line);
	}
		

	// Make list of needles and directions in tube order
	QList<QPair<QString, QString>> sts;
	if (dir == "clockwise") {
		for (int i = front.size() - 1; i >= 0; --i) {
			sts.append({ "-", "f" + QString::number(front[i]) });
		}
		for (int i = 0; i < back.size(); ++i) {
			sts.append({ "+", "b" + QString::number(back[i]) });
		}
	}
	else {
		Q_ASSERT(dir == "anticlockwise");
		for (int i = back.size() - 1; i >= 0; --i) {
			sts.append({ "-", "b" + QString::number(back[i]) });
		}
		for (int i = 0; i < front.size(); ++i) {
			sts.append({ "+", "f" + QString::number(front[i]) });
		}
	}

	// Alternating tuck cast on
	qDebug() << "Alternating";
	output = "x-stitch-number " + QString(QChar(CastOnStitchNumber));
	out(output);
	for (int i = 0; i < sts.size(); ++i) {
		if (i % 2 == 0) knit(sts[i].first, sts[i].second);
	}
	for (int i = 0; i < sts.size(); ++i) {
		if (i % 2 == 1) knit(sts[i].first, sts[i].second);
	}

	// Drop everything in 'toDrop' that wasn't part of alternating tucks
	for (QString& bn : toDrop) {
		if (!bns.contains(bn)) {
			this->drop(bn);
		}
	}

	// Knit some plain rows
	output = "x-stitch-number " + QString(QChar(PlainStitchNumber));
	this->out(output);
	qDebug() << "plain rows";
	for (int row = 0; row < StartRows; ++row) {
		for (auto& dbn : sts) {
			knit(dbn.first, dbn.second);
		}
	}

	int first = 0;
	while (first < sts.size() && !bns.contains(sts[first].second)) ++first;
	Q_ASSERT(first < sts.size());

	// Knit a bit extra to get aligned to the input bns
	qDebug() << "a bit extra";
	for (int i = 0; i < first; ++i) {
		QPair<QString, QString> st = sts.takeFirst();
		this->knit(st.first, st.second);
		sts.append(st);
	}

	// Alternating stitches to separate starting tube from knitting
	output = "x-stitch-number " + QString(QChar(CastOnStitchNumber));
	out(output);
	qDebug() << "awwooge";
	for (int i = 0; i < sts.size(); ++i) {
		if (i % 2 == 0) this->knit(sts[i].first, sts[i].second);
	}
	for (int i = 0; i < sts.size(); ++i) {
		if (i % 2 == 1) this->knit(sts[i].first, sts[i].second);
	}

	output = "x-stitch-number " + QString(QChar(PlainStitchNumber));
	this->out(output);
}






bool KnitoutScheduler::embed_DAG(
	std::vector< DAGNode > const& nodes,
	std::vector< DAGEdge > const& edges,
	std::vector< uint32_t >* node_options,
	std::vector< int32_t >* node_positions, //positions give total left-to-right order of edges/nodes
	std::vector< int32_t >* edge_positions
) {
	for (auto const& node : nodes) {
		if (node.options.empty()) {
			std::cerr << "WARNING: embed_DAG will fail because a node has no options." << std::endl;
			return false;
		}
	}

	{ //PARANOIA:
		for (auto const& node : nodes) {
			assert(!node.options.empty());
			std::set< DAGEdgeIndex > ins(node.options[0].in_order.begin(), node.options[0].in_order.end());
			std::set< DAGEdgeIndex > outs(node.options[0].out_order.begin(), node.options[0].out_order.end());
			//ins and outs reference node properly:
			for (auto i : ins) {
				assert(i < edges.size());
				assert(edges[i].to == &node - &nodes[0]);
			}
			for (auto o : outs) {
				assert(o < edges.size());
				assert(edges[o].from == &node - &nodes[0]);
			}

			for (auto const& option : node.options) {
				assert(option.in_order.size() == ins.size());
				assert(option.in_shapes.size() == ins.size());
				assert(option.out_order.size() == outs.size());
				assert(option.out_shapes.size() == outs.size());
				//node orders all talk about the same edges:
				for (auto i : option.in_order) assert(ins.count(i));
				for (auto o : option.out_order) assert(outs.count(o));
				//node shapes don't index outside edge cost matrices:
				for (uint32_t i = 0; i < option.in_shapes.size(); ++i) {
					assert(option.in_shapes[i] < edges[option.in_order[i]].to_shapes);
				}
				for (uint32_t o = 0; o < option.out_shapes.size(); ++o) {
					assert(option.out_shapes[o] < edges[option.out_order[o]].from_shapes);
				}
			}
		}

		for (auto const& edge : edges) {
			//edges have properly-sized cost matrix:
			assert(edge.costs.size() == edge.from_shapes * edge.to_shapes);
			//edges reference nodes that reference them:
			{
				assert(edge.from < nodes.size());
				DAGNode const& node = nodes[edge.from];
				assert(!node.options.empty());
				auto f = std::find(node.options[0].out_order.begin(), node.options[0].out_order.end(), &edge - &edges[0]);
				assert(f != node.options[0].out_order.end());
			}
			{
				assert(edge.to < nodes.size());
				DAGNode const& node = nodes[edge.to];
				assert(!node.options.empty());
				auto f = std::find(node.options[0].in_order.begin(), node.options[0].in_order.end(), &edge - &edges[0]);
				assert(f != node.options[0].in_order.end());
			}
			//nodes are sorted:
			assert(edge.from < edge.to);
		}
	} //end PARANOIA

	std::vector< uint32_t > select_order;

	//TODO: ~smart ordering~
	for (uint32_t i = 0; i < nodes.size(); ++i) {
		select_order.emplace_back(i);
	}
	std::stable_sort(select_order.begin(), select_order.end(), [&nodes](uint32_t a, uint32_t b) -> bool {
		return nodes[a].options.size() < nodes[b].options.size();
		});

	struct State {
		std::vector< uint32_t > selected; //selected options for each node
		//derived:
		uint32_t step = 0;
		std::vector< bool > left_of; //partial order on the edges
		std::vector< uint32_t > from_shape, to_shape; //selected shapes on edges

		bool operator==(State const& o) const {
			return selected == o.selected;
		}
	};

	struct HashState {
		size_t operator()(State const& state) const {
			static std::hash< std::string > hash;
			return hash(std::string(
				reinterpret_cast<const char*>(&state.selected[0]),
				state.selected.size() * sizeof(uint32_t)));
		}
	};


	std::vector< std::pair< DAGCost, const State* > > to_expand;
	std::unordered_map< State, DAGCost, HashState > visited;
	//std::unordered_set< State, HashState > expanded;

	const uint32_t Initial = 8000000;
	to_expand.reserve(Initial);
	visited.reserve(Initial);
	//expanded.reserve(Initial);

	std::greater< std::pair< DAGCost, const State* > > CompareCost;

	auto queue_state = [&to_expand, &visited, &CompareCost](State const& state, DAGCost const& cost) {
		auto ret = visited.insert(std::make_pair(state, DAGCost::max()));
		if (cost < ret.first->second) {
			ret.first->second = cost;
			to_expand.emplace_back(cost, &(ret.first->first));
			std::push_heap(to_expand.begin(), to_expand.end(), CompareCost);
		}
	};


	//Any edge that crosses a node must be left_of *all* in/out edges of that node:
	std::vector< std::vector< std::vector< DAGEdgeIndex > > > edge_excludes;
	edge_excludes.reserve(edges.size());
	for (auto const& edge : edges) {
		edge_excludes.emplace_back();
		std::vector< std::vector< DAGEdgeIndex > >& excl = edge_excludes.back();
		for (uint32_t n = edge.from + 1; n + 1 < edge.to; ++n) {
			DAGNode const& node = nodes[n];
			excl.emplace_back();
			excl.back().insert(excl.back().begin(), node.options[0].in_order.begin(), node.options[0].in_order.end());
			excl.back().insert(excl.back().begin(), node.options[0].out_order.begin(), node.options[0].out_order.end());
			std::sort(excl.back().begin(), excl.back().end());
		}
	}

	std::function< bool(std::vector< bool >&, uint32_t, uint32_t) > add_left_of;
	auto run_excludes = [&add_left_of, &edges, &edge_excludes](std::vector< bool >& left_of, uint32_t a) -> bool {
		for (std::vector< uint32_t > const& excl : edge_excludes[a]) {
			bool is_left = false;
			bool is_right = false;
			for (uint32_t b : excl) {
				if (left_of[a * edges.size() + b]) is_left = true;
				if (left_of[b * edges.size() + a]) is_right = true;
			}
			if (is_left && is_right) return false;
			if (is_left) {
				for (uint32_t b : excl) {
					if (!add_left_of(left_of, a, b)) return false;
				}
			}
			if (is_right) {
				for (uint32_t b : excl) {
					if (!add_left_of(left_of, b, a)) return false;
				}
			}
		}
		return true;
	};

	add_left_of = [&edges, &add_left_of, &run_excludes](std::vector< bool >& left_of, uint32_t a, uint32_t b) -> bool {
		if (left_of[a * edges.size() + b]) return true;
		if (left_of[b * edges.size() + a]) return false;
		left_of[a * edges.size() + b] = true;
		for (uint32_t c = 0; c < edges.size(); ++c) {
			if (left_of[c * edges.size() + a] && !add_left_of(left_of, c, b)) return false;
			if (left_of[b * edges.size() + c] && !add_left_of(left_of, a, c)) return false;
		}
		if (!run_excludes(left_of, a)) return false;
		if (!run_excludes(left_of, b)) return false;
		return true;
	};

	auto expand_state = [&queue_state, &select_order, &add_left_of, &nodes, &edges](State const& state, DAGCost const& cost) {
		assert(state.step < select_order.size());
		uint32_t n = select_order[state.step];
		DAGNode const& node = nodes[n];
		for (auto const& option : node.options) {
			State next = state;
			next.step += 1;
			next.selected[n] = &option - &node.options[0];

			bool bad = false;
			for (uint32_t i = 1; i < option.in_order.size(); ++i) {
				if (!add_left_of(next.left_of, option.in_order[i - 1], option.in_order[i])) {
					bad = true;
					break;
				}
			}
			if (bad) continue;
			for (uint32_t o = 1; o < option.out_order.size(); ++o) {
				if (!add_left_of(next.left_of, option.out_order[o - 1], option.out_order[o])) {
					bad = true;
					break;
				}
			}
			if (bad) continue;

			DAGCost next_cost = cost;
			next_cost += option.cost;
			for (uint32_t i = 0; i < option.in_order.size(); ++i) {
				uint32_t e = option.in_order[i];
				assert(next.to_shape[e] == -1U);
				next.to_shape[e] = option.in_shapes[i];
				if (next.from_shape[e] != -1U) {
					//add cost from edge if edge got completed:
					auto const& edge_cost = edges[e].costs[next.from_shape[e] * edges[e].to_shapes + next.to_shape[e]];
					if (edge_cost == DAGCost::max()) {
						bad = true;
						break;
					}
					next_cost += edge_cost;
				}
			}
			if (bad) continue;

			for (uint32_t o = 0; o < option.out_order.size(); ++o) {
				uint32_t e = option.out_order[o];
				assert(next.from_shape[e] == -1U);
				next.from_shape[e] = option.out_shapes[o];
				if (next.to_shape[e] != -1U) {
					//add cost from edge if edge got completed:
					auto const& edge_cost = edges[e].costs[next.from_shape[e] * edges[e].to_shapes + next.to_shape[e]];
					if (edge_cost == DAGCost::max()) {
						bad = true;
						break;
					}
					next_cost += edge_cost;
				}
			}
			if (bad) continue;


			queue_state(next, next_cost);

		}
	};
	auto set_output = [&](State const& state) {
		assert(state.step == select_order.size());
		assert(state.selected.size() == nodes.size());
		if (node_options) {
			*node_options = state.selected;
		}
		//come up with a total ordering:
		std::vector< bool > left_of = state.left_of;
		for (uint32_t a = 0; a < edges.size(); ++a) {
			for (uint32_t b = 0; b < edges.size(); ++b) {
				//break ordering ties in an arbitrary way:
				if (a != b && (!left_of[a * edges.size() + b] && !left_of[b * edges.size() + a])) {
					bool ret = add_left_of(left_of, a, b);
					assert(ret);
				}
			}
		}
		//store edge and node positions (easy to compute from total order):
		if (node_positions) {
			node_positions->assign(nodes.size(), std::numeric_limits< int32_t >::max());
			node_positions->reserve(nodes.size());
		}
		if (edge_positions) {
			edge_positions->clear();
			edge_positions->reserve(edges.size());
		}
		for (uint32_t a = 0; a < edges.size(); ++a) {
			int32_t pos = edges.size();
			for (uint32_t b = 0; b < edges.size(); ++b) {
				if (left_of[a * edges.size() + b]) {
					pos -= 1;
				}
			}
			if (edge_positions) {
				edge_positions->emplace_back(pos);
			}
			if (node_positions) {
				if (edges[a].from != -1U) (*node_positions)[edges[a].from] = std::min((*node_positions)[edges[a].from], pos);
				if (edges[a].to != -1U) (*node_positions)[edges[a].to] = std::min((*node_positions)[edges[a].to], pos);
			}
		}

	};

	{
		State start;
		start.selected.resize(nodes.size(), -1U);
		start.step = 0;
		start.left_of.resize(edges.size() * edges.size(), false);
		start.from_shape.resize(edges.size(), -1U);
		start.to_shape.resize(edges.size(), -1U);
		queue_state(start, DAGCost::zero());
	}

	uint32_t step = 0;

	while (!to_expand.empty()) {
		std::pop_heap(to_expand.begin(), to_expand.end(), CompareCost);
		DAGCost cost = to_expand.back().first;
		const State& state = *to_expand.back().second;
		to_expand.pop_back();
		auto f = visited.find(state);
		assert(f != visited.end());
		assert(!(cost < f->second));
		if (cost == f->second) {
			//auto res = expanded.insert(state);
			//assert(res.second);
			if ((++step) % 10000 == 0) {
				//DEBUG:

				QString emitString = "\t(" + QString::number(visited.size()) + "/" + QString::number(to_expand.size())
					+ ") [step " + QString::number(state.step) + "]: ";
				for (auto s : state.selected) emitString += " " + (s == -1U ? "." : QString::number(s));
				emit helpBoxCommunication(emitString);
				//qDebug() << /*expanded.size() << "/" <<*/ visited.size() << "/" << to_expand.size() << "   ";
				//qDebug() << "[" << state.step << "]";
				//for (auto s : state.selected) qDebug() << ' ' << (s == -1U ? std::string(".") : std::to_string(s));
			}

			if (state.step < select_order.size()) {
				expand_state(state, cost);
			}
			else {
				//Found the cheapest selected state!
				set_output(state);

				return true;
			}
		}
	}


	return false;
}
