#pragma once
#include "KnitGrapher.h"



//helper global variable defined in originak ak::trim_model as 'edges'
//in KnitGrapher it is used in trimModel() function
std::vector<KnitGrapher::Edge> KnitGrapher::globedges;

/*
*	Constructor: set default parameters of the KnitGrapher during construction
*		all 'values' are in mm
*/
KnitGrapher::KnitGrapher(QObject* parent) : QObject(parent)
{
	stitchWidth = 3.66f;
	stitchHeight = 1.73f;
	modelUnitLength = 10.0f;
	stepCount = 0;

}

/*
*	Function: set the original mesh that will be used to create the knit graph
*		this mesh will be remeshed and interpolated to create a new mesh that the AutoKnit algorithms will run on
*	Parameters: ObjectMesh mesh - the mesh that will be used to create the knit graph
* 	Return: void
*	Called: called when the user selects a .obj file through the GUI
*/
void KnitGrapher::setOriginalMesh(ObjectMesh mesh)
{	
	AutoKnitMesh step(mesh);
	originalMesh = step;
}


/*
*	Function: get the max edge length allowed in the remeshing algorithm
*		the mathematical function is copied from ak::pipeline.hpp, where its part of the Parameters class
* 	Return: float - the max edge length
* 	Called: called at various stages of the algorithm, mainly remesh() and peelSlice()
*/
float KnitGrapher::getMaxEdgeLength() const
{
	return 0.5f * std::min(stitchWidth, 2.0f * stitchHeight) / modelUnitLength;
}

/* 
*	Function: Dijkstra's algorithm 'visit' function
*	Return: void, but modifies the 'todo' and 'visited' vectors which can be accessed by the calling function
* 	Called: called in the divide() function, during the constraint embedding of remesh()
*/ 
void visit(std::vector< std::pair< float, uint32_t > >& todo,
	std::vector< std::pair< float, uint32_t > >& visited,
	uint32_t vertex, float distance, uint32_t from) 
{
	if (distance < visited[vertex].first) {
		visited[vertex] = std::make_pair(distance, from);
		todo.emplace_back(distance, vertex);
		std::push_heap(todo.begin(), todo.end(), std::greater<std::pair<float, uint32_t>>());
	}
}

/*
*	Function: check if a pair of vertices is in the marked set
*	Return: the second value of the pair if it exists, -1 otherwise
*	Called: called during the divide() function as part of remesh()
*/
uint32_t KnitGrapher::lookup(uint32_t a, uint32_t b, std::unordered_map< glm::uvec2, uint32_t >& marked_verts) const
{
	auto f = marked_verts.find((a < b ? glm::uvec2(a, b) : glm::uvec2(b, a)));
	if (f != marked_verts.end()) return f->second;
	else return -1U;
}

/*
*	Function: split a quad in a mesh into 2 triangles
*	Return: void, but emplaces the newly created triangles into an argument vector
*	Called: called during the divide() function as part of remesh()
*/
void KnitGrapher::quad(std::vector<glm::uvec3>& new_tris, std::vector<glm::vec3> const& verts, 
	uint32_t a, uint32_t b, uint32_t c, uint32_t d) 
{
	float ac = glm::length2(verts[c] - verts[a]);
	float bd = glm::length2(verts[d] - verts[b]);
	if (ac < bd) {
		new_tris.emplace_back(glm::uvec3(a, b, c));
		new_tris.emplace_back(glm::uvec3(c, d, a));
	}
	else {
		new_tris.emplace_back(glm::uvec3(a, b, d));
		new_tris.emplace_back(glm::uvec3(b, c, d));
	}
}

/*
*	Function: divide the mesh triangles based on their marked edges in various ways
*	Return: void, but generates a new set of triangles for the mesh
*	Called: called during the remesh() step
*/
void KnitGrapher::divide(std::unordered_set< glm::uvec2 > const &marked,       
	std::vector< glm::vec3 >& verts,										
	std::vector<std::vector<uint32_t>>& paths,							
	std::vector<glm::uvec3>& tris)										
{	
	//if no marked edges to divide by, return
	if (marked.empty()) {
		qDebug() << "Marked is empty!";
		return;
	}

	std::unordered_map< glm::uvec2, uint32_t > marked_verts;
	marked_verts.reserve(marked.size());

	std::vector< glm::ivec2 > edges(marked.begin(), marked.end());

	//lambda sort  
	std::sort(edges.begin(), edges.end(), [](glm::uvec2 const& a, glm::uvec2 const& b) {
		if (a.x != b.x) return a.x < b.x;
		else return a.y < b.y;
	});

	for (auto const& e : edges) {
		marked_verts.insert(std::make_pair(e, verts.size()));
		verts.emplace_back((verts[e.x] + verts[e.y]) / 2.0f);
	}

	for (auto& path : paths) {
		std::vector<uint32_t> new_path;
		new_path.emplace_back(path[0]);
		for (uint32_t i = 1; i < path.size(); ++i) {
			uint32_t v = lookup(path[i - 1], path[i], marked_verts);
			if (v != -1U) new_path.emplace_back(v);
			new_path.emplace_back(path[i]);
		}
		path = std::move(new_path);
	}

	std::vector<glm::uvec3> new_tris;

	for (auto const& tri : tris) {
		uint32_t a = tri.x;
		uint32_t b = tri.y;
		uint32_t c = tri.z;
		uint32_t ab = lookup(a, b, marked_verts);
		uint32_t bc = lookup(b, c, marked_verts);
		uint32_t ca = lookup(c, a, marked_verts);

		if (ab != -1U && bc != -1U && ca != -1U) {
			//1 -> 4 subdiv!
			new_tris.emplace_back(glm::uvec3(a, ab, ca));
			new_tris.emplace_back(glm::uvec3(b, bc, ab));
			new_tris.emplace_back(glm::uvec3(c, ca, bc));
			new_tris.emplace_back(glm::uvec3(ab, bc, ca));
		}
		else if (ab != -1U && bc != -1U && ca == -1U) {
			quad(new_tris, verts, a, ab, bc, c);
			new_tris.emplace_back(glm::uvec3(ab, b, bc));
		}
		else if (ab != -1U && bc == -1U && ca != -1U) {
			new_tris.emplace_back(glm::uvec3(a, ab, ca));
			quad(new_tris, verts, ab, b, c, ca);
		}
		else if (ab == -1U && bc != -1U && ca != -1U) {
			quad(new_tris, verts, a, b, bc, ca);
			new_tris.emplace_back(glm::uvec3(bc, c, ca));
		}
		else if (ab != -1U && bc == -1U && ca == -1U) {
			//1 -> 2 subdiv!
			new_tris.emplace_back(glm::uvec3(a, ab, c));
			new_tris.emplace_back(glm::uvec3(b, c, ab));
		}
		else if (ab == -1U && bc != -1U && ca == -1U) {
			new_tris.emplace_back(glm::uvec3(a, b, bc));
			new_tris.emplace_back(glm::uvec3(bc, c, a));
		}
		else if (ab == -1U && bc == -1U && ca != -1U) {
			new_tris.emplace_back(glm::uvec3(a, b, ca));
			new_tris.emplace_back(glm::uvec3(b, c, ca));
		}
		else {
			new_tris.emplace_back(glm::uvec3(a, b, c));
		}
	}
	//make the tris argument equal to the newly created triangle data
	tris = std::move(new_tris);
}


/*
*	Function: helper function that checks if the 3 points are in a counter clockwise orientation on a 2D plane
*	Return: bool, are the points coutner clockwise or not
*	Called: called during the unfold() function, as part of remesh()
*/
bool KnitGrapher::isCcw(glm::vec2 const& a, glm::vec2 const& b, glm::vec2 const& c) 
{
	return glm::dot(glm::vec2(-(b.y - a.y), (b.x - a.x)), c - a) > 0.0f;
};

/*
*	Function: get the distance of 2 vertices 'a' and 'b'
*	Return: float&, the reference to the distance, which will be infinity if the distance is not yet inside 'min_dis'
*	Called: called during the unfold() function, as part of remesh()
*/
float& KnitGrapher::get_dis(uint32_t a, uint32_t b, std::unordered_map< glm::uvec2, float >& min_dis) 
{
	if (a > b) std::swap(a, b);
	return min_dis.insert(std::make_pair(glm::uvec2(a, b), std::numeric_limits< float >::infinity())).first->second;
};


/*
*	Function: 'unfold' the 3D mesh onto a 2D plane, computing vertex distances and creating new triangles if the distance is too large
*	Return: void, but modifies the min_dis distance saving argument
*	Called: called during the remesh() function
*/
void KnitGrapher::unfold(uint32_t depth, uint32_t root, glm::vec2 const& flat_root,
	uint32_t ai, glm::vec2 const& flat_a,
	uint32_t bi, glm::vec2 const& flat_b,
	glm::vec2 const& limit_a, 
	glm::vec2 const& limit_b, std::unordered_map< glm::uvec2, uint32_t >  const& opposite,
	std::vector<glm::vec3> const& newVertices, std::unordered_map< glm::uvec2, float >& min_dis) 
{

	uint32_t ci;
	glm::vec2 flat_c;

	//if there is a triangle over the ai->bi edge, find other vertex and flatten it:
	{ 
		auto f = opposite.find(glm::uvec2(bi, ai));
		if (f == opposite.end()) return;
		ci = f->second;

		//figure out c's position along ab and distance from ab:
		glm::vec3 const& a = newVertices[ai];
		glm::vec3 const& b = newVertices[bi];
		glm::vec3 const& c = newVertices[ci];

		glm::vec3 ab = glm::normalize(b - a);
		float along = glm::dot(c - a, ab);
		float perp = -glm::length(c - a - ab * along);

		glm::vec2 flat_ab = glm::normalize(flat_b - flat_a);
		glm::vec2 flat_perp_ab = glm::vec2(-flat_ab.y, flat_ab.x);

		flat_c = flat_a + flat_ab * along + flat_perp_ab * perp;
	}

	bool ccw_rac = isCcw(flat_root, limit_a, flat_c) && isCcw(flat_root, flat_a, flat_c);
	bool ccw_rcb = isCcw(flat_root, flat_c, limit_b) && isCcw(flat_root, flat_c, flat_b);

	if (ccw_rac && ccw_rcb) {
		float& dis = get_dis(root, ci, min_dis);
		dis = std::min(dis, glm::length(flat_root - flat_c));

		if (depth > 1) {
			assert(isCcw(flat_root, flat_a, flat_c));
			unfold(depth - 1, root, flat_root, ai, flat_a, ci, flat_c, limit_a, flat_c, opposite, newVertices, min_dis);
			assert(isCcw(flat_root, flat_c, flat_b));
			unfold(depth - 1, root, flat_root, ci, flat_c, bi, flat_b, flat_c, limit_b, opposite, newVertices, min_dis);
		}
	}

	else if (ccw_rac && !ccw_rcb) {
		if (depth > 1) {
			unfold(depth - 1, root, flat_root, ai, flat_a, ci, flat_c, limit_a, limit_b, opposite, newVertices, min_dis);
		}
	}
	else if (!ccw_rac && ccw_rcb) {
		if (depth > 1) {
			//assert(is_ccw(flat_root, flat_c, flat_b));
			unfold(depth - 1, root, flat_root, ci, flat_c, bi, flat_b, limit_a, limit_b, opposite, newVertices, min_dis);
		}
	}
}

/*
*	Function: the first step of the 3D mesh -> KnitGraph algorithm, remeshes the original model into a higher detail one
*		done because the parameters may lead to edges of the model being too long, thus the triangles need to be split
*	Return: void, but modifies the newMesh class variable, along with some others
*	Called: called when the user selects the 'interpolate' option in the GUI
*/
void KnitGrapher::remesh()
{
	qDebug() << "Starting remesh()";

	//check if constraints are empty, no way to remesh without them
	if (constraints.size() == 0) {
		qDebug() << "No Constraints found! exiting remesh()...";
		return;
	}

	//extract edges from the model
	std::vector<std::vector<std::pair<uint32_t, float>>> adj(originalMesh.vertices.size());
	QSet <QPair<GLuint, GLuint>> edges;
	{ 
		std::set< std::pair< uint32_t, uint32_t > > edges;
		for (auto const& tri : originalMesh.triangles) {
			edges.insert(std::minmax(tri.x, tri.y));
			edges.insert(std::minmax(tri.y, tri.z));
			edges.insert(std::minmax(tri.z, tri.x));
		}
		for (auto const& e : edges) {
			float len = glm::length(originalMesh.vertices[e.second] - originalMesh.vertices[e.first]);
			adj[e.first].emplace_back(e.second, len);
			adj[e.second].emplace_back(e.first, len);
		}
	}


	//find chain paths on original model
	std::vector<std::vector<uint32_t>> paths;

	//altough done here, technically my constraints chain has nice 'single jump' vertex chains, so not really necessary
	//could be coopted in Visualizer for constraint by dijkstra adding...
	for (Constraint* constraint : constraints)
	{
		if (constraint->vertices.empty()) {
			qDebug() << "Constraint has no vertices, skipping...";
			continue;
		}
		std::vector<uint32_t> path;
		for (uint32_t goal : constraint->vertices) {
			if (path.empty()) {
				path.push_back(goal);
				continue;
			}
			std::vector< std::pair< float, uint32_t > > todo;
			std::vector< std::pair< float, uint32_t > > visited(originalMesh.vertices.size(), std::make_pair(std::numeric_limits< float >::infinity(), -1U));
			visit(todo, visited, goal, 0.0f, -1U);

			while (!todo.empty()) {
				std::pop_heap(todo.begin(), todo.end(), std::greater< std::pair< float, uint32_t > >());
				auto at = todo.back();
				todo.pop_back();
				if (at.first > visited[at.second].first) continue;
				if (at.second == path.back()) break;
				for (auto const &a : adj[at.second]) {
					visit(todo, visited, a.first, at.first + a.second, at.second);
				}
			}
			while (path.back() != goal) {
				if (visited[path.back()].second == -1) {
					qDebug() << "ERROR: constraint chain moves between connected components.";
					break;
				}
				path.emplace_back(visited[path.back()].second);
			}

		}
		paths.emplace_back(path);
	}


	//create a higher resolution mesh according to the parameters set by the user
	//in algorithm is part of the 'Parameters' class, i have those values inside this KnitGrapher class
	float maxEdgeLength = getMaxEdgeLength();
	float minEdgeRatio = 0.3f;											//allow user to set? probably not wise...
	float minEdgeRatioSquared = minEdgeRatio * minEdgeRatio;
	float maxEdgeLengthSquared = maxEdgeLength * maxEdgeLength;

	std::vector< glm::vec3 > verts = originalMesh.vertices;
	std::vector< glm::uvec3 > tris = originalMesh.triangles;

	//there is degenarate triangle checking in the original algorithm, will do the same here
	if (degenerateCheck(tris)) {
		qDebug() << "Degenerate triangle found, error...";
		return;
	}

	//lambda functions in original code were made into their own functions...
	while (true) {
		std::unordered_set< glm::uvec2 > marked;
		auto mark = [&marked](uint32_t a, uint32_t b) {
			if (b < a) std::swap(a, b);
			marked.insert(glm::uvec2(a, b));
		};
		auto is_marked = [&marked](uint32_t a, uint32_t b) {
			if (b < a) std::swap(a, b);
			return marked.find(glm::uvec2(a, b)) != marked.end();
		};
		(void)is_marked;
		(void)minEdgeRatioSquared;

		for (auto const& tri : tris) {
			float len_ab2 = glm::length2(verts[tri.y] - verts[tri.x]);
			float len_bc2 = glm::length2(verts[tri.z] - verts[tri.y]);
			float len_ca2 = glm::length2(verts[tri.x] - verts[tri.z]);
			if (len_ab2 > maxEdgeLengthSquared) mark(tri.x, tri.y);
			if (len_bc2 > maxEdgeLengthSquared) mark(tri.y, tri.z);
			if (len_ca2 > maxEdgeLengthSquared) mark(tri.z, tri.x);
		}
		if (marked.empty()) {
			break;
		}
		divide(marked, verts, paths, tris);
	}
	
	
	//once again the degenerate triangle check in original algortihm



	//extract edges from subdivided model:
	adj.assign(verts.size(), std::vector< std::pair< uint32_t, float > >());
	{ 
		std::set< std::pair< uint32_t, uint32_t > > edges;
		for (auto const& tri : tris) {
			edges.insert(std::minmax(tri.x, tri.y));
			edges.insert(std::minmax(tri.y, tri.z));
			edges.insert(std::minmax(tri.z, tri.x));
		}

		for (auto const& e : edges) {
			float len = glm::length(verts[e.second] - verts[e.first]);
			adj[e.first].emplace_back(e.second, len);
			adj[e.second].emplace_back(e.first, len);
		}
	}

	std::unordered_map< glm::uvec2, uint32_t > opposite; //vertex opposite each [oriented] triangle edge
	opposite.reserve(tris.size() * 3);
	for (auto const& tri : tris) {
		auto ret_xy = opposite.insert(std::make_pair(glm::uvec2(tri.x, tri.y), tri.z));
		assert(ret_xy.second);
		auto ret_yz = opposite.insert(std::make_pair(glm::uvec2(tri.y, tri.z), tri.x));
		assert(ret_yz.second);
		auto ret_zx = opposite.insert(std::make_pair(glm::uvec2(tri.z, tri.x), tri.y));
		assert(ret_zx.second);
	}

	
	{ //build (+ add to adj) extra "shortcut" edges by unwrapping triangle neighborhoods:
		std::unordered_map< glm::uvec2, float > min_dis;
		for (auto const& tri : tris) {
			glm::vec2 flat_x, flat_y, flat_z; //original verts
			glm::vec3 const& x = verts[tri.x];
			glm::vec3 const& y = verts[tri.y];
			glm::vec3 const& z = verts[tri.z];
			flat_x = glm::vec2(0.0f, 0.0f);
			flat_y = glm::vec2(glm::length(y - x), 0.0f);

			glm::vec3 xy = glm::normalize(y - x);
			glm::vec3 perp_xy = glm::normalize(glm::cross(glm::cross(y - x, z - x), y - x));
			float along = glm::dot(z - x, xy);
			float perp = glm::dot(z - x, perp_xy);

			flat_z = glm::vec2(along, perp);

			const constexpr GLuint D = 3; //depth to unfold triangles to for more adjacency information; makes slightly nicer geodesics at the expense of increased compute time.
			//THIS COULD BE A SETTABLE PARAMETER!
			//TODO: make this a parameter?

			if (D > 0) {
				unfold(D, tri.x, flat_x, tri.y, flat_y, tri.z, flat_z, flat_y, flat_z, opposite, verts, min_dis);
				unfold(D, tri.y, flat_y, tri.z, flat_z, tri.x, flat_x, flat_z, flat_x, opposite, verts, min_dis);
				unfold(D, tri.z, flat_z, tri.x, flat_x, tri.y, flat_y, flat_x, flat_y, opposite, verts, min_dis);
			}
		}
		for (uint32_t x = 0; x < verts.size(); ++x) {
			for (auto const& yd : adj[x]) {
				float& dis = get_dis(x, yd.first, min_dis);
				dis = std::min(dis, yd.second);
			}
		}

		//clear adj + re-create from min_dis:
		uint32_t old_adj = 0;
		for (auto const& a : adj) {
			old_adj += a.size();
		}

		adj.assign(verts.size(), std::vector< std::pair< uint32_t, float > >());

		
		for (auto const& xyd : min_dis) {
			assert(xyd.first.x != xyd.first.y);
			adj[xyd.first.x].emplace_back(xyd.first.y, xyd.second);
			adj[xyd.first.y].emplace_back(xyd.first.x, xyd.second);
		}

		uint32_t new_adj = 0;
		for (auto const& a : adj) {
			new_adj += a.size();
		}

		qDebug() << "Went from " << old_adj << " to " << new_adj << " by unfolding triangles.";

		//for consistency:
		for (auto& a : adj) {
			std::sort(a.begin(), a.end());
		}
	}
	

	//start creating embedded vertices
	std::vector< std::vector< EmbeddedVertex > > embedded_chains;
	
	//was a big for cycle, but that was because of radius checking, which is not part of my program
	for (auto const& cons : constraints) {
		if (cons->vertices.size() != 0) {   //checking for dummy constraints, dummy...
			embedded_chains.emplace_back();

			auto const& path = paths[&cons - &constraints[0]];
			//add directly to embedded constrained edges.
			for (auto v : path) {
				//assert(v < verts.size());
				embedded_chains.back().emplace_back(EmbeddedVertex::on_vertex(v));
			}
			continue;
		}
	}

	qDebug() << "embedded chains size" << embedded_chains.size() << " constraints size (with empty dummy):" << constraints.size();
	
	EmbeddedPlanarMap< float, SameValue< float >, ReplaceValue< float > > epm;
	
	uint32_t total_chain_edges = 0;
	for (uint32_t c = 0; c < embedded_chains.size(); ++c) {
		uint32_t first = 0;
		uint32_t last = 0;
		for (uint32_t i = 0; i + 1 < embedded_chains[c].size(); ++i) {
			uint32_t a = epm.add_vertex(embedded_chains[c][i]);
			uint32_t b = epm.add_vertex(embedded_chains[c][i+1]);
			epm.add_edge(a,b,constraints[c]->timeValue);
			++total_chain_edges;
			if (i == 0) first = a;
			if (i + 2 == embedded_chains[c].size()) last = b;
		}
		if (first != last) qDebug() << "NOTE: have open chain.";
	}
	uint32_t total_simplex_edges = 0;
	for (const auto &edges : epm.simplex_edges) {
		total_simplex_edges += edges.second.size();
	}

	qDebug() << "EPM has " << epm.vertices.size() << " vertices." ;
	qDebug() << "EPM has " << epm.simplex_vertices.size() << " simplices with vertices.";
	qDebug() << "EPM has " << epm.simplex_edges.size() << " simplices with edges (" << total_simplex_edges << " edges from " << total_chain_edges << " chain edges).";


	//Build a mesh that is split at the embedded edges:
	std::vector< EmbeddedVertex > split_evs;
	std::vector< glm::uvec3 > split_tris;
	std::vector< uint32_t > epm_to_split;

	epm.split_triangles(verts, tris, &split_evs, &split_tris, &epm_to_split);

	std::vector< glm::vec3 > split_verts;

	split_verts.reserve(split_evs.size());
	for (auto const& ev : split_evs) {
		split_verts.emplace_back(ev.interpolate(verts));
	}

	//record constrained edges in terms of split_verts:
	std::unordered_map< glm::uvec2, float > constrained_edges;
	std::vector< float > split_values(split_verts.size(), std::numeric_limits< float >::quiet_NaN());
	for (const auto& se : epm.simplex_edges) {
		for (auto const& e : se.second) {
			glm::uvec2 ab = glm::uvec2(epm_to_split[e.first], epm_to_split[e.second]);
			if (ab.x > ab.y) std::swap(ab.x, ab.y);
			constrained_edges.insert(std::make_pair(ab, e.value));
			//also grab vertex values:
			split_values[epm_to_split[e.first]] = e.value;
			split_values[epm_to_split[e.second]] = e.value;
		}
	}
	qDebug() << constrained_edges.size() << " constrained edges.";


	std::vector< uint32_t > tri_component(split_tris.size(), -1U);
	std::vector< bool > component_keep;
	{ //mark connected components + delete the "wrong" ones
		std::unordered_map< glm::uvec2, uint32_t > over;
		for (const auto& tri : split_tris) {
			uint32_t ti = &tri - &split_tris[0];
			auto res = over.insert(std::make_pair(glm::uvec2(tri.x, tri.y), ti));
			//assert(res.second);
			res = over.insert(std::make_pair(glm::uvec2(tri.y, tri.z), ti));
			//assert(res.second);
			res = over.insert(std::make_pair(glm::uvec2(tri.z, tri.x), ti));
			//assert(res.second);
		}
		for (uint32_t seed = 0; seed < split_tris.size(); ++seed) {
			if (tri_component[seed] != -1U) continue;
			//std::cout << "Doing CC with seed " << seed << std::endl; //DEBUG
			uint32_t component = component_keep.size();
			tri_component[seed] = component;
			std::set< float > values;
			std::vector< uint32_t > todo;
			todo.emplace_back(seed);
			auto do_edge = [&](uint32_t a, uint32_t b) {
				{ //if edge is constrained, don't traverse over:
					glm::uvec2 e(a, b);
					if (e.x > e.y) std::swap(e.x, e.y);
					auto v = constrained_edges.find(e);
					if (v != constrained_edges.end()) {
						values.insert(v->second);
						return;
					}
				}
				//otherwise, traverse over:
				auto f = over.find(glm::uvec2(b, a));
				if (f != over.end()) {
					if (tri_component[f->second] != component) {
						assert(tri_component[f->second] == -1U);
						tri_component[f->second] = component;
						todo.emplace_back(f->second);
					}
				}
			};
			while (!todo.empty()) {
				uint32_t at = todo.back();
				todo.pop_back();
				assert(tri_component[at] == component);
				do_edge(split_tris[at].x, split_tris[at].y);
				do_edge(split_tris[at].y, split_tris[at].z);
				do_edge(split_tris[at].z, split_tris[at].x);
			}
			component_keep.emplace_back(values.size() > 1);
		}
		qDebug() << "Have " << component_keep.size() << " connected components.";
	}

	//remove any split_verts that aren't used:
	std::vector< glm::vec3 > compressed_verts;
	std::vector< float > compressed_values;
	std::vector< glm::uvec3 > compressed_tris;
	for (uint32_t ti = 0; ti < split_tris.size(); ++ti) {
		if (component_keep[tri_component[ti]]) {
			compressed_tris.emplace_back(split_tris[ti]);
		}
	}
	compressed_verts.reserve(split_verts.size());
	std::vector< uint32_t > to_compressed(split_verts.size(), -1U);
	auto add_vert = [&](uint32_t vi) {
		if (to_compressed[vi] == -1U) {
			to_compressed[vi] = compressed_verts.size();
			compressed_verts.emplace_back(split_verts[vi]);
			compressed_values.emplace_back(split_values[vi]);
		}
		return to_compressed[vi];
	};
	for (auto& tri : compressed_tris) {
		tri.x = add_vert(tri.x);
		tri.y = add_vert(tri.y);
		tri.z = add_vert(tri.z);
	}
	qDebug() << "Went from " << tris.size() << " to (via split) " << split_tris.size() << " to (via discard) " << compressed_tris.size() << " triangles."; //DEBUG

	newMesh.vertices = compressed_verts;
	newMesh.triangles = compressed_tris;
	
	qDebug() << compressed_values.size() << " constrained values and " << compressed_verts.size() << " vertices";

	constrained_values = compressed_values;
}



/*
*	Function: find the first active chains that will start the knitgraph construction on the mesh
*		basically find the 'origin point' from which the knitting should start on the mesh
*	Return: void, but modifies many variables defined in the class, mainly creates a RowColGraph
*	Called: called when the user selects the 'step' button in the GUI, but only during the first click
*/
void KnitGrapher::findFirstActiveChains(std::vector< std::vector< EmbeddedVertex > >* active_chains_,
	std::vector< std::vector< Stitch > >* active_stitches_,
	RowColGraph* graph_) 
{
	
	assert(active_chains_);
	auto& active_chains = *active_chains_;
	active_chains.clear();

	assert(active_stitches_);
	auto& active_stitches = *active_stitches_;
	active_stitches.clear();

	//triangles must reference valid time values:
	for (glm::uvec3 const& tri : newMesh.triangles) {
		assert(tri.x < constrained_values.size());
		assert(tri.y < constrained_values.size());
		assert(tri.z < constrained_values.size());
	}


	//find boundary loops:
	//build half-edge -> other vertex map:
	std::unordered_map< glm::uvec2, uint32_t > next_vertex;
	{
		auto do_edge = [&next_vertex](uint32_t a, uint32_t b, uint32_t c) {
			assert(a != b && a != c && b != c);
			auto ret = next_vertex.insert(std::make_pair(glm::uvec2(a, b), c));
			if (!ret.second) {
				qDebug() << "ERROR: non-manifold mesh [or inconsistent orientation] -- directed edge appears twice.";
			}
		};
		for (glm::uvec3 const& tri : newMesh.triangles) {
			do_edge(tri.x, tri.y, tri.z);
			do_edge(tri.y, tri.z, tri.x);
			do_edge(tri.z, tri.x, tri.y);
		}
	}

	//extract boundary (== unpaired half-edges):
	std::unordered_map< uint32_t, uint32_t > boundary;
	{
		for (const auto& ev : next_vertex) {
			if (next_vertex.count(glm::uvec2(ev.first.y, ev.first.x))) continue; //skip paired edges
			auto ret = boundary.insert(std::make_pair(ev.first.x, ev.first.y)); //half-edge is on the boundary
			if (!ret.second) {
				qDebug() << "ERROR: non-manifold mesh [or inconsistent orientation] -- vertex is source of multiple boundary edges.";
			}
		}
	}

	while (!boundary.empty()) {
		std::vector< uint32_t > chain;
		chain.emplace_back(boundary.begin()->first);
		chain.emplace_back(boundary.begin()->second);
		boundary.erase(boundary.begin());
		do {
			auto f = boundary.find(chain.back());
			assert(f != boundary.end());
			chain.emplace_back(f->second);
			boundary.erase(f);
		} while (chain.back() != chain[0]);

		//to be marked active, chain must be constant value and < adjacent time values.
		float chain_min = std::numeric_limits< float >::infinity();
		float chain_max = -std::numeric_limits< float >::infinity();
		float adj_min = std::numeric_limits< float >::infinity();
		float adj_max = -std::numeric_limits< float >::infinity();
		for (uint32_t i = 0; i + 1 < chain.size(); ++i) {
			chain_min = std::min(chain_min, constrained_values[chain[i]]);
			chain_max = std::max(chain_max, constrained_values[chain[i]]);
			auto f = next_vertex.find(glm::uvec2(chain[i], chain[i + 1]));
			assert(f != next_vertex.end());
			adj_min = std::min(adj_min, constrained_values[f->second]);
			adj_max = std::max(adj_max, constrained_values[f->second]);
		}

		if (chain_min != chain_max) {
			qDebug() << "WARNING: discarding chain with non-constant value range [" << chain_min << ", " << chain_max << "].";
			continue;
		}
		//the 1e-3 is to add some slop in case of somewhat noisy interpolation
		if (!(adj_min > chain_max - 1e-3)) {
			if (adj_max < chain_min) {
				//this is a maximum chain
			}
			else { //this is a mixed min/max chain; weird
				qDebug() << "WARNING: discarding chain with value range [" << chain_min << ", " << chain_max << "] because neighbors have value range [" << adj_min << ", " << adj_max << "].";
			}
			continue;
		}

		//sample chain to get output active chain:
		std::vector< EmbeddedVertex > embedded_chain;
		embedded_chain.reserve(chain.size());

		for (auto c : chain) {
			embedded_chain.emplace_back(EmbeddedVertex::on_vertex(c));
		}

		//subdivide the chain to make sure there are enough samples for later per-sample computations:
		std::vector< EmbeddedVertex > divided_chain;
		sampleChain(getChainSampleSpacing(), embedded_chain, &divided_chain);

		//further subdivide and place stitches:
		float total_length = 0.0f;
		for (uint32_t ci = 1; ci < divided_chain.size(); ++ci) {
			glm::vec3 a = divided_chain[ci - 1].interpolate(newMesh.vertices);
			glm::vec3 b = divided_chain[ci].interpolate(newMesh.vertices);
			total_length += glm::length(b - a);
		}

		float stitch_width = stitchWidth / modelUnitLength;
		uint32_t stitches = std::max(3, int32_t(std::round(total_length / stitch_width)));

		active_chains.emplace_back(divided_chain);
		active_stitches.emplace_back();
		active_stitches.back().reserve(stitches);
		for (uint32_t s = 0; s < stitches; ++s) {
			active_stitches.back().emplace_back((s + 0.5f) / float(stitches), Stitch::FlagLinkAny);
		}
	}

	assert(active_chains.size() == active_stitches.size());
	qDebug() << "Found " << active_chains.size() << " first active chains.";

	if (graph_) {
		for (uint32_t ci = 0; ci < active_chains.size(); ++ci) {
			auto const& chain = active_chains[ci];

			std::vector< float > lengths;
			lengths.reserve(chain.size());
			lengths.emplace_back(0.0f);
			for (uint32_t i = 1; i < chain.size(); ++i) {
				glm::vec3 a = chain[i - 1].interpolate(newMesh.vertices);
				glm::vec3 b = chain[i].interpolate(newMesh.vertices);
				lengths.emplace_back(lengths.back() + glm::length(b - a));
			}
			assert(lengths.size() == chain.size());

			auto li = lengths.begin();
			for (auto& s : active_stitches[ci]) {
				float l = lengths.back() * s.t;

				while (li != lengths.end() && *li <= l) ++li;
				assert(li != lengths.begin());
				assert(li != lengths.end());

				float m = (l - *(li - 1)) / (*li - *(li - 1));
				uint32_t i = li - lengths.begin();

				assert(s.vertex == -1U);
				s.vertex = graph_->vertices.size();
				graph_->vertices.emplace_back();
				graph_->vertices.back().at = EmbeddedVertex::mix(
					chain[i - 1], chain[i], m
				);
			}

			uint32_t prev = (chain[0] == chain.back() ? active_stitches[ci].back().vertex : -1U);
			for (auto& s : active_stitches[ci]) {
				if (prev != -1U) {
					assert(prev < graph_->vertices.size());
					assert(s.vertex < graph_->vertices.size());
					assert(graph_->vertices[prev].row_out == -1U);
					graph_->vertices[prev].row_out = s.vertex;
					assert(graph_->vertices[s.vertex].row_in == -1U);
					graph_->vertices[s.vertex].row_in = prev;
				}
				prev = s.vertex;
			}
		}
	}

}

/*
*	Function: sample points along the edges of a chain of embedded vertices at regular intervals based on 'spacing'
*	Return: void, but also modifies the sampled_chain_ variable
*	Called: called during peelSlice() as well as findFirstActiveChains()
*/
void KnitGrapher::sampleChain(float spacing,
	std::vector< EmbeddedVertex > const& chain, //in: chain to be sampled
	std::vector< EmbeddedVertex >* sampled_chain_)
{
	
	auto& sampled_chain = *sampled_chain_;
	sampled_chain.clear();

	for (uint32_t ci = 0; ci + 1 < chain.size(); ++ci) {
		sampled_chain.emplace_back(chain[ci]);

		glm::vec3 a = chain[ci].interpolate(newMesh.vertices);
		glm::vec3 b = chain[ci + 1].interpolate(newMesh.vertices);
		glm::uvec3 common = EmbeddedVertex::common_simplex(chain[ci].simplex, chain[ci + 1].simplex);
		glm::vec3 wa = chain[ci].weights_on(common);
		glm::vec3 wb = chain[ci + 1].weights_on(common);

		float length = glm::length(b - a);
		int32_t insert = std::floor(length / spacing);
		for (int32_t i = 0; i < insert; ++i) {
			float m = float(i + 1) / float(insert + 1);
			sampled_chain.emplace_back(common, glm::mix(wa, wb, m));
		}
	}
	sampled_chain.emplace_back(chain.back());
}


/*
*	Function: linearly interpolate time values on the mesh as given by constrained_values
*	Return: void, but modifies constrained_values time values
*	Called: called after the remesh() function remeshes the 3D mesh and generates the time values
*/
void KnitGrapher::interpolateValues() 
{
	if (constrained_values.size() != newMesh.vertices.size()) {
		qDebug() << "Constraints size does not match model vertices size" << constrained_values.size() << newMesh.vertices.size();
		return;
	}

	qDebug() << newMesh.vertices.size() << " vertices and " << newMesh.triangles.size() << " triangles.";
	auto& values = constrained_values;

	std::vector< uint32_t > dofs;
	dofs.reserve(constrained_values.size());
	uint32_t total_dofs = 0;
	for (auto c : constrained_values) {
		if (c == c) dofs.emplace_back(-1U);
		else dofs.emplace_back(total_dofs++);
	}

	qDebug() << "Have " << total_dofs << " degrees of freedom and " << (constrained_values.size() - total_dofs) << " constraints.";

	if (total_dofs == constrained_values.size()) {
		qDebug() << "Cannot interpolate from no constraints.";
	}

	std::map< std::pair< uint32_t, uint32_t >, float > edge_weights;

	for (const auto& tri : newMesh.triangles) {
		const glm::vec3& a = newMesh.vertices[tri.x];
		const glm::vec3& b = newMesh.vertices[tri.y];
		const glm::vec3& c = newMesh.vertices[tri.z];

		float weight_ab = glm::dot(a - c, b - c) / glm::length(glm::cross(a - c, b - c));
		float weight_bc = glm::dot(b - a, c - a) / glm::length(glm::cross(b - a, c - a));
		float weight_ca = glm::dot(c - b, a - b) / glm::length(glm::cross(c - b, a - b));

		edge_weights.insert(std::make_pair(std::minmax(tri.x, tri.y), 0.0f)).first->second += weight_ab;
		edge_weights.insert(std::make_pair(std::minmax(tri.y, tri.z), 0.0f)).first->second += weight_bc;
		edge_weights.insert(std::make_pair(std::minmax(tri.z, tri.x), 0.0f)).first->second += weight_ca;
	}

	//turn edge weights vector into adjacency lists:
	std::vector< std::vector< std::pair< uint32_t, float > > > adj(newMesh.vertices.size());
	for (const auto& ew : edge_weights) {
		adj[ew.first.first].emplace_back(ew.first.second, ew.second);
		adj[ew.first.second].emplace_back(ew.first.first, ew.second);
	}


	std::vector< Eigen::Triplet< double > > coefficients;
	coefficients.reserve(newMesh.triangles.size() * 3);
	Eigen::VectorXd rhs(total_dofs);

	for (uint32_t i = 0; i < dofs.size(); ++i) {
		if (dofs[i] == -1U) continue;
		float sum = 0.0f;
		float one = 0.0f;
		for (auto a : adj[i]) {
			if (dofs[a.first] == -1U) {
				one += a.second * constrained_values[a.first];
			}
			else {
				coefficients.emplace_back(dofs[i], dofs[a.first], a.second);
			}
			sum += a.second;
		}
		coefficients.emplace_back(dofs[i], dofs[i], -sum);
		rhs[dofs[i]] = -one;
	}

	Eigen::SparseMatrix< double > A(total_dofs, total_dofs);
	A.setFromTriplets(coefficients.begin(), coefficients.end());
	A.makeCompressed(); 

	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > solver;	
	solver.compute(A);
	if (solver.info() != Eigen::Success) {
		qDebug() << "Decomposition failed.";
		return;
	}
	Eigen::VectorXd x = solver.solve(rhs);
	if (solver.info() != Eigen::Success) {
		qDebug() << "Solving failed.";
		return;
	}
	
	values = constrained_values;
	for (uint32_t i = 0; i < dofs.size(); ++i) {
		if (dofs[i] != -1U) values[i] = x[dofs[i]];
	}

	qDebug() << "Interpolation complete!";
	for (float val : values) {
		qDebug() << val;
	}
}

/*
*	Function: helper function that converts an AutoKnitMesh's 'triangles' class variable into ObjectMesh's 'indices' class variable
*	Return: a vector array of int values, representing different vertices on the triangle
*	Called: called whenever an ObjectMesh should be created and passed to the Visualizer class, such as during peelSlice()
*/
std::vector<GLuint> KnitGrapher::toIntArray(std::vector<glm::uvec3> triangles) 
{
	std::vector<GLuint> result;

	for (auto t : triangles) {
		result.push_back(t.x);
		result.push_back(t.y);
		result.push_back(t.z);
	}
	return result;
}

/*
*	Function: a checking function that finds if a there is a degenerate triangle in a set
*	Return: bool, is there such a triangle or not
*	Called: called during remesh, implies a fault of the loaded .obj file...
*/
bool KnitGrapher::degenerateCheck(std::vector<glm::uvec3> tris) 
{
	for (auto const& tri : tris) {
		glm::vec3 const& x = originalMesh.vertices[tri.x];
		glm::vec3 const& y = originalMesh.vertices[tri.y];
		glm::vec3 const& z = originalMesh.vertices[tri.z];

		if (tri.x == tri.y || tri.y == tri.z || tri.z == tri.x)
		{
			return true;
		}
		// check if indices are degenerate
		if (x == y || x == z || y == z)
		{
			return true;
		}

	}
	return false;
}

/*
*	Function: extract chains of edges at a given level from the mesh based on vertex values. 
*	Return: void, but generates a vector of EmbeddedVertex objects
*	Called: called duing peelSlice()
*/
void KnitGrapher::extractLevelChains(
	AutoKnitMesh const& model, //in: model on which to embed vertices
	std::vector< float > const& values, //in: values at vertices
	float const level, //in: level at which to extract chains
	std::vector< std::vector< EmbeddedVertex > >* chains_ ) //chains of edges at given level
{
	assert(chains_);
	auto& chains = *chains_;
	chains.clear();

	//embed points along all edges that start below level and end at or above it:
	std::vector< glm::vec3 > const& verts = model.vertices;
	std::vector< glm::uvec3 > const& tris = model.triangles;

	std::unordered_map< glm::uvec2, EmbeddedVertex > embedded_pts;
	std::unordered_map< glm::uvec2, glm::vec3 > pts;
	auto add = [&](uint32_t a, uint32_t b) {
		assert(values[a] < level && values[b] >= level);
		float mix = (level - values[a]) / (values[b] - values[a]);
		pts[glm::uvec2(a, b)] = glm::mix(verts[a], verts[b], mix);
		embedded_pts[glm::uvec2(a, b)] = EmbeddedVertex::on_edge(a, b, mix);
		return glm::uvec2(a, b);
	};
	std::unordered_map< glm::uvec2, glm::uvec2 > links;
	std::unordered_map< glm::uvec2, glm::uvec2 > back_links;
	auto link = [&links, &back_links](glm::uvec2 f, glm::uvec2 t) {
		auto res = links.insert(std::make_pair(f, t));
		assert(res.second);
		auto res2 = back_links.insert(std::make_pair(t, f));
		assert(res2.second);
	};
	for (auto const& tri : tris) {
		uint32_t a = tri.x;
		uint32_t b = tri.y;
		uint32_t c = tri.z;
		//spin triangle until 'a' is the minimum distance value:
		for (uint32_t i = 0; i < 3; ++i) {
			if (values[a] <= values[b] && values[a] <= values[c]) break;
			uint32_t t = a; a = b; b = c; c = t;
		}

		if (values[a] >= level) continue; //all above border
		assert(values[a] < level);

		if (values[b] >= level && values[c] >= level) {
			//edge is from ca to ab
			link(add(a, c), add(a, b));
		}
		else if (values[b] >= level && values[c] < level) {
			//edge is from bc to ab
			link(add(c, b), add(a, b));
		}
		else if (values[b] < level && values[c] >= level) {
			//edge is from ca to bc
			link(add(a, c), add(b, c));
		}
		else {
			assert(values[b] < level && values[c] < level);
			//all below border, nothing to do.
		}
	}

	uint32_t found_chains = 0;
	uint32_t found_loops = 0;

	//read back path from links:
	while (!links.empty()) {
		std::deque< glm::uvec2 > loop;
		loop.emplace_back(links.begin()->first);
		loop.emplace_back(links.begin()->second);

		//remove seed link:
		links.erase(links.begin());
		{
			auto b = back_links.find(loop.back());
			assert(b != back_links.end());
			assert(b->second == loop[0]);
			back_links.erase(b);
		}

		//extend forward:
		while (true) {
			auto f = links.find(loop.back());
			if (f == links.end()) break;
			loop.emplace_back(f->second);
			//remove link:
			auto b = back_links.find(loop.back());
			assert(b != back_links.end());
			assert(b->second == loop[loop.size() - 2]);
			links.erase(f);
			back_links.erase(b);
		}
		//extend backward:
		while (true) {
			auto b = back_links.find(loop[0]);
			if (b == back_links.end()) break;
			loop.emplace_front(b->second);
			//remove link:
			auto f = links.find(loop.front());
			assert(f != links.end());
			assert(f->second == loop[1]);
			back_links.erase(b);
			links.erase(f);
		}

		if (loop.front() == loop.back()) ++found_loops;
		else ++found_chains;

		chains.emplace_back();
		chains.back().reserve(loop.size());
		for (glm::uvec2 e : loop) {
			auto f = embedded_pts.find(e);
			assert(f != embedded_pts.end());
			chains.back().emplace_back(f->second);
		}
	}

	qDebug() << "extract_level_chains found " << found_loops << " loops and " << found_chains << " chains." ;
}

/*
*	Function: trims the mesh so that only relevant portion of it is left
*	Return: void, but generates the trimmed mesh and passes it through the 'clipped' argument, also generates the set of Embedded vertices
*		optionally can also return the vertices left of, and right of the clipped mesh
*	Called: called during the peelSlice() function
*/
void KnitGrapher::trimModel(std::vector< std::vector< EmbeddedVertex > >& left_of,
	std::vector< std::vector< EmbeddedVertex > >& right_of,
	AutoKnitMesh* clipped_,
	std::vector< EmbeddedVertex >* clipped_vertices_,
	std::vector< std::vector< uint32_t > >* left_of_vertices_, 
	std::vector< std::vector< uint32_t > >* right_of_vertices_) 
{
	assert(clipped_);
	auto& clipped = *clipped_;
	clipped.clear();

	assert(clipped_vertices_);
	auto& clipped_vertices = *clipped_vertices_;
	clipped_vertices.clear();

	{ //PARANOIA: make sure all chains are loops or edge-to-edge:
		std::unordered_set< glm::uvec2 > edges;
		for (auto const& tri : newMesh.triangles) {
			auto do_edge = [&edges](uint32_t a, uint32_t b) {
				if (a > b) std::swap(a, b);
				auto ret = edges.insert(glm::uvec2(a, b));
				if (!ret.second) edges.erase(ret.first);
			};
			do_edge(tri.x, tri.y);
			do_edge(tri.y, tri.z);
			do_edge(tri.z, tri.x);
		}
		std::unordered_set< uint32_t > edge_verts;
		for (auto const& e : edges) {
			edge_verts.insert(e.x);
			edge_verts.insert(e.y);
		}

		auto on_edge = [&edges, &edge_verts](EmbeddedVertex const& ev) -> bool {
			if (ev.simplex.z != -1U) {
				return false;
			}
			else if (ev.simplex.y != -1U) {
				return edges.count(glm::uvec2(ev.simplex.x, ev.simplex.y)) != 0;
			}
			else {
				return edge_verts.count(ev.simplex.x) != 0;
			}
		};

		for (auto const& chain : left_of) {
			assert(chain.size() >= 2);
			if (chain[0] == chain.back()) continue;
			assert(on_edge(chain[0]));
			assert(on_edge(chain.back()));
		}

		for (auto const& chain : right_of) {
			assert(chain.size() >= 2);
			if (chain[0] == chain.back()) continue;
			assert(on_edge(chain[0]));
			assert(on_edge(chain.back()));
		}
	}

	globedges.clear(); //global list used for edge tracking in epm; awkward but should work.
	EmbeddedPlanarMap< Value, Value::Reverse, Value::Combine, Value::Split > epm;
	std::unordered_set< glm::uvec2 > empty_edges; //when it's the same vertex after rounding
	uint32_t total_chain_edges = 0;
	uint32_t fresh_id = 0;
	for (auto & chain : left_of) {
		uint32_t prev = epm.add_vertex(chain[0]);
		uint32_t prev_id = fresh_id++;
		for (uint32_t i = 1; i < chain.size(); ++i) {
			uint32_t cur = epm.add_vertex(chain[i]);
			uint32_t cur_id = fresh_id++;
			if (prev == cur) {
				//std::cout << "NOTE: vertex " << chain[i - 1] << " and " << chain[i] << " (in a left_of chain) round to the same value." << std::endl; //DEBUG
				empty_edges.insert(glm::uvec2(prev_id, cur_id));
			}
			Value value;
			value.sum = 1;
			globedges.emplace_back(Edge::Initial, prev_id, cur_id);
			value.edge = globedges.size() - 1;
			epm.add_edge(prev, cur, value);
			prev = cur;
			prev_id = cur_id;
			++total_chain_edges;
		}
	}



	for (auto & chain : right_of) {
		uint32_t prev = epm.add_vertex(chain[0]);
		uint32_t prev_id = fresh_id++;
		for (uint32_t i = 1; i < chain.size(); ++i) {
			uint32_t cur = epm.add_vertex(chain[i]);
			uint32_t cur_id = fresh_id++;
			if (prev == cur) {
				//std::cout << "NOTE: vertex " << chain[i - 1] << " and " << chain[i] << " (in a right_of chain) round to the same value." << std::endl; //DEBUG
				empty_edges.insert(glm::uvec2(prev_id, cur_id));
			}
			Value value;
			value.sum = (1 << 8);
			globedges.emplace_back(Edge::Initial, cur_id, prev_id);
			value.edge = globedges.size() - 1;
			epm.add_edge(cur, prev, value);
			prev = cur;
			prev_id = cur_id;
			++total_chain_edges;
		}
	}

	uint32_t total_simplex_edges = 0;
	for (const auto& edges : epm.simplex_edges) {
		total_simplex_edges += edges.second.size();
	}
	qDebug() << "EPM has " << epm.vertices.size() << " vertices.";
	qDebug() << "EPM has " << epm.simplex_vertices.size() << " simplices with vertices.";
	qDebug() << "EPM has " << epm.simplex_edges.size() << " simplices with edges (" << total_simplex_edges << " edges from " << total_chain_edges << " chain edges).";

	std::vector< std::vector< uint32_t > > left_of_epm, right_of_epm;
	{
		//for each original id->id edge, extract the chain of vertices that it expands to:
		struct Source {
			Source(uint32_t a_, uint32_t b_, uint32_t num_, uint32_t den_) : a(a_), b(b_), num(num_), den(den_) { }
			uint32_t a, b; //a / b vertices
			uint32_t num, den; //position along (subdivided?) edge
		};
		std::vector< std::vector< Source > > sources;
		sources.reserve(globedges.size());
		for (uint32_t e = 0; e < globedges.size(); ++e) {
			assert(e == sources.size());
			sources.emplace_back();
			if (globedges[e].type == Edge::Initial) {
				assert(globedges[e].a != globedges[e].b);
				sources.back().emplace_back(globedges[e].a, globedges[e].b, 1, 2);
			}
			else if (globedges[e].type == Edge::Reverse) {
				assert(globedges[e].a == globedges[e].b);
				assert(globedges[e].a < e);
				for (auto const& s : sources[globedges[e].a]) {
					sources.back().emplace_back(s.b, s.a, s.den - s.num, s.den);
				}
			}
			else if (globedges[e].type == Edge::Combine) {
				assert(globedges[e].a < e);
				assert(globedges[e].b < e);
				for (auto const& s : sources[globedges[e].a]) {
					sources.back().emplace_back(s.a, s.b, s.num, s.den);
				}
				for (auto const& s : sources[globedges[e].b]) {
					sources.back().emplace_back(s.a, s.b, s.num, s.den);
				}
			}
			else if (globedges[e].type == Edge::SplitFirst) {
				assert(globedges[e].a == globedges[e].b);
				assert(globedges[e].a < e);
				for (auto const& s : sources[globedges[e].a]) {
					sources.back().emplace_back(s.a, s.b, s.num * 2 - 1, s.den * 2);
				}
			}
			else if (globedges[e].type == Edge::SplitSecond) {
				assert(globedges[e].a == globedges[e].b);
				assert(globedges[e].a < e);
				for (auto const& s : sources[globedges[e].a]) {
					assert(uint32_t(s.den * 2) > s.den); //make sure we're not overflowing
					sources.back().emplace_back(s.a, s.b, s.num * 2 + 1, s.den * 2);
				}
			}
			else {
				assert(false);
			}
		}

		struct SubEdge {
			SubEdge(uint32_t a_, uint32_t b_, uint32_t num_, uint32_t den_) : a(a_), b(b_), num(num_), den(den_) { }
			uint32_t a, b; //a / b vertices (epm indices)
			uint32_t num, den; //position of center
		};

		//now iterate actual epm edges and check where they end up mapping.
		std::unordered_map< glm::uvec2, std::vector< SubEdge > > edge_subedges; //<-- indexed by 'vertex_id' values, contains epm.vertices indices
		for (auto const& se : epm.simplex_edges) {
			for (auto const& ee : se.second) {
				assert(ee.value.edge < sources.size());
				assert(!sources[ee.value.edge].empty());
				//copy subedge(s) corresponding to this edge to their source edges:
				for (auto const& s : sources[ee.value.edge]) {
					assert(s.a != s.b);
					if (s.a < s.b) {
						edge_subedges[glm::uvec2(s.a, s.b)].emplace_back(ee.first, ee.second, s.num, s.den);
					}
					else {
						edge_subedges[glm::uvec2(s.b, s.a)].emplace_back(ee.second, ee.first, s.den - s.num, s.den);
					}
				}
			}
		}
		//make sure subedge chains are logically sorted:
		for (auto& ese : edge_subedges) {
			std::vector< SubEdge >& subedges = ese.second;
			assert(!subedges.empty());
			std::stable_sort(subedges.begin(), subedges.end(), [](SubEdge const& a, SubEdge const& b) {
				return uint64_t(a.num) * uint64_t(b.den) < uint64_t(b.num) * uint64_t(a.den);
				});
			for (uint32_t i = 0; i + 1 < subedges.size(); ++i) {
				assert(subedges[i].b == subedges[i + 1].a);
			}
		}

		//okay, now read back vertex chains by querying by ID values:

		uint32_t fresh_id = 0;
		left_of_epm.reserve(left_of.size());
		for (auto const& chain : left_of) {
			left_of_epm.emplace_back();
			std::vector< uint32_t >& epm_chain = left_of_epm.back();
			uint32_t prev_id = fresh_id++;
			for (uint32_t i = 1; i < chain.size(); ++i) {
				uint32_t cur_id = fresh_id++;
				auto f = edge_subedges.find(glm::uvec2(prev_id, cur_id));
				if (empty_edges.count(glm::uvec2(prev_id, cur_id))) {
					assert(f == edge_subedges.end());
				}
				else {
					assert(f != edge_subedges.end());
					for (auto const& se : f->second) {
						if (epm_chain.empty()) epm_chain.emplace_back(se.a);
						else assert(epm_chain.back() == se.a);
						epm_chain.emplace_back(se.b);
					}
				}
				prev_id = cur_id;
			}
			assert((chain[0] == chain.back()) == (epm_chain[0] == epm_chain.back()));
		}
		right_of_epm.reserve(right_of.size());
		for (auto const& chain : right_of) {
			right_of_epm.emplace_back();
			std::vector< uint32_t >& epm_chain = right_of_epm.back();
			uint32_t prev_id = fresh_id++;
			for (uint32_t i = 1; i < chain.size(); ++i) {
				uint32_t cur_id = fresh_id++;
				auto f = edge_subedges.find(glm::uvec2(prev_id, cur_id));
				if (empty_edges.count(glm::uvec2(prev_id, cur_id))) {
					assert(f == edge_subedges.end());
				}
				else {
					assert(f != edge_subedges.end());
					for (auto const& se : f->second) {
						if (epm_chain.empty()) epm_chain.emplace_back(se.a);
						else assert(epm_chain.back() == se.a);
						epm_chain.emplace_back(se.b);
					}
				}
				prev_id = cur_id;
			}
			assert((chain[0] == chain.back()) == (epm_chain[0] == epm_chain.back()));
		}
	}
	auto cleanup_chain = [&](std::vector< uint32_t >& epm_chain, int32_t value) {
		//want chain to be loop-free --> look for loops!
		float length_removed = 0.0f;
		uint32_t verts_removed = 0;

		//useful:
		auto remove_from_epm = [&epm_chain, &epm, value](uint32_t first, uint32_t last) {
			assert(first <= last);
			assert(last < epm_chain.size());
			//remove value of all segments [first,last] from planar map:
			for (uint32_t i = first; i + 1 <= last; ++i) {
				auto& a = epm.vertices[epm_chain[i]];
				auto& b = epm.vertices[epm_chain[i + 1]];
				glm::uvec3 common = IntegerEmbeddedVertex::common_simplex(a.simplex, b.simplex);
				auto f = epm.simplex_edges.find(common);
				assert(f != epm.simplex_edges.end());
				bool found = false;
				for (auto& e : f->second) {
					if (e.first == epm_chain[i] && e.second == epm_chain[i + 1]) {
						e.value.sum -= value;
						found = true;
						break;
					}
					else if (e.second == epm_chain[i] && e.first == epm_chain[i + 1]) {
						e.value.sum += value;
						found = true;
						break;
					}
				}
				assert(found);
			}
		};

		float initial_length = std::numeric_limits< float >::quiet_NaN();

		bool again = true;
		while (again) {
			again = false;

			std::vector< float > lengths;
			lengths.reserve(epm_chain.size());
			lengths.emplace_back(0.0f);
			for (uint32_t i = 1; i < epm_chain.size(); ++i) {
				glm::vec3 a = epm.vertices[epm_chain[i - 1]].interpolate(newMesh.vertices);
				glm::vec3 b = epm.vertices[epm_chain[i]].interpolate(newMesh.vertices);
				lengths.emplace_back(lengths.back() + glm::length(b - a));
			}
			if (!(initial_length == initial_length)) initial_length = lengths.back();

			std::unordered_map< uint32_t, uint32_t > visited;
			visited.reserve(epm.vertices.size());

			uint32_t first_i = (epm_chain[0] == epm_chain.back() ? 1 : 0);
			for (uint32_t i = first_i; i < epm_chain.size(); ++i) {
				auto ret = visited.insert(std::make_pair(epm_chain[i], i));
				if (ret.second) continue;

				uint32_t first = ret.first->second;
				uint32_t last = i;
				assert(first < last);
				assert(epm_chain[first] == epm_chain[last]);

				float length = lengths[last] - lengths[first];
				float length_outer = (lengths[first] - lengths[0]) + (lengths.back() - lengths[last]);


				if (epm_chain[0] == epm_chain.back() && length_outer < length) {
					//remove the 'outer' loop -- (last,back] + [0,first)
					verts_removed += first + (epm_chain.size() - (last + 1));
					length_removed += length_outer;

					remove_from_epm(0, first);
					remove_from_epm(last, epm_chain.size() - 1);
					epm_chain.erase(epm_chain.begin() + last + 1, epm_chain.end());
					epm_chain.erase(epm_chain.begin(), epm_chain.begin() + first);
					assert(epm_chain[0] == epm_chain.back()); //preserve circularity, right?
				}
				else {
					//remove the inner loop -- (first, last)
					verts_removed += (last + 1) - (first + 1);
					length_removed += length;

					remove_from_epm(first, last);
					epm_chain.erase(epm_chain.begin() + first + 1, epm_chain.begin() + last + 1);
				}
				again = true;
				break;
			}
		}

		if (verts_removed) {
			qDebug() << "Removed " << verts_removed << " vertices (that's " << length_removed << " units; " << length_removed / initial_length * 100.0 << "% of the initial length of " << initial_length << " units).";
		}

	};



	for (auto& epm_chain : left_of_epm) {
		cleanup_chain(epm_chain, 1);
	}

	for (auto& epm_chain : right_of_epm) {
		cleanup_chain(epm_chain, -(1 << 8));
	}

	//build split mesh:
	std::vector< EmbeddedVertex > split_verts;
	std::vector< glm::uvec3 > split_tris;
	std::vector< uint32_t > epm_to_split;
	epm.split_triangles(newMesh.vertices, newMesh.triangles, &split_verts, &split_tris, &epm_to_split);


	//transfer edge values to split mesh:
	std::unordered_map< glm::uvec2, int32_t > edge_values;
	for (auto const& se : epm.simplex_edges) {
		for (auto const& ee : se.second) {
			uint32_t a = epm_to_split[ee.first];
			uint32_t b = epm_to_split[ee.second];
			int32_t value = ee.value.sum;
			auto ret = edge_values.insert(std::make_pair(glm::uvec2(a, b), value));
			assert(ret.second);

			ret = edge_values.insert(std::make_pair(glm::uvec2(b, a), -value));
			assert(ret.second);
		}
	}

	//tag triangles with values:
	std::unordered_map< glm::uvec2, uint32_t > edge_to_tri;
	edge_to_tri.reserve(split_tris.size() * 3);
	for (auto const& tri : split_tris) {
		uint32_t ti = &tri - &split_tris[0];
		auto ret = edge_to_tri.insert(std::make_pair(glm::uvec2(tri.x, tri.y), ti)); assert(ret.second);
		ret = edge_to_tri.insert(std::make_pair(glm::uvec2(tri.y, tri.z), ti)); assert(ret.second);
		ret = edge_to_tri.insert(std::make_pair(glm::uvec2(tri.z, tri.x), ti)); assert(ret.second);
	}

	constexpr int32_t const Unvisited = std::numeric_limits< int32_t >::max();
	std::vector< int32_t > values(split_tris.size(), Unvisited);

	std::vector< bool > keep(split_tris.size(), false);

	for (uint32_t seed = 0; seed < split_tris.size(); ++seed) {
		if (values[seed] != Unvisited) continue;

		std::vector< uint32_t > component;
		component.reserve(split_tris.size());

		values[seed] = (128 << 8) | (128);
		component.emplace_back(seed);

		for (uint32_t ci = 0; ci < component.size(); ++ci) {
			uint32_t ti = component[ci];
			uint32_t value = values[ti];
			assert(value != Unvisited);
			glm::uvec3 tri = split_tris[ti];
			auto over = [&](uint32_t a, uint32_t b) {
				auto f = edge_to_tri.find(glm::uvec2(b,a));
				if (f == edge_to_tri.end()) return;
				int32_t nv = value;
				auto f2 = edge_values.find(glm::uvec2(b,a));
				if (f2 != edge_values.end()) nv += f2->second;

				if (values[f->second] == Unvisited) {
					values[f->second] = nv;
					component.emplace_back(f->second);
				} else {
					assert(values[f->second] == nv);
				}
			};
			over(tri.x, tri.y);
			over(tri.y, tri.z);
			over(tri.z, tri.x);
		}

		int32_t max_right = 128;
		int32_t max_left = 128;
		for (auto ti : component) {
			int32_t right = values[ti] >> 8;
			int32_t left = values[ti] & 0xff;
			max_right = std::max(max_right, right);
			max_left = std::max(max_left, left);
		}

		int32_t keep_value = (max_right << 8) | max_left;
		for (auto ti : component) {
			if (values[ti] == keep_value) {
				keep[ti] = true;
			}
		}
	}

	std::vector< uint32_t > split_vert_to_clipped_vertex(split_verts.size(), -1U);
	auto use_vertex = [&](uint32_t v) {
		if (split_vert_to_clipped_vertex[v] == -1U) {
			assert(v < split_verts.size());
			split_vert_to_clipped_vertex[v] = clipped_vertices.size();
			clipped_vertices.emplace_back(split_verts[v]);
		}
		return split_vert_to_clipped_vertex[v];
	};

	for (auto const& tri : split_tris) {
		if (!keep[&tri - &split_tris[0]]) continue;
		clipped.triangles.emplace_back(
			use_vertex(tri.x),
			use_vertex(tri.y),
			use_vertex(tri.z)
		);
	}
	clipped.vertices.reserve(clipped_vertices.size());
	for (auto const& v : clipped_vertices) {
		clipped.vertices.emplace_back(v.interpolate(newMesh.vertices));
	}

	qDebug() << "Trimmed model from " << newMesh.triangles.size() << " triangles on " << newMesh.vertices.size() << " vertices to " << clipped.triangles.size() << " triangles on " << clipped.vertices.size() << " vertices.";
	
	auto transform_chain = [&](std::vector< uint32_t > const& epm_chain) {
		assert(!epm_chain.empty());
		std::vector< uint32_t > split_chain;
		split_chain.reserve(epm_chain.size());
		for (auto v : epm_chain) {
			assert(v < epm_to_split.size());
			v = epm_to_split[v];
			assert(v != -1U);
			assert(v < split_vert_to_clipped_vertex.size());
			v = split_vert_to_clipped_vertex[v];
			//assert(v != -1U); //<-- sometimes chains don't include the edge loops (like when they cross)
			split_chain.emplace_back(v);
		}
		assert(!split_chain.empty());
		assert((epm_chain[0] == epm_chain.back()) == (split_chain[0] == split_chain.back()));
		return split_chain;
	};
	if (left_of_vertices_) {
		left_of_vertices_->clear();
		left_of_vertices_->reserve(left_of_epm.size());
		for (auto const& chain : left_of_epm) {
			left_of_vertices_->emplace_back(transform_chain(chain));
		}
	}
	if (right_of_vertices_) {
		right_of_vertices_->clear();
		right_of_vertices_->reserve(right_of_epm.size());
		for (auto const& chain : right_of_epm) {
			right_of_vertices_->emplace_back(transform_chain(chain));
		}
	}
}


/*
*	Function: peel a new portion of the mesh that represents the next step of KnitGraph creation
*		in essence, starting from the current active_chains, work within the length constraint to find the ending active_chains
*		and create a model that only contains vertices and triangles within this boundary
*	Return: void, but generates the sliced model and relevant information which is passed through the arguments
*	Called: called when the user selects the 'step' button in the GUI, is the second step after finding the active chains
*/
void KnitGrapher::peelSlice(std::vector< std::vector< EmbeddedVertex > > & active_chains,
	AutoKnitMesh* slice_,
	std::vector< EmbeddedVertex >* slice_on_model_,
	std::vector< std::vector< uint32_t > >* slice_active_chains_,
	std::vector< std::vector< uint32_t > >* slice_next_chains_,
	std::vector< bool >* used_boundary_) 
{

	//hope to god it works now;
	assert(slice_);
	auto& slice = *slice_;
	slice.clear();

	assert(slice_on_model_);
	auto& slice_on_model = *slice_on_model_;
	slice_on_model.clear();

	assert(slice_active_chains_);
	auto& slice_active_chains = *slice_active_chains_;
	slice_active_chains.clear();

	assert(slice_next_chains_);
	auto& slice_next_chains = *slice_next_chains_;
	slice_next_chains.clear();

	{
		uint32_t loops = 0;
		uint32_t lines = 0;
		for (auto const& chain : active_chains) {
			if (chain[0] == chain.back()) ++loops;
			else ++lines;
		}
		qDebug() << "---- peel slice on [" << loops << " loops and " << lines << " lines] ----";
	}

	AutoKnitMesh clipped;
	std::vector< EmbeddedVertex > clipped_on_model;
	std::vector< std::vector< EmbeddedVertex > > dummy;


	trimModel(active_chains, dummy, &clipped, &clipped_on_model, nullptr, nullptr);


	//This version of the code just uses the 3D distance to the curve.
	//might have problems with models that get really close to themselves.
	std::vector< float > values(clipped.vertices.size(), std::numeric_limits< float >::infinity());

	auto do_seg = [&values, &clipped](glm::vec3 const& a, glm::vec3 const& b) {
		if (a == b) return;
		glm::vec3 ab = b - a;
		float limit = glm::dot(ab, ab);
		float inv_limit = 1.0f / limit;
		for (auto const& v : clipped.vertices) {
			float amt = glm::dot(v - a, ab);
			amt = std::max(0.0f, std::min(limit, amt));
			glm::vec3 pt = (amt * inv_limit) * (b - a) + a;
			float dis2 = glm::length2(v - pt);
			float& best2 = values[&v - &clipped.vertices[0]];
			best2 = std::min(best2, dis2);
		}
	};

	for (auto const& chain : active_chains) {
		for (uint32_t i = 0; i + 1 < chain.size(); ++i) {
			do_seg(chain[i].interpolate(newMesh.vertices), chain[i + 1].interpolate(newMesh.vertices));
		}
	}

	for (auto& v : values) {
		v = std::sqrt(v);
	}

	
	std::vector< std::vector< EmbeddedVertex > > next_chains;
	{
		float level = 2.0f * stitchHeight / modelUnitLength;

		std::vector< std::vector< EmbeddedVertex > > level_chains;

		extractLevelChains(clipped, values, level, &level_chains);

		{ //(sort-of) hack: make all chains into loops by including portions of the boundary if needed:
			std::unordered_map< glm::uvec2, uint32_t > next;
			auto do_edge = [&next](uint32_t a, uint32_t b, uint32_t c) {
				auto ret = next.insert(std::make_pair(glm::uvec2(a, b), c));
				assert(ret.second);
			};
			for (auto const& tri : clipped.triangles) {
				do_edge(tri.x, tri.y, tri.z);
				do_edge(tri.y, tri.z, tri.x);
				do_edge(tri.z, tri.x, tri.y);
			}
			struct ChainEnd {
				ChainEnd(float along_, uint32_t chain_, bool is_start_) : along(along_), chain(chain_), is_start(is_start_) { }
				float along;
				uint32_t chain;
				bool is_start;
			};
			std::unordered_map< glm::uvec2, std::vector< ChainEnd > > on_edge;
			auto do_end = [&](EmbeddedVertex const& ev, uint32_t chain, bool is_start) {
				assert(ev.simplex.x != -1U);
				assert(ev.simplex.y != -1U);
				assert(ev.simplex.z == -1U);

				glm::uvec2 e = glm::uvec2(ev.simplex.x, ev.simplex.y);
				float amt = ev.weights.y;
				if (next.count(e)) {
					e = glm::uvec2(ev.simplex.y, ev.simplex.x);
					amt = ev.weights.x;
				}

				assert(next.count(e) + next.count(glm::uvec2(e.y, e.x)) == 1); //e should be a boundary edge
				assert(next.count(e) == 0);

				on_edge[e].emplace_back(amt, chain, is_start);
			};
			for (auto& chain : level_chains) {
				if (chain[0] == chain.back()) continue; //ignore loops
				do_end(chain[0], &chain - &level_chains[0], true);
				do_end(chain.back(), &chain - &level_chains[0], false);
			}
			std::vector< uint32_t > append(level_chains.size(), -1U);
			std::vector< bool > used_boundary(level_chains.size(), true);

			//loops marked as such before connection-making:
			for (uint32_t c = 0; c < level_chains.size(); ++c) {
				if (level_chains[c][0] == level_chains[c].back()) {
					append[c] = c;
					used_boundary[c] = false;
				}
			}

			auto chase_path = [&](glm::uvec2 begin_e, ChainEnd begin_ce) {
				assert(!begin_ce.is_start); //start at the end of a path, chase to the start of another
				std::vector< EmbeddedVertex > path;

				ChainEnd const* end_ce = nullptr;

				glm::uvec2 e = begin_e;
				ChainEnd ce = begin_ce;
				while (true) {
					{ //find next point along edge, if it exists:
						ChainEnd const* found_ce = nullptr;
						auto f = on_edge.find(e);
						if (f != on_edge.end()) {
							for (auto const& nce : f->second) {
								if (nce.along <= ce.along) continue;
								if (found_ce == nullptr || nce.along < found_ce->along) {
									found_ce = &nce;
								}
							}
						}
						if (found_ce) {
							assert(found_ce->is_start);
							end_ce = found_ce;
							break;
						}
					}
					//next point is end of edge.
					//add vertex:
					path.emplace_back(EmbeddedVertex::on_vertex(e.y));
					//circulate to next edge:
					glm::uvec2 old_e = e;
					while (true) {
						auto f = next.find(glm::uvec2(e.y, e.x));
						if (f == next.end()) break;
						e = glm::uvec2(f->second, e.y);
					}
					assert(e.y == old_e.y);
					assert(e.x != old_e.x);
					e = glm::uvec2(e.y, e.x);
					ce.chain = -1U;
					ce.along = 0.0f;
				}
				assert(end_ce);
				assert(end_ce->is_start); //start of (another?) path.
				path.emplace_back(level_chains[end_ce->chain][0]);
				level_chains[begin_ce.chain].insert(level_chains[begin_ce.chain].end(), path.begin(), path.end());
				//flag that append code should append other chain:
				assert(append[begin_ce.chain] == -1U);
				append[begin_ce.chain] = end_ce->chain;
			};

			for (auto const& seed_ece : on_edge) {
				for (auto const& seed : seed_ece.second) {
					if (!seed.is_start) {
						chase_path(seed_ece.first, seed);
					}
				}
			}

			for (uint32_t c = 0; c < level_chains.size(); ++c) {
				if (append[c] == -1U) continue; //marked for discard
				while (append[c] != c) {
					uint32_t a = append[c];
					assert(a < level_chains.size());
					assert(!level_chains[a].empty());
					assert(level_chains[c].back() == level_chains[a][0]); //already have first vertex
					level_chains[c].insert(level_chains[c].end(), level_chains[a].begin() + 1, level_chains[a].end());
					append[c] = append[a];

					append[a] = -1U; //mark for discard
					level_chains[a].clear();
				}
			}
			for (uint32_t c = 0; c < level_chains.size(); /* later */) {
				if (append[c] == -1U) {
					assert(level_chains[c].empty());
					used_boundary[c] = used_boundary.back();
					used_boundary.pop_back();
					level_chains[c] = level_chains.back();
					level_chains.pop_back();
				}
				else {
					assert(!level_chains[c].empty());
					assert(level_chains[c][0] == level_chains[c].back());
					++c;
				}
			}

			if (used_boundary_) *used_boundary_ = used_boundary;

		}

		uint32_t loops = 0;
		uint32_t lines = 0;

		next_chains.reserve(level_chains.size());
		for (auto& chain : level_chains) {
			//chain is embedded on 'clipped' which is embedded on 'model'; re-embed on just 'model':
			for (auto& v : chain) {
				glm::uvec3 simplex = clipped_on_model[v.simplex.x].simplex;
				if (v.simplex.y != -1U) simplex = EmbeddedVertex::common_simplex(simplex, clipped_on_model[v.simplex.y].simplex);
				if (v.simplex.z != -1U) simplex = EmbeddedVertex::common_simplex(simplex, clipped_on_model[v.simplex.z].simplex);
				glm::vec3 weights = v.weights.x * clipped_on_model[v.simplex.x].weights_on(simplex);
				if (v.simplex.y != -1U) weights += v.weights.y * clipped_on_model[v.simplex.y].weights_on(simplex);
				if (v.simplex.z != -1U) weights += v.weights.z * clipped_on_model[v.simplex.z].weights_on(simplex);
				v.simplex = simplex;
				v.weights = weights;
			}

			//subdivide chain and add to outputs:
			if (chain[0] == chain.back()) ++loops;
			else ++lines;
			next_chains.emplace_back();
			sampleChain(getChainSampleSpacing(),  chain, &next_chains.back());
		}
		qDebug() << "  extracted " << loops << " loops and " << lines << " lines.";
	}

	for (auto const& chain : active_chains) {
		for (auto const& v : chain) {
			assert(v.simplex.x < newMesh.vertices.size());
			assert(v.simplex.y == -1U || v.simplex.y < newMesh.vertices.size());
			assert(v.simplex.z == -1U || v.simplex.z < newMesh.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i - 1] != chain[i]);
		}
	}

	for (auto const& chain : next_chains) {
		for (auto const& v : chain) {
			assert(v.simplex.x < newMesh.vertices.size());
			assert(v.simplex.y == -1U || v.simplex.y < newMesh.vertices.size());
			assert(v.simplex.z == -1U || v.simplex.z < newMesh.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i - 1] != chain[i]);
		}
	}

	trimModel(active_chains, next_chains, &slice, &slice_on_model, &slice_active_chains, &slice_next_chains);

	//sometimes this can combine vertices, in which case the output chains should be trimmed:
	uint32_t trimmed = 0;
	for (auto& chain : slice_active_chains) {
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); /* later */) {
			if (chain[i - 1] == chain[i]) {
				chain.erase(chain.begin() + i);
				++trimmed;
			}
			else {
				++i;
			}
		}
	}

	for (auto& chain : slice_next_chains) {
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			if (chain[i - 1] == chain[i]) {
				chain.erase(chain.begin() + i);
				++trimmed;
			}
			else {
				++i;
			}
		}
	}

	if (trimmed) {
		qDebug() << "Trimmed " << trimmed << " too-close-for-epm vertices from slice chains.";
	}
}




/*
*	Function: fill in unassigned elements in closest vertex vector based on the weights provided
*	Return: bool, has the function assigned the values or not
*	Called: called during the linkChains() function
*/
bool KnitGrapher::fillUnassigned(std::vector< uint32_t >& closest, std::vector< float > const& weights, bool is_loop) 
{

	bool have_assigned = false;
	for (auto c : closest) {
		if (c != -1U) {
			have_assigned = true;
			break;
		}
	}
	if (!have_assigned) return false;

	auto do_range = [&](uint32_t first, uint32_t last) {
		uint32_t before = closest[first > 0 ? first - 1 : closest.size() - 1];
		uint32_t after = closest[last + 1 < closest.size() ? last + 1 : 0];
		if (!is_loop) {
			if (first == 0) {
				assert(last + 1 < closest.size());
				before = after;
			}
			if (last + 1 == closest.size()) {
				assert(first > 0);
				after = before;
			}
		}
		assert(closest[first] == -1U);
		assert(closest[last] == -1U);
		assert(before != -1U);
		assert(after != -1U);

		//---------- fill in approximately equal weight sums -----------

		float total = 0.0f;
		{ //first pass: figure out total weight sum:
			uint32_t i = first;
			while (true) {
				assert(closest[i] == -1U);
				total += weights[i];
				if (i == last) break;
				++i;
			}
		}
		float sum = 0.0f;
		{ //second pass: assign before/after based on weight sum:
			uint32_t i = first;
			while (true) {
				if (sum + 0.5f * weights[i] < 0.5f * total) {
					closest[i] = before;
				}
				else {
					closest[i] = after;
				}
				sum += weights[i];
				if (i == last) break;
				++i;
			}
		}

	};

	for (uint32_t seed = 0; seed < closest.size(); ++seed) {
		if (closest[seed] != -1U) continue;
		uint32_t first = seed;
		while ((first > 0 || is_loop) && closest[first > 0 ? first - 1 : closest.size() - 1] == -1U) {
			first = (first > 0 ? first - 1 : closest.size() - 1);
		}
		uint32_t last = seed;
		while ((last + 1 < closest.size() || is_loop) && closest[last + 1 < closest.size() ? last + 1 : 0] == -1U) {
			last = (last + 1 < closest.size() ? last + 1 : 0);
		}

		do_range(first, last);
	}

	for (auto c : closest) {
		assert(c != -1U);
	}

	return true;
}


/*
*	Function: find the links between vertices such that the total distance between the points is minimized
*		in essence is linking the different vectors with a 'knit'
*	Return: void, but modifies the passed 'links_' argument, which is the resulting links between vertices
*	Called: called during the linkChains() function 
*/
void KnitGrapher::optimalLink(
	float target_distance, bool do_roll,
	std::vector< glm::vec3 > const& source,
	std::vector< bool > const& source_linkone,
	std::vector< glm::vec3 > const& target,
	std::vector< bool > const& target_linkone,
	std::vector< std::pair< uint32_t, uint32_t > >* links_) 
{

	assert(source.size() == source_linkone.size());
	assert(target.size() == target_linkone.size());

	assert(links_);
	auto& links = *links_;
	links.clear();

	//options:
	// s[i]  s[i] s[i+1]   s[i]
	//  |       \ /       /  |  
	// t[i]     t[i]    t[i] t[i+1]

	struct State {
		uint16_t source_idx = 0xffff;
		uint16_t target_idx = 0xffff;
		uint16_t source_remain = 0x0; //unlinked source
		uint16_t target_remain = 0x0; //unlinked target
		uint64_t pack() {
			uint64_t ret;
			memcpy(reinterpret_cast<char*>(&ret), this, sizeof(uint64_t));
			return ret;
		}
		static State unpack(uint64_t packed) {
			State ret;
			memcpy(reinterpret_cast<char*>(&ret), &packed, sizeof(uint64_t));
			return ret;
		}
	};
	static_assert(sizeof(State) == sizeof(uint64_t), "State should be packed");

	struct Action {
		Action(uint8_t take_source_, uint8_t take_target_) : take_source(take_source_), take_target(take_target_) { }
		Action() = default;
		uint8_t take_source = 0;
		uint8_t take_target = 0;
	};

	std::unordered_map< uint64_t, std::pair< float, Action > > distance_via;
	std::vector< std::pair< float, uint64_t > > heap;

	auto visit = [&distance_via, &heap](uint64_t state, float distance, Action const& action) {
		auto f = distance_via.insert(std::make_pair(state, std::make_pair(std::numeric_limits< float >::infinity(), Action()))).first;
		if (distance < f->second.first) {
			f->second = std::make_pair(distance, action);
			heap.emplace_back(-distance, state);
			std::push_heap(heap.begin(), heap.end());
		}
	};

	if (do_roll) {
		for (uint32_t t = 0; t < target.size(); ++t) {
			State s;
			s.source_idx = 0;
			s.target_idx = t;
			s.source_remain = source.size();
			s.target_remain = target.size();
			visit(s.pack(), 0.0f, Action());

			s.source_idx = 1;
			visit(s.pack(), 0.0f, Action());
		}
	}
	else {
		State s;
		s.source_idx = 0;
		s.target_idx = 0;
		s.source_remain = source.size();
		s.target_remain = target.size();
		visit(s.pack(), 0.0f, Action());
	}

	auto penalty = [&source, &target, &target_distance](uint32_t si, uint32_t ti) {
		assert(si < source.size());
		assert(ti < target.size());
		float dis = glm::length(source[si] - target[ti]) - target_distance;
		return dis * dis;
	};

	auto is_possible = [](State const& state) {
		if (state.source_remain * 2 < state.target_remain) return false;
		if (state.target_remain * 2 < state.source_remain) return false;
		return true;
	};

	State best;
	while (!heap.empty()) {
		std::pop_heap(heap.begin(), heap.end());
		State state = State::unpack(heap.back().second);
		float distance = -heap.back().first;
		heap.pop_back();
		auto f = distance_via.find(state.pack());
		assert(f != distance_via.end());
		assert(f->second.first <= distance);
		if (f->second.first < distance) continue;
		if (state.source_remain == 0 && state.target_remain == 0) {
			best = state;
			break;
		}
		assert(is_possible(state));
		//try the actions:
		{ // 1-1
			float next_distance = distance + penalty(state.source_idx, state.target_idx);
			State next = state;
			next.source_idx = (next.source_idx + 1) % source.size();
			next.target_idx = (next.target_idx + 1) % target.size();
			next.source_remain -= 1;
			next.target_remain -= 1;
			if (is_possible(next)) visit(next.pack(), next_distance, Action(1, 1));
		}
		// 1-2
		if (state.target_remain >= 2
			&& !source_linkone[state.source_idx]
			&& !target_linkone[state.target_idx]
			&& !target_linkone[(state.target_idx + 1) % target.size()]) {
			float next_distance = distance
				+ penalty(state.source_idx, state.target_idx)
				+ penalty(state.source_idx, (state.target_idx + 1) % target.size());
			State next = state;
			next.source_idx = (next.source_idx + 1) % source.size();
			next.target_idx = (next.target_idx + 2) % target.size();
			next.source_remain -= 1;
			next.target_remain -= 2;
			if (is_possible(next)) visit(next.pack(), next_distance, Action(1, 2));
		}
		// 2-1
		if (state.source_remain >= 2
			&& !source_linkone[state.source_idx]
			&& !source_linkone[(state.source_idx + 1) % source.size()]
			&& !target_linkone[state.target_idx]) {
			float next_distance = distance
				+ penalty(state.source_idx, state.target_idx)
				+ penalty((state.source_idx + 1) % source.size(), state.target_idx);
			State next = state;
			next.source_idx = (next.source_idx + 2) % source.size();
			next.target_idx = (next.target_idx + 1) % target.size();
			next.source_remain -= 2;
			next.target_remain -= 1;
			if (is_possible(next)) visit(next.pack(), next_distance, Action(2, 1));
		}
	}

	if (best.source_idx == 0xffff) {
		//Failed!
		throw std::runtime_error("Failed to link");
	}

	//Read back:
	State at = best;
	while (true) {
		auto f = distance_via.find(at.pack());
		assert(f != distance_via.end());
		Action action = f->second.second;
		if (action.take_source == 0 && action.take_target == 0) break;
		if (action.take_source == 1) {
			at.source_idx = (at.source_idx + source.size() - 1) % source.size();
			at.source_remain += 1;
			while (action.take_target > 0) {
				--action.take_target;
				at.target_idx = (at.target_idx + target.size() - 1) % target.size();
				at.target_remain += 1;
				links.emplace_back(at.source_idx, at.target_idx);
			}
		}
		else if (action.take_target == 1) {
			at.target_idx = (at.target_idx + target.size() - 1) % target.size();
			at.target_remain += 1;
			while (action.take_source > 0) {
				--action.take_source;
				at.source_idx = (at.source_idx + source.size() - 1) % source.size();
				at.source_remain += 1;
				links.emplace_back(at.source_idx, at.target_idx);
			}
		}
		else {
			assert(action.take_source == 1 || action.take_target == 1);
		}
	}
	assert(at.source_remain == source.size() && at.target_remain == target.size());
	assert(at.source_idx == 0 || at.source_idx == 1);

	std::reverse(links.begin(), links.end());

}


/*
*	Function: try to 'flatten' the mesh such that it can be knit on a V-bed machine, but preserving the original mesh structure
*	Return: void, but modifies the 'closest' argument
*	Called: called during the linkChains() function
*/
void KnitGrapher::flatten(std::vector< uint32_t >& closest, std::vector< float > const& weights, bool is_loop) 
{
	assert(closest.size() == weights.size());
	if (closest.empty()) return;

	//make sure that 'closest' looks like:
	//  a a a b b b c c c
	//   a a b b b b c c
	// that is, can be flattened to the knitting machine bed while preserving constituent cycles
	// One view of this: if you start at some stitch A, then the left side should be a mirror of the right side
	//  (with the possible exception that some symbols may be skipped on the left or right)

	//(a) condense closest into short list of symbols:
	std::vector< std::pair< uint32_t, float > > symbols;
	symbols.reserve(closest.size()); //certainly no more symbols than closest
	for (uint32_t i = 0; i < closest.size(); ++i) {
		uint32_t symb = closest[i];
		if (symbols.empty() || symbols.back().first != symb) {
			symbols.emplace_back(std::make_pair(symb, 0.0f));
		}
		assert(symbols.back().first == symb);
		symbols.back().second += weights[i];
	}
	assert(!symbols.empty());

	//(b) early-out in certain easy-to-check conditions:
	if (symbols.size() == 1) return; //single symbol

	//symbols -> bits
	std::vector< std::pair< uint16_t, float > > bit_symbols;
	uint32_t bits;
	{
		bit_symbols.reserve(symbols.size());
		std::map< uint32_t, uint16_t > symbol_bit;
		for (auto& sw : symbols) {
			auto ret = symbol_bit.insert(std::make_pair(sw.first, uint16_t(1 << symbol_bit.size())));
			assert(ret.first->second && "Only have enough bits for 16-symbol flattening");
			bit_symbols.emplace_back(ret.first->second, sw.second);
		}
		assert(bit_symbols.size() == symbols.size());
		bits = symbol_bit.size();
	}

	//helper struct inside the function, could technically be moved to the class definition but unecessary
	struct State {
		uint16_t used = 0; //have seen all symbols with bits in used
		uint8_t min = 0; //symbols that are strictly between min and max have been processed
		uint8_t max = 0;
		uint16_t current = 0; //most recently kept symbol
		uint16_t padding = 0; //padding to make state 64 bits long

		typedef uint64_t Packed;
		Packed pack() const {
			return *reinterpret_cast<Packed const*>(this);
		}
		static State unpack(Packed packed) {
			State ret;
			memcpy(reinterpret_cast<char*>(&ret), &packed, sizeof(State));
			return ret;
		}

		std::string to_string(uint32_t bits = 16) const {
			std::string ret = "(" + std::to_string(min) + "," + std::to_string(max) + ")";
			for (uint32_t i = bits - 1; i < bits; --i) {
				if (used & (1 << i)) {
					if ((1 << i) == current) {
						ret += "*";
					}
					else {
						ret += "x";
					}
				}
				else {
					ret += ".";
				}
			}
			return ret;
		}
	};
	static_assert(sizeof(State) == 8, "packed state");


	struct {
		State::Packed state;
		float cost = std::numeric_limits< float >::infinity();
		State::Packed from;
	} finished;
	std::unordered_map< State::Packed, std::pair< float, State::Packed > > visited;
	std::vector< std::pair< float, State::Packed > > todo;
	static std::greater< std::pair< float, State::Packed > > const TODOCompare;


	auto queue_state = [&visited, &finished, &todo, bits](State const state, float const cost, State const from) {
		assert(state.min != from.min || state.max != from.max); //must have done *something*

		(void)bits;

		if ((state.min != from.min && state.min == from.max)
			|| (state.max != from.max && state.max == from.min)) {
			//pointers crossed or met -> state is finished!
			assert(state.min == state.max || (state.min == from.max && state.max == from.min));
			if (cost < finished.cost) {
				finished.state = state.pack();
				finished.cost = cost;
				finished.from = from.pack();
			}
			return;
		}

		//queue/indicate regular
		auto ret = visited.insert(std::make_pair(state.pack(), std::make_pair(cost, from.pack())));
		if (ret.second || ret.first->second.first > cost) {
			ret.first->second = std::make_pair(cost, from.pack());
			todo.emplace_back(std::make_pair(cost, state.pack()));
			std::push_heap(todo.begin(), todo.end(), TODOCompare);
		}

	};

	auto expand_state = [&bit_symbols, &queue_state](State const state, float const cost) {
		auto const used = state.used;
		auto const min = state.min;
		auto const max = state.max;
		auto const current = state.current;
		assert(state.padding == 0);

		//*_next_symbol is the symbol that is advanced over when moving min/max,
		// leading to some asymmetry in indexing:
		// a(bc)d -> (abc)d <-- min_next_symbol is 'a' (at index of min_next)
		// a(bc)d -> a(bcd) <-- max_next_symbol is 'd' (at index of max)

		auto min_next = (min == 0 ? bit_symbols.size() - 1 : min - 1);
		auto min_next_symbol = bit_symbols[min_next];
		auto max_next = (max + 1U < bit_symbols.size() ? max + 1 : 0);
		auto max_next_symbol = bit_symbols[max];

		//actions:
		//no reason not to keep if symbol is current:
		if (min_next_symbol.first == current || max_next_symbol.first == current) {
			assert(used & current);
			State next;
			next.used = used;
			next.min = (min_next_symbol.first == current ? min_next : min);
			next.max = (max_next_symbol.first == current ? max_next : max);
			next.current = current;

			float next_cost = cost;

			queue_state(next, next_cost, state);

			return; //no other actions worth taking; this one was free!
		}

		//keep min (symbol must be unused):
		if (!(used & min_next_symbol.first)) {
			State next;
			next.used = used | min_next_symbol.first;
			next.min = min_next;
			next.max = max;
			next.current = min_next_symbol.first;

			float next_cost = cost;

			queue_state(next, next_cost, state);
		}

		//keep max (symbol must be unused):
		if (!(used & max_next_symbol.first)) {
			State next;
			next.used = used | max_next_symbol.first;
			next.min = min;
			next.max = max_next;
			next.current = max_next_symbol.first;

			float next_cost = cost;

			queue_state(next, next_cost, state);
		}

		//discard min:
		{
			State next;
			next.used = used;
			next.min = min_next;
			next.max = max;
			next.current = current;

			float next_cost = cost + min_next_symbol.second;

			queue_state(next, next_cost, state);
		}

		//discard max:
		{
			State next;
			next.used = used;
			next.min = min;
			next.max = max_next;
			next.current = current;

			float next_cost = cost + max_next_symbol.second;

			queue_state(next, next_cost, state);
		}
	};


	for (uint32_t s = 0; s < bit_symbols.size(); ++s) {
		State init;
		init.used = 0;
		init.min = s;
		init.max = s;
		init.current = 0;
		expand_state(init, 0.0f);
		if (!is_loop) break;
	}


	while (!todo.empty()) {
		std::pop_heap(todo.begin(), todo.end(), TODOCompare);
		auto state = State::unpack(todo.back().second);
		float cost = todo.back().first;
		todo.pop_back();
		//if the cheapest thing is more expensive than the finish, we're done:
		if (cost >= finished.cost) break;

		{ //cost should either be stale or what is stored in 'visited':
			auto f = visited.find(state.pack());
			assert(f != visited.end());
			if (cost > f->second.first) continue;
			assert(cost == f->second.first);
		}
		expand_state(state, cost);
	}
	assert(finished.cost != std::numeric_limits< float >::infinity()); //found ~some~ path

	
	std::vector< State::Packed > path;
	path.emplace_back(finished.state);
	path.emplace_back(finished.from);
	while (true) {
		auto f = visited.find(path.back());
		if (f == visited.end()) break;
		path.emplace_back(f->second.second);
	}
	std::reverse(path.begin(), path.end());
	//std::cout << "----" << std::endl; //DEBUG

	std::vector< int8_t > keep(bit_symbols.size(), -1);
	for (uint32_t i = 1; i < path.size(); ++i) {
		State state = State::unpack(path[i - 1]);
		State next = State::unpack(path[i]);

		if (state.min != next.min && state.max != next.max) {
			//a(bc)d -> (abcd), keep 'a' (next.min), 'd' (state.max)
			assert(bit_symbols[next.min].first == bit_symbols[state.max].first);
			assert(next.current == bit_symbols[state.min].first);
			assert(keep[next.min] == -1);
			assert(keep[state.max] == -1);
			keep[next.min] = keep[state.max] = 1;
		}
		else if (state.min != next.min) {
			//a(bc)d -> (abc)d, keep/discard next.min
			if (bit_symbols[next.min].first == next.current) {
				assert(keep[next.min] == -1);
				keep[next.min] = 1;
			}
			else {
				assert(keep[next.min] == -1);
				keep[next.min] = 0;
			}
		}
		else {
			assert(state.max != next.max);
			//a(bc)d -> a(bcd), keep/discard state.max
			if (bit_symbols[state.max].first == next.current) {
				assert(keep[state.max] == -1);
				keep[state.max] = 1;
			}
			else {
				assert(keep[state.max] == -1);
				keep[state.max] = 0;
			}
		}
	}

	//DEBUG: was at least one thing kept?
	bool kept_at_least_one = false;
	for (auto k : keep) {
		assert(k != -1);
		if (k == 1) kept_at_least_one = true;
	}
	assert(kept_at_least_one);


	//use keep to figure out which elements of closest to re-label.
	std::vector< bool > relabel; relabel.reserve(closest.size());
	{
		//std::cout << "relabel:"; //DEBUG
		auto si = symbols.begin();
		for (auto c : closest) {
			assert(si != symbols.end());
			if (si->first != c) ++si;
			assert(si != symbols.end());
			assert(si->first == c);
			relabel.emplace_back(keep[si - symbols.begin()] == 0);
		}
		assert(relabel.size() == closest.size());
		assert(si != symbols.end());
		++si;
		assert(si == symbols.end());

		bool have_keep = false;
		for (uint32_t seed = 0; seed < closest.size(); ++seed) {
			if (relabel[seed] == false) have_keep = true;
		}
		assert(have_keep);
	}

	{ //do relabelling:
		auto relabel_range = [&closest, &weights, &relabel](uint32_t first, uint32_t last) {
			uint32_t before = (first == 0 ? closest.back() : closest[first - 1]);
			uint32_t after = (last + 1 == closest.size() ? closest[0] : closest[last + 1]);

			assert(!relabel[(first == 0 ? closest.size() : first) - 1]);
			assert(!relabel[(last + 1 == closest.size() ? 0 : last) + 1]);
			assert(relabel[first]);
			assert(relabel[last]);

			float total = 0.0f;
			{ //first pass: figure out total weight sum:
				uint32_t i = first;
				while (true) {
					assert(relabel[i]);
					total += weights[i];
					if (i == last) break;
					++i;
				}
			}
			float sum = 0.0f;
			{ //second pass: assign before/after based on weight sum:
				uint32_t i = first;
				while (true) {
					if (sum + 0.5f * weights[i] < 0.5f * total) {
						closest[i] = before;
					}
					else {
						closest[i] = after;
					}
					relabel[i] = false;
					sum += weights[i];
					if (i == last) break;
					++i;
				}
			}
		};
		for (uint32_t seed = 0; seed < closest.size(); ++seed) {
			if (!relabel[seed]) continue;
			uint32_t first = seed;
			while (relabel[first > 0 ? first - 1 : closest.size() - 1]) {
				first = (first > 0 ? first - 1 : closest.size() - 1);
			}
			uint32_t last = seed;
			while (relabel[last + 1 < closest.size() ? last + 1 : 0]) {
				last = (last + 1 < closest.size() ? last + 1 : 0);
			}
			relabel_range(first, last);
		}
	}
}



//huge function from autoknit, below is the copied comments
//NOTE:
// in the output of link_chains, links are ordered so that if a vertex appears
// in more than one link, the links in which it appears are adjacent, and the
// order is increasing along the opposite chain.
// e.g.
//  -------a--->
//       / |
//  ----b--c--->
//  the links array will contain
//        ... (a,b) (a,c) ...
//     or ... (b,a) (c,a) ...
//  (depending on which chain is 'from' and which is 'to')
/*
*	Function: link the 2 chains (current and next) with new stitches
*		in essence defines how the knitting operations should continue to get from current row to the next row
*	Return: void, but modifies a sleu of passed arguments
*	Called: called when the user selects the 'step' button on the GUI, is the 3rd step after peelSlice()
*/
void KnitGrapher::linkChains(
	AutoKnitMesh const& slice, //in: slice on which the chains reside
	std::vector< float > const& slice_times, //in: time field (times @ vertices), for slice
	std::vector< std::vector< uint32_t > > const& active_chains, //in: current active chains (slice vertex #'s)
	std::vector< std::vector< Stitch > > const& active_stitches, //in: current active stitches, sorted by time
	std::vector< std::vector< uint32_t > > const& next_chains, //in: current next chains (slice vertex #'s)
	std::vector< bool > const& next_used_boundary, //in: did next chain use boundary? (forces no discard)
	//need this or slice_times (above) std::vector< std::vector< bool > > const &discard_segments,
	std::vector< std::vector< Stitch > >* next_stitches_, //out: next active stitches
	std::vector< Link >* links_ //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
) 
{


	assert(slice_times.size() == slice.vertices.size());

	for (auto const& chain : active_chains) {
		assert(chain.size() >= 2);
		assert(chain[0] != chain.back() || chain.size() >= 3);
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
	}
	assert(active_stitches.size() == active_chains.size());

	for (auto const& stitches : active_stitches) {
		for (auto const& s : stitches) {
			assert(s.t >= 0.0f);
			assert(s.t < 1.0f);
			if (&s > &stitches[0]) {
				assert(s.t > (&s - 1)->t);
			}
		}
	}

	for (auto const& chain : next_chains) {
		assert(chain.size() >= 2);
		assert(chain[0] != chain.back() || chain.size() >= 3);
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
	}

	assert(next_used_boundary.size() == next_chains.size());

	assert(next_stitches_);
	auto& next_stitches = *next_stitches_;
	next_stitches.clear();

	assert(links_);
	auto& links = *links_;
	links.clear();

	//figure out the time to trim after:
	float active_max_time = -std::numeric_limits< float >::infinity();
	for (auto const& chain : active_chains) {
		for (auto v : chain) {
			active_max_time = std::max(active_max_time, slice_times[v]);
		}
	}

	//first, write down length to each vertex (makes it easy to move to/from parameter space):
	auto make_lengths = [&slice](std::vector< std::vector< uint32_t > > const& chains) {
		std::vector< std::vector< float > > all_lengths;
		all_lengths.reserve(chains.size());
		for (auto const& chain : chains) {
			all_lengths.emplace_back();
			std::vector< float >& lengths = all_lengths.back();
			lengths.reserve(chain.size());
			float total_length = 0.0f;
			lengths.emplace_back(total_length);
			for (uint32_t vi = 1; vi < chain.size(); ++vi) {
				glm::vec3 const& a = slice.vertices[chain[vi - 1]];
				glm::vec3 const& b = slice.vertices[chain[vi]];
				total_length += glm::length(b - a);
				lengths.emplace_back(total_length);
			}
			assert(lengths.size() == chain.size());
		}
		return all_lengths;
	};
	std::vector< std::vector< float > > active_lengths = make_lengths(active_chains);
	std::vector< std::vector< float > > next_lengths = make_lengths(next_chains);

	
	std::vector< std::vector< std::pair< float, bool > > > next_discard_after;
	next_discard_after.reserve(next_chains.size());
	for (auto const& chain : next_chains) {
		uint32_t ci = &chain - &next_chains[0];
		std::vector< float > const& lengths = next_lengths[ci];
		assert(lengths.size() == chain.size());
		next_discard_after.emplace_back();
		auto& discard_after = next_discard_after.back();
		discard_after.emplace_back(0.0f, (slice_times[chain[0]] > active_max_time));
		for (uint32_t vi = 0; vi + 1 < chain.size(); ++vi) {
			float ta = slice_times[chain[vi]];
			float tb = slice_times[chain[vi + 1]];
			float la = lengths[vi];
			float lb = lengths[vi + 1];
			//note: for the purposes of assigning ranges, treat time t == active_max_time as t - epsilon.
			if ((ta <= active_max_time && tb > active_max_time)
				|| (ta > active_max_time && tb <= active_max_time)) {
				float m = (active_max_time - ta) / (tb - ta);
				float l = m * (lb - la) + la;
				discard_after.emplace_back(l / lengths.back(), (tb > active_max_time));
			}
			else {
				assert((ta > active_max_time) == (tb > active_max_time));
				//do nothing!
			}
		}

		bool is_loop = (chain[0] == chain.back());

		//remove any non-discard segments shorter than 1.5 stitches:
		if (discard_after.size() > 1) {
			float const MinLength = 1.5f * stitchWidth / modelUnitLength;

			float first_len = discard_after[1].first - 0.0f;
			float last_len = 1.0f - discard_after.back().first;
			if (is_loop) {
				first_len = last_len = first_len + last_len;
			}

			for (uint32_t s = 0; s < discard_after.size(); ++s) {
				if (discard_after[s].second != false) continue;
				float len;
				if (s == 0) len = first_len;
				else if (s + 1 == discard_after.size()) len = last_len;
				else len = discard_after[s + 1].first - discard_after[s].first;
				if (len * lengths.back() < MinLength) {
					discard_after[s].second = true;
				}
			}
			if (is_loop) {
				assert(discard_after[0].second == discard_after.back().second);
			}
			for (uint32_t s = 0; s + 1 < discard_after.size(); /* later */) {
				if (discard_after[s].second == discard_after[s + 1].second) {
					discard_after.erase(discard_after.begin() + (s + 1));
				}
				else {
					++s;
				}
			}
		}

		//remove any discard segments shorter than 0.5 stitches:
		if (discard_after.size() > 1) {
			float const MinLength = 0.5f * stitchWidth / modelUnitLength;

			float first_len = discard_after[1].first - 0.0f;
			float last_len = 1.0f - discard_after.back().first;
			if (is_loop) {
				first_len = last_len = first_len + last_len;
			}

			for (uint32_t s = 0; s < discard_after.size(); ++s) {
				if (discard_after[s].second != true) continue;
				float len;
				if (s == 0) len = first_len;
				else if (s + 1 == discard_after.size()) len = last_len;
				else len = discard_after[s + 1].first - discard_after[s].first;
				if (len * lengths.back() < MinLength) {
					discard_after[s].second = false;
				}
			}
			if (is_loop) {
				assert(discard_after[0].second == discard_after.back().second);
			}
			for (uint32_t s = 0; s + 1 < discard_after.size(); /* later */) {
				if (discard_after[s].second == discard_after[s + 1].second) {
					discard_after.erase(discard_after.begin() + (s + 1));
				}
				else {
					++s;
				}
			}
		}
	}

	{ //for any next chains that touch a boundary, mark as 'accept':
		uint32_t marked = 0;
		for (uint32_t ni = 0; ni < next_chains.size(); ++ni) {
			if (next_used_boundary[ni] == false) continue;
			if (next_discard_after[ni].size() != 1 || next_discard_after[ni][0] != std::make_pair(0.0f, false)) {
				next_discard_after[ni].assign(1, std::make_pair(0.0f, false));
				++marked;
			}
		}
		if (marked) {
			qDebug() << "NOTE: marked " << marked << " next chains as all-accept because they touch boundaries.";
		}
	}

	{ //if all segments are marked 'discard', then mark everything 'accept':
		bool only_discard = true;
		for (auto const& discard_after : next_discard_after) {
			for (auto const& d : discard_after) {
				if (!d.second) {
					only_discard = false;
					break;
				}
			}
			if (only_discard == false) break;
		}

		if (only_discard) {
			qDebug() << "Marking everything accept because it was all marked discard.";
			for (auto& discard_after : next_discard_after) {
				assert(discard_after.size() == 1);
				assert(discard_after[0] == std::make_pair(0.0f, true));
				discard_after[0].second = false;
			}
		}
		else {
			qDebug() << "Have a mix of discard and accept.";
		}
	}



	std::vector< std::vector< uint32_t > > active_closest;
	std::vector< std::vector< uint32_t > > next_closest;

	
	{ //find closest pairs:

		std::vector< std::vector< uint32_t > > adj(slice.vertices.size());
		{ //adjacency matrix -- always handy:
			std::unordered_set< glm::uvec2 > edges;
			for (auto const& tri : slice.triangles) {
				auto do_edge = [&edges](uint32_t a, uint32_t b) {
					if (b > a) std::swap(a, b);
					edges.insert(glm::uvec2(a, b));
				};
				do_edge(tri.x, tri.y);
				do_edge(tri.y, tri.z);
				do_edge(tri.z, tri.x);
			}
			for (auto const& e : edges) {
				adj[e.x].emplace_back(e.y);
				adj[e.y].emplace_back(e.x);
			}
		}

		auto closest_source_chain = [&slice, &adj](
			std::vector< std::vector< uint32_t > > const& sources,
			std::vector< std::vector< uint32_t > > const& targets) {

				std::vector< float > dis(slice.vertices.size(), std::numeric_limits< float >::infinity());
				std::vector< uint32_t > from(slice.vertices.size(), -1U);
				std::vector< std::pair< float, uint32_t > > todo;

				auto queue = [&dis, &todo, &from](uint32_t n, float d, uint32_t f) {
					assert(d < dis[n]);
					dis[n] = d;
					from[n] = f;
					todo.emplace_back(std::make_pair(-d, n));
					std::push_heap(todo.begin(), todo.end());
				};

				for (auto const& chain : sources) {
					uint32_t ci = &chain - &sources[0];
					for (auto const& v : chain) {
						if (v == -1U) continue;
						if (0.0f < dis[v]) queue(v, 0.0f, ci); //some verts appear twice
					}
				}

				while (!todo.empty()) {
					std::pop_heap(todo.begin(), todo.end());
					uint32_t at = todo.back().second;
					float d = -todo.back().first;
					todo.pop_back();
					if (d > dis[at]) continue;
					assert(d == dis[at]);
					for (auto n : adj[at]) {
						float nd = d + glm::length(slice.vertices[n] - slice.vertices[at]);
						if (nd < dis[n]) queue(n, nd, from[at]);
					}
				}

				std::vector< std::vector< uint32_t > > closest;
				closest.reserve(targets.size());
				for (auto const& chain : targets) {
					closest.emplace_back();
					closest.back().reserve(targets.size());
					for (auto v : chain) {
						if (v == -1U) {
							closest.back().emplace_back(-1U);
						}
						else {
							closest.back().emplace_back(from[v]);
						}
					}
				}
				return closest;
		};

		active_closest = closest_source_chain(next_chains, active_chains);
		next_closest = closest_source_chain(active_chains, next_chains);

		//HACK: treat these as 'segment' values instead of 'vertex' (but really every segment just gets first vertex's value):
		for (auto& c : active_closest) {
			if (!c.empty()) c.pop_back();
		}
		for (auto& c : next_closest) {
			if (!c.empty()) c.pop_back();
		}
	}

	auto discard_nonmutual = [&]() {
		std::vector< std::set< uint32_t > > active_refs; active_refs.reserve(active_closest.size());
		for (auto const& ac : active_closest) {
			active_refs.emplace_back(ac.begin(), ac.end());
		}
		std::vector< std::set< uint32_t > > next_refs; next_refs.reserve(next_closest.size());
		for (auto const& nc : next_closest) {
			next_refs.emplace_back(nc.begin(), nc.end());
		}

		uint32_t discarded = 0;

		for (auto& ac : active_closest) {
			uint32_t ai = &ac - &active_closest[0];
			for (auto& c : ac) {
				if (c != -1U && !next_refs[c].count(ai)) {
					c = -1U;
					++discarded;
				}
			}
		}

		for (auto& nc : next_closest) {
			uint32_t ni = &nc - &next_closest[0];
			for (auto& c : nc) {
				if (c != -1U && !active_refs[c].count(ni)) {
					c = -1U;
					++discarded;
				}
			}
		}

		if (discarded) qDebug()<< "Discarded " << discarded << " non-mutual segment matches.";

		return discarded > 0;
	};

	//start by removing any non-mutual links:
	discard_nonmutual();

	//fill discarded areas with adjacent information, and process further to make sure the result can be flattened:
	while (true) {
		for (auto& closest : active_closest) {
			uint32_t ai = &closest - &active_closest[0];
			auto const& lengths = active_lengths[ai];
			assert(lengths.size() == closest.size() + 1);

			bool is_loop = active_chains[ai].empty() || active_chains[ai][0] == active_chains[ai].back();

			std::vector< float > weights; weights.reserve(closest.size());
			for (uint32_t i = 1; i < lengths.size(); ++i) {
				weights.emplace_back(lengths[i] - lengths[i - 1]);
			}
			fillUnassigned(closest, weights, is_loop);
			flatten(closest, weights, is_loop);
		}

		for (auto& closest : next_closest) {
			uint32_t ni = &closest - &next_closest[0];
			auto const& lengths = next_lengths[ni];
			assert(lengths.size() == closest.size() + 1);

			bool is_loop = next_chains[ni].empty() || next_chains[ni][0] == next_chains[ni].back();

			std::vector< float > weights; weights.reserve(closest.size());
			for (uint32_t i = 1; i < lengths.size(); ++i) {
				weights.emplace_back(lengths[i] - lengths[i - 1]);
			}
			fillUnassigned(closest, weights, is_loop);
			flatten(closest, weights, is_loop);
		}

		if (discard_nonmutual()) {
			qDebug() << "NOTE: doing another pass through fill/flatten because additional non-mutual links were discarded.";
		}
		else {
			break;
		}
	}

	//helper structs
	struct BeginEnd {
		BeginEnd(float begin_, float end_) : begin(begin_), end(end_) { }
		float begin;
		float end;
	};

	struct BeginEndStitches : public BeginEnd {
		BeginEndStitches(float begin_, float end_, uint32_t next_) : BeginEnd(begin_, end_), next(next_) { }
		std::vector< uint32_t > stitches;
		uint32_t next; //convenience field for split balancing
	};

	struct BeginEndStitches2 : public BeginEnd {
		BeginEndStitches2(float begin_, float end_, uint32_t active_) : BeginEnd(begin_, end_), active(active_) { }
		uint32_t stitches = 0; //std::vector< uint32_t > stitches;
		uint32_t active; //convenience field for merge balancing
	};

	//sort active and next into matching parametric ranges:
	struct Match {
		std::vector< BeginEndStitches > active; //at most two (after post-processing)
		std::vector< BeginEndStitches2 > next; //at most two (after post-processing)
	};

	std::map< std::pair< uint32_t, uint32_t >, Match > matches;


	//build parametric segments into matches:
	for (auto& closest : active_closest) {
		uint32_t ai = &closest - &active_closest[0];
		auto const& lengths = active_lengths[ai];
		assert(lengths.size() == closest.size() + 1);

		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			assert(end < lengths.size());
			matches[std::make_pair(ai, closest[begin])].active.emplace_back(lengths[begin] / lengths.back(), lengths[end] / lengths.back(), closest[begin]);
			auto m = matches[std::make_pair(ai, closest[begin])].active;
			begin = end;
		}
	}
	qDebug() << "matches size " << matches.size();

	for (auto& closest : next_closest) {
		uint32_t ni = &closest - &next_closest[0];

		assert(ni < next_lengths.size());
		auto const& lengths = next_lengths[ni];
		assert(lengths.size() == closest.size() + 1);

		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			assert(end < lengths.size());
			matches[std::make_pair(closest[begin], ni)].next.emplace_back(lengths[begin] / lengths.back(), lengths[end] / lengths.back(), closest[begin]);
			begin = end;
		}
	}


	{ //do stitch assignments:
		//sort ranges from matches back to actives:
		std::vector< std::vector< BeginEndStitches* > > active_segments(active_chains.size());
		for (auto& anm : matches) {
			if (anm.first.first == -1U) {
				assert(anm.second.active.empty());
				continue;
			}
			for (auto& bse : anm.second.active) {
				assert(bse.next == anm.first.second); //'next' should be set correctly!
				active_segments[anm.first.first].emplace_back(&bse);
			}
		}

		//assign stitches to segments:
		for (uint32_t ai = 0; ai < active_chains.size(); ++ai) {
			auto& segments = active_segments[ai];
			assert(!segments.empty());

			std::sort(segments.begin(), segments.end(), [](BeginEndStitches const* a, BeginEndStitches const* b) {
				return a->begin < b->begin;
				});


			//segments should partition [0.0, 1.0):
			assert(!segments.empty());
			assert(segments[0]->begin == 0.0f);
			assert(segments.back()->end == 1.0f);
			for (uint32_t i = 1; i < segments.size(); ++i) {
				assert(segments[i - 1]->end == segments[i]->begin);
			}

			//copy stitches to segments:
			auto const& stitches = active_stitches[ai];
			auto si = stitches.begin();
			for (auto seg : segments) {
				if (si == stitches.end()) break;
				assert(si->t >= seg->begin); //segments should cover everything.
				while (si != stitches.end() && (si->t < seg->end)) {
					seg->stitches.emplace_back(si - stitches.begin());
					++si;
				}
			}
			assert(si == stitches.end()); //segments should cover everything.
		}
	}


	uint32_t active_merges = 0;
	for (auto& anm : matches) {
		if (anm.second.active.size() <= 1) continue;
		assert(anm.first.first != -1U);
		uint32_t ai = anm.first.first;
		bool is_loop = (active_chains[ai].empty() || active_chains[ai][0] == active_chains[ai].back());
		if (!is_loop) continue;
		std::sort(anm.second.active.begin(), anm.second.active.end(), [](BeginEndStitches const& a, BeginEndStitches const& b) {
			return a.begin < b.begin;
			});
		if (anm.second.active[0].begin == 0.0f && anm.second.active.back().end == 1.0f) {
			++active_merges;
			anm.second.active.back().end = anm.second.active[0].end;
			anm.second.active.back().stitches.insert(anm.second.active.back().stitches.end(), anm.second.active[0].stitches.begin(), anm.second.active[0].stitches.end());
			anm.second.active.erase(anm.second.active.begin());
		}
		assert(anm.second.active.size() <= 2); //<-- should be guaranteed by flatten
	}
	if (active_merges) qDebug() << "Merged " << active_merges << " active segments.";

	uint32_t next_merges = 0;
	for (auto& anm : matches) {
		if (anm.second.next.size() <= 1) continue;
		assert(anm.first.second != -1U);
		uint32_t ni = anm.first.second;
		bool is_loop = (next_chains[ni].empty() || next_chains[ni][0] == next_chains[ni].back());
		if (!is_loop) continue;
		std::sort(anm.second.next.begin(), anm.second.next.end(), [](BeginEnd const& a, BeginEnd const& b) {
			return a.begin < b.begin;
			});
		if (anm.second.next[0].begin == 0.0f && anm.second.next.back().end == 1.0f) {
			++next_merges;
			anm.second.next.back().end = anm.second.next[0].end;
			anm.second.next.erase(anm.second.next.begin());
		}
		assert(anm.second.next.size() <= 2); //<-- should be guaranteed by flatten
	}
	if (next_merges) qDebug() << "Merged " << next_merges << " next segments.";


	{ //balance stitch assignments for splits:
		//sort ranges from matches back to actives: (doing again because of merging)
		std::vector< std::vector< BeginEndStitches* > > active_segments(active_chains.size());
		for (auto& anm : matches) {
			if (anm.first.first == -1U) {
				assert(anm.second.active.empty());
				continue;
			}
			for (auto& bse : anm.second.active) {
				assert(bse.next == anm.first.second); //'next' should be set correctly!
				active_segments[anm.first.first].emplace_back(&bse);
			}
		}

		for (uint32_t ai = 0; ai < active_chains.size(); ++ai) {
			auto& segments = active_segments[ai];
			assert(!segments.empty());

			std::sort(segments.begin(), segments.end(), [](BeginEndStitches const* a, BeginEndStitches const* b) {
				return a->begin < b->begin;
				});

			//Note: because of merging, last segment may wrap around weirdly;
			// -- this is okay!
			bool is_loop = (active_chains[ai].empty() || active_chains[ai][0] == active_chains[ai].back());
			assert(
				(segments[0]->begin == 0.0f && segments.back()->end == 1.0f)
				|| (is_loop && (segments[0]->begin == segments.back()->end))
			);
			//segments should still partition [0,1):
			for (uint32_t i = 1; i < segments.size(); ++i) {
				assert(segments[i - 1]->end == segments[i]->begin);
			}

			//find counts (on the way to finding a singleton segment):
			std::unordered_map< uint32_t, uint32_t > next_counts;
			for (auto seg : segments) {
				next_counts.insert(std::make_pair(seg->next, 0)).first->second += 1;
			}

			uint32_t singles = 0;
			uint32_t doubles = 0;
			uint32_t multis = 0;
			for (auto const& nc : next_counts) {
				if (nc.second == 1) ++singles;
				else if (nc.second == 2) ++doubles;
				else ++multis;
			}
			assert(multis == 0);

			if (singles == 1 && doubles == 0 && multis == 0) {
				//just one segment, nothing to re-assign.
				continue;
			}
			else if (singles == 2 && doubles == 0 && multis == 0) {
				//just two segments, *also* always balanced
				continue;
			}
			else if (!(singles == 2 && multis == 0)) {
				throw std::runtime_error("Unhandled split situation with " + std::to_string(singles) + " singles, " + std::to_string(doubles) + " doubles, and " + std::to_string(multis) + " multis.");
			}
			assert(singles == 2 && multis == 0);

			//rotate until a single is in the first position:
			for (uint32_t s = 0; s < segments.size(); ++s) {
				if (next_counts[segments[s]->next] == 1) {
					std::rotate(segments.begin(), segments.begin() + s, segments.end());
					break;
				}
			}
			assert(next_counts[segments[0]->next] == 1);


			//DEBUG: show counts (and adjustments)
			std::string old_back;
			std::string old_front;
			auto pad = [](std::string s) -> std::string {
				while (s.size() < 3) s = ' ' + s;
				return s;
			};
			old_back += pad(std::to_string(segments[0]->stitches.size()));
			old_front += pad("");
			for (uint32_t i = 1; i < segments.size() / 2; ++i) {
				uint32_t io = segments.size() - i;
				assert(segments[i]->next == segments[io]->next);
				old_back += pad(std::to_string(segments[i]->stitches.size()));
				old_front += pad(std::to_string(segments[io]->stitches.size()));
			}
			old_back += pad("");
			old_front += pad(std::to_string(segments[segments.size() / 2]->stitches.size()));
			//end DEBUG



			//expecting things to look like this now (doubles, in order, then another single, then the doubles, reversed):
			// a b c d c b
			assert(segments.size() % 2 == 0);
			assert(segments.size() >= 4);
			assert(next_counts[segments[segments.size() / 2]->next] == 1);
			for (uint32_t i = 1; i < segments.size() / 2; ++i) {
				assert(segments[i]->next == segments[segments.size() - i]->next);
			}

			//now push extra stitches from the center outward:
			uint32_t m = (segments.size() / 2) / 2; //so a b c b -> m = 1 ; a b c d e d c b -> m = 2;
			assert(0 < m && m < segments.size() / 2);
			{ //push from the middle segment outward (alternating sides):
				uint32_t mo = segments.size() - m;
				assert(segments.size() / 2 < mo && mo < segments.size());
				uint32_t mo_p = mo - 1; assert(mo_p < segments.size());
				uint32_t mo_n = (mo + 1 < segments.size() ? mo + 1 : 0); assert(mo_n < segments.size());
				//layout:
				//   m-1  m  m+1
				//   mo_n mo mo_p

				assert(segments[mo]->next == segments[m]->next);
				while (segments[m]->stitches.size() > segments[mo]->stitches.size()) {
					segments[m + 1]->stitches.insert(segments[m + 1]->stitches.begin(), segments[m]->stitches.back());
					segments[m]->stitches.pop_back();
					if (segments[m]->stitches.size() > segments[mo]->stitches.size()) {
						segments[m - 1]->stitches.push_back(segments[m]->stitches[0]);
						segments[m]->stitches.erase(segments[m]->stitches.begin());
					}
				}
				while (segments[mo]->stitches.size() > segments[m]->stitches.size()) {
					segments[mo_p]->stitches.push_back(segments[mo]->stitches[0]);
					segments[mo]->stitches.erase(segments[mo]->stitches.begin());
					if (segments[mo]->stitches.size() > segments[m]->stitches.size()) {
						segments[mo_n]->stitches.insert(segments[mo_n]->stitches.begin(), segments[mo]->stitches.back());
						segments[mo]->stitches.pop_back();
					}
				}
			}
			//push from left-side segments leftward:
			for (uint32_t l = m - 1; l > 0; --l) {
				uint32_t lo = segments.size() - l;
				uint32_t lo_n = (lo + 1 < segments.size() ? lo + 1 : 0); assert(lo_n < segments.size());
				assert(segments[lo]->next == segments[l]->next);
				while (segments[l]->stitches.size() > segments[lo]->stitches.size()) {
					segments[l - 1]->stitches.push_back(segments[l]->stitches[0]);
					segments[l]->stitches.erase(segments[l]->stitches.begin());
				}
				while (segments[lo]->stitches.size() > segments[l]->stitches.size()) {
					segments[lo_n]->stitches.insert(segments[lo_n]->stitches.begin(), segments[lo]->stitches.back());
					segments[lo]->stitches.pop_back();
				}
			}
			//push from right-side segments leftward:
			for (uint32_t r = m + 1; r < segments.size() / 2; ++r) {
				uint32_t ro = segments.size() - r;
				uint32_t ro_p = ro - 1; assert(ro_p < segments.size());
				assert(segments[ro]->next == segments[r]->next);
				while (segments[r]->stitches.size() > segments[ro]->stitches.size()) {
					segments[r + 1]->stitches.insert(segments[r + 1]->stitches.begin(), segments[r]->stitches.back());
					segments[r]->stitches.pop_back();
				}
				while (segments[ro]->stitches.size() > segments[r]->stitches.size()) {
					segments[ro_p]->stitches.push_back(segments[ro]->stitches[0]);
					segments[ro]->stitches.erase(segments[ro]->stitches.begin());
				}
			}

			//check for balance:
			for (uint32_t i = 1; i < segments.size() / 2; ++i) {
				uint32_t io = segments.size() - i;
				assert(segments[i]->stitches.size() == segments[io]->stitches.size());
			}

			//DEBUG: show counts (and adjustments?)
			std::string new_back;
			std::string new_front;
			new_back += pad(std::to_string(segments[0]->stitches.size()));
			new_front += pad("");
			for (uint32_t i = 1; i < segments.size() / 2; ++i) {
				uint32_t io = segments.size() - i;
				assert(segments[i]->next == segments[io]->next);
				new_back += pad(std::to_string(segments[i]->stitches.size()));
				new_front += pad(std::to_string(segments[io]->stitches.size()));
			}
			new_back += pad("");
			new_front += pad(std::to_string(segments[segments.size() / 2]->stitches.size()));
			//end DEBUG

			if (old_back != new_back || old_front != new_front) {
				qDebug() << "Balanced a merge: old: " << old_back << old_front << " new: " << new_back  << new_front;
			}

		}
	}


	//Look for merges/splits:
	std::vector< uint32_t > active_matches(active_chains.size(), 0);
	std::vector< uint32_t > next_matches(next_chains.size(), 0);
	uint32_t empty_matches = 0;
	for (auto const& anm : matches) {
		if (anm.first.first != -1U) active_matches[anm.first.first] += 1;
		if (anm.first.second != -1U) next_matches[anm.first.second] += 1;
		if (anm.first.first == -1U || anm.first.second == -1U) ++empty_matches;
	}
	if (empty_matches) qDebug() << "NOTE: have " << empty_matches << " segments that match with nothing.";

	{ //If there are any merges or splits, all participating next cycles are marked 'accept':
		std::set< uint32_t > to_mark;
		for (auto const& anm : matches) {
			if (anm.first.first == -1U || anm.first.second == -1U) continue;
			bool is_split_or_merge = (active_matches[anm.first.first] > 1 || next_matches[anm.first.second] > 1);
			if (is_split_or_merge) {
				to_mark.insert(anm.first.second);
			}
		}
		uint32_t were_marked = 0;
		for (auto ni : to_mark) {
			for (auto const& td : next_discard_after[ni]) {
				if (td.second) ++were_marked;
			}
			next_discard_after[ni].assign(1, std::make_pair(0.0f, false));
		}
		if (!to_mark.empty() && were_marked != 0) {
			qDebug() << "Marked " << were_marked << " segments on " << to_mark.size() << " next cycles as 'accept' based on participating in a merge/split.";
		}
	}
	//allocate next stitches:

	next_stitches.assign(next_chains.size(), std::vector< Stitch >());

	//allocate stitch counts based on segment lengths + source stitches:
	//(and limit based on active stitch counts)
	//then make next stitches
	for (auto& anm : matches) {
		Match& match = anm.second;

		if (match.active.empty()) {
			qDebug() << "Ignoring match with empty active chain.";
			continue;
		}
		else if (match.next.empty()) {
			if (active_matches[anm.first.first] == 1) {
				qDebug() << "WARNING: active chain matches nothing at all; will not be linked and will thus be discarded.";
			}
			else {
				qDebug() << "Ignoring match with empty next chain.";
			}

			continue;
		}	
		assert(!match.active.empty());
		assert(!match.next.empty());
		bool is_split_or_merge = (active_matches[anm.first.first] > 1 || next_matches[anm.first.second] > 1);
		//if is_split_or_merge will only link 1-1
		
		//compute min/max totals from assigned stitches:
		uint32_t active_ones = 0;
		uint32_t active_anys = 0;
		{
			std::vector< Stitch > const& stitches = active_stitches[anm.first.first];
			for (auto const& be : match.active) {
				for (auto si : be.stitches) {
					assert(si < stitches.size());
					if (stitches[si].flag == Stitch::FlagLinkOne) ++active_ones;
					else if (stitches[si].flag == Stitch::FlagLinkAny) ++active_anys;
					else assert(stitches[si].flag == Stitch::FlagLinkOne || stitches[si].flag == Stitch::FlagLinkAny);
				}
			}
		}
		
		//compute number of ones based on discard segments:
		//Ideally we'd like this to be true (saves on discarded stitches):
		//  discarded segments contain at least one stitch marked 'LinkOne'
		//  while kept segments contain at least two stitches marked 'LinkOne'
		//What we're doing below actually puts two stitches in discarded segments as well.
		uint32_t next_ones = 0;
		
		bool next_is_loop = (next_chains[anm.first.second][0] == next_chains[anm.first.second].back());
		auto const& discard_after = next_discard_after[anm.first.second];
		{ //compute next_ones:
			if (next_is_loop) {
				assert(discard_after[0].second == discard_after.back().second);
			}
			//every discard/non-discard edge needs at least one stitch next to it on each side:
			for (auto const& be : match.next) {
				for (auto tdi = discard_after.begin(); tdi != discard_after.end(); ++tdi) {
					if (next_is_loop && tdi == discard_after.begin()) continue;
					if (be.begin <= tdi->first && tdi->first < be.end) ++next_ones;
					if (be.begin < tdi->first && tdi->first <= be.end) ++next_ones;
				}
			}

			if (next_ones > 2 * active_anys + active_ones) {
				qDebug() << "ERROR: more discard/non-discard ends are required (" << next_ones << ") than are permitted by the current active flags (" << active_anys << "*2 + " << active_ones << "); code to fix this (by removing shortest same-discard segment) not yet implemented.";
				assert(next_ones <= 2 * active_anys + active_ones);
			}
		}

		//compute desired stitch count based on segment lengths:
		assert(anm.first.second < next_lengths.size());
		std::vector< float > const& lengths = next_lengths[anm.first.second];
		float total_length = 0.0f;
		for (auto const& be : match.next) {
			if (be.begin <= be.end) {
				total_length += be.end - be.begin;
			}
			else {
				assert(be.begin < be.end + 1.0f);
				total_length += (be.end + 1.0f) - be.begin;
			}
		}
		//std::cout << std::endl; //DEBUG
		total_length *= lengths.back();

		float stitch_width = stitchWidth / modelUnitLength;
		uint32_t stitches = std::max(1, int32_t(std::round(total_length / stitch_width)));

		{ //adjust for possible links:
			//least is to link 1-1 for every next_ones and then link everything else 2-1:
			uint32_t lower = next_ones //next ones to link to one stitch
				+ (std::max(0, int32_t(active_ones + active_anys) - int32_t(next_ones)) + 1) / 2 //other stitches link 2-1
				;
			uint32_t upper = active_ones + 2 * active_anys; //most is to increase from all anys
			assert(lower <= upper);

			if (is_split_or_merge) {
				assert(lower <= active_ones + active_anys && active_ones + active_anys <= upper);
				qDebug() << "NOTE: setting stitches from " << stitches << " to ";
				stitches = active_ones + active_anys;
				qDebug() << stitches << " to make split/merge 1-1.";
				//stitches = active_ones + active_anys;
			}

			if (stitches < lower || stitches > upper) {
				qDebug() << "NOTE: stitches (" << stitches << ") will be clamped to possible range [" << lower << ", " << upper << "], which might cause some shape distortion.";
				stitches = std::max(lower, std::min(upper, stitches));
			}
			qDebug() << "Will make " << stitches << " stitches, given active with " << active_ones << " ones, " << active_anys << " anys; next with " << next_ones << " ones."; //DEBUG
		}

		std::vector< Stitch > new_stitches;
		if (stitches > 0) {
			//spread stitches among "allocation ranges" (same discard status)
			struct Alloc {
				Alloc(float begin_, float end_, bool first_one_, bool last_one_, uint32_t bi_) : begin(begin_), end(end_), first_one(first_one_), last_one(last_one_), bi(bi_) { assert(begin < end); }
				float begin, end;
				bool first_one, last_one;
				uint32_t bi; //<-- BeginEnd range this came from
				uint32_t stitches = 0;
				float length = std::numeric_limits< float >::quiet_NaN();
			};
			std::vector< Alloc > alloc;
			//NOTE: care is taken so that when stitches are added to the match.next.stitches[] arrays during creation, they will be in CCW order:
			for (auto const& be : match.next) {
				uint32_t bi = &be - &match.next[0];
				auto split_back = [&]() {
					//split allocation range on discards:
					for (auto tdi = discard_after.begin(); tdi != discard_after.end(); ++tdi) {
						if (next_is_loop && tdi == discard_after.begin()) continue;
						if (tdi->first < alloc.back().begin) {
						}
						else if (tdi->first == alloc.back().begin) {
							alloc.back().first_one = true;
						}
						else if (tdi->first < alloc.back().end) {
							float end = alloc.back().end;
							alloc.back().last_one = true;
							alloc.back().end = tdi->first;
							alloc.emplace_back(tdi->first, end, true, false, alloc.back().bi);
						}
						else if (tdi->first == alloc.back().end) {
							alloc.back().last_one = true;
						}
						else {
							assert(tdi->first > alloc.back().end);
						}
					}
				};
				if (be.begin <= be.end) {
					alloc.emplace_back(be.begin, be.end, false, false, bi);
					split_back();
				}
				else {
					alloc.emplace_back(be.begin, 1.0f, false, false, bi);
					split_back();
					alloc.emplace_back(0.0f, be.end, false, false, bi);
					split_back();
				}
			}
			uint32_t total_ones = 0;
			for (auto& a : alloc) {
				a.stitches = (a.first_one ? 1 : 0) + (a.last_one ? 1 : 0);
				total_ones += a.stitches;
				a.length = (a.end - a.begin) * lengths.back();
				assert(a.length >= 0.0f);
			}

			assert(total_ones == next_ones); //better have the same number of ones as we accounted for previously

			//add remaining stitches to alloc ranges based on which range has the least-dense stitches:
			for (uint32_t s = total_ones; s < stitches; ++s) {
				uint32_t best = -1U;
				float best_density = 0.0f;
				for (auto const& a : alloc) {
					float d = a.length / float(a.stitches + 1);
					if (d > best_density) {
						best = &a - &alloc[0];
						best_density = d;
					}
				}
				assert(best < alloc.size());
				alloc[best].stitches += 1;
			}

			//actually create stitches from allocation ranges:
			for (auto const& a : alloc) {
				if (a.stitches == 0) continue;
				for (uint32_t s = 0; s < a.stitches; ++s) {
					float t = (s + 0.5f) / float(a.stitches) * (a.end - a.begin) + a.begin;
					Stitch::Flag flag = Stitch::FlagLinkAny;
					if (s == 0 && a.first_one) flag = Stitch::FlagLinkOne;
					if (s + 1 == a.stitches && a.last_one) flag = Stitch::FlagLinkOne;
					new_stitches.emplace_back(t, flag);
					//make sure it's in the range it's being assigned to:
					if (match.next[a.bi].begin < match.next[a.bi].end) {
						assert(match.next[a.bi].begin <= t && t < match.next[a.bi].end);
					}
					else {
						assert(t < match.next[a.bi].end || match.next[a.bi].begin <= t);
					}
					match.next[a.bi].stitches += 1; 
				}
			}
		}
		assert(new_stitches.size() == stitches);

		next_stitches[anm.first.second].insert(next_stitches[anm.first.second].end(), new_stitches.begin(), new_stitches.end());
	} //end stitch allocation


	{ //balance new stitch allocations for merges:

		std::vector< std::vector< BeginEndStitches2* > > next_segments(next_chains.size());

		for (auto& anm : matches) {
			if (anm.first.second == -1U) {
				assert(anm.second.next.empty());
				continue;
			}
			for (auto& bse : anm.second.next) {
				assert(bse.active == anm.first.first);
				next_segments[anm.first.second].emplace_back(&bse);
			}
		}
		for (uint32_t ni = 0; ni < next_chains.size(); ++ni) {
			auto& segments = next_segments[ni];
			assert(!segments.empty());

			std::sort(segments.begin(), segments.end(), [](BeginEndStitches2 const* a, BeginEndStitches2 const* b) {
				return a->begin < b->begin;
				});

			//Note: because of merging, last segment may wrap around weirdly;
			// -- this is okay!
			bool is_loop = (next_chains[ni].empty() || next_chains[ni][0] == next_chains[ni].back());
			assert(
				(segments[0]->begin == 0.0f && segments.back()->end == 1.0f)
				|| (is_loop && (segments[0]->begin == segments.back()->end))
			);
			//segments should still partition [0,1):
			for (uint32_t i = 1; i < segments.size(); ++i) {
				assert(segments[i - 1]->end == segments[i]->begin);
			}

			//find counts (on the way to finding a singleton segment):
			std::unordered_map< uint32_t, uint32_t > active_counts;
			for (auto seg : segments) {
				active_counts.insert(std::make_pair(seg->active, 0)).first->second += 1;
			}

			uint32_t singles = 0;
			uint32_t doubles = 0;
			uint32_t multis = 0;
			for (auto const& nc : active_counts) {
				if (nc.second == 1) ++singles;
				else if (nc.second == 2) ++doubles;
				else ++multis;
			}
			assert(multis == 0);

			if (singles == 1 && doubles == 0 && multis == 0) {
				//just one segment, nothing to re-assign.
				continue;
			}
			else if (singles == 2 && doubles == 0 && multis == 0) {
				//just two segments, *also* always balanced
				continue;
			}
			else if (!(singles == 2 && multis == 0)) {
				throw std::runtime_error("Unhandled merge situation with " + std::to_string(singles) + " singles, " + std::to_string(doubles) + " doubles, and " + std::to_string(multis) + " multis.");
			}
			assert(singles == 2 && multis == 0);

			//rotate until a single is in the first position:
			for (uint32_t s = 0; s < segments.size(); ++s) {
				if (active_counts[segments[s]->active] == 1) {
					std::rotate(segments.begin(), segments.begin() + s, segments.end());
					break;
				}
			}
			assert(active_counts[segments[0]->active] == 1);

			//expecting things to look like this now (doubles, in order, then another single, then the doubles, reversed):
			// a b c d c b
			assert(segments.size() % 2 == 0);
			assert(segments.size() >= 4);
			assert(active_counts[segments[segments.size() / 2]->active] == 1);
			for (uint32_t i = 1; i < segments.size() / 2; ++i) {
				assert(segments[i]->active == segments[segments.size() - i]->active);
			}

			//DEBUG: show counts (and adjustments?)
			std::string old_back;
			std::string old_front;
			auto pad = [](std::string s) -> std::string {
				while (s.size() < 3) s = ' ' + s;
				return s;
			};
			old_back += pad(std::to_string(segments[0]->stitches));
			old_front += pad("");
			for (uint32_t i = 1; i < segments.size() / 2; ++i) {
				uint32_t io = segments.size() - i;
				assert(segments[i]->active == segments[io]->active);
				old_back += pad(std::to_string(segments[i]->stitches));
				old_front += pad(std::to_string(segments[io]->stitches));
			}
			old_back += pad("");
			old_front += pad(std::to_string(segments[segments.size() / 2]->stitches));
			//end DEBUG

			//now walk from left-to-right along internal segments, redoing stitch counts as needed:
			uint32_t sum = 0;
			uint32_t sumo = 0;
			for (uint32_t i = 1; i < segments.size() / 2; ++i) {
				uint32_t io = segments.size() - i;
				assert(segments[i]->active == segments[io]->active);
				uint32_t total = segments[i]->stitches + segments[io]->stitches;
				segments[i]->stitches = total / 2;
				segments[io]->stitches = (total + 1) / 2;
				if (sumo > sum) {
					std::swap(segments[i]->stitches, segments[io]->stitches);
				}
				sum += segments[i]->stitches;
				sumo += segments[io]->stitches;
				assert(std::abs(int32_t(sumo) - int32_t(sum)) <= 1);
			}
			//REMEMBER: as per Appendix A, this algorithm isn't perfect (can fail in very large merge/split chains).

			//Now delete the existing stitches for this chain and re-allocate!
			uint32_t old_count = next_stitches[ni].size();

			next_stitches[ni].clear();
			for (auto seg : segments) {
				for (uint32_t s = 0; s < seg->stitches; ++s) {
					float end = (seg->begin <= seg->end ? seg->end : seg->end + 1.0f);
					assert(seg->begin <= end);
					float t = (s + 0.5f) / float(seg->stitches) * (end - seg->begin) + seg->begin;
					if (t >= 1.0f) t -= 1.0f;
					assert(t < 1.0f);

					next_stitches[ni].emplace_back(t, Stitch::FlagLinkAny);
					//make sure it's in the range it's being assigned to:
					if (seg->begin < seg->end) {
						assert(seg->begin <= t && t < seg->end);
					}
					else {
						assert(t < seg->end || seg->begin <= t);
					}
				}
			}
			assert(old_count == next_stitches[ni].size());

			//DEBUG: show counts (and adjustments?)
			std::string new_back;
			std::string new_front;
			new_back += pad(std::to_string(segments[0]->stitches));
			new_front += pad("");
			for (uint32_t i = 1; i < segments.size() / 2; ++i) {
				uint32_t io = segments.size() - i;
				assert(segments[i]->active == segments[io]->active);
				new_back += pad(std::to_string(segments[i]->stitches));
				new_front += pad(std::to_string(segments[io]->stitches));
			}
			new_back += pad("");
			new_front += pad(std::to_string(segments[segments.size() / 2]->stitches));
			//end DEBUG

			if (old_back != new_back || old_front != new_front) {
				qDebug() << "Balanced a split: " << "old: " << old_back << old_front  << " new: " << new_back << new_front;
			}


		}
	}


	for (auto& stitches : next_stitches) {
		std::stable_sort(stitches.begin(), stitches.end(), [](Stitch const& a, Stitch const& b) {
			return a.t < b.t;
			});
	}
	auto make_stitch_info = [&slice](
		std::vector< uint32_t > const& chain,
		std::vector< float > const& lengths,
		std::vector< Stitch > const& stitches,
		std::vector< glm::vec3 >* stitch_locations_,
		std::vector< bool >* stitch_linkones_) {

			assert(chain.size() == lengths.size());

			assert(stitch_locations_);
			auto& stitch_locations = *stitch_locations_;
			stitch_locations.clear();

			assert(stitch_linkones_);
			auto& stitch_linkones = *stitch_linkones_;
			stitch_linkones.clear();

			auto li = lengths.begin();
			for (auto si = stitches.begin(); si != stitches.end(); ++si) {
				float l = si->t * lengths.back();
				assert(si->t >= 0.0f && si->t <= 1.0f);
				while (li != lengths.end() && *li <= l) ++li;
				assert(li != lengths.end());
				assert(li != lengths.begin());
				uint32_t i = li - lengths.begin();

				stitch_locations.emplace_back(
					glm::mix(
						slice.vertices[chain[i - 1]],
						slice.vertices[chain[i]],
						(l - *(li - 1)) / (*li - *(li - 1))
					)
				);
				stitch_linkones.emplace_back(si->flag == Stitch::FlagLinkOne);
			}
	};
	std::vector< std::vector< glm::vec3 > > all_active_stitch_locations(active_chains.size());
	std::vector< std::vector< bool > > all_active_stitch_linkones(active_chains.size());

	std::vector< std::vector< glm::vec3 > > all_next_stitch_locations(next_chains.size());
	std::vector< std::vector< bool > > all_next_stitch_linkones(next_chains.size());

	for (uint32_t ai = 0; ai < active_chains.size(); ++ai) {
		assert(ai < active_stitches.size());
		make_stitch_info(
			active_chains[ai],
			active_lengths[ai],
			active_stitches[ai],
			&all_active_stitch_locations[ai],
			&all_active_stitch_linkones[ai]);
	}
	for (uint32_t ni = 0; ni < next_chains.size(); ++ni) {
		assert(ni < next_stitches.size());
		make_stitch_info(
			next_chains[ni],
			next_lengths[ni],
			next_stitches[ni],
			&all_next_stitch_locations[ni],
			&all_next_stitch_linkones[ni]);
	}

	//PARANOIA:
	std::vector< std::unordered_set< uint32_t > > all_next_claimed(next_chains.size());
	std::vector< std::unordered_set< uint32_t > > all_active_claimed(active_chains.size());

	std::vector< std::unordered_set< uint32_t > > empty_dummy(2);
	qDebug() << "check matches " << matches.size();

	for (auto const& anm : matches) {
		qDebug() << "start of loop";

		Match const& match = anm.second;

		std::vector< uint32_t > next_stitch_indices;
		std::vector< glm::vec3 > next_stitch_locations;
		std::vector< bool > next_stitch_linkones;
		qDebug() << "vars inited" << anm.first.second;

		/*
		if (match.active.empty()) {
			qDebug() << "Ignoring match with empty active chain.";
			continue;
		}
		else if (match.next.empty()) {
			qDebug() << "Ignoring match with empty next chain.";
			continue;
		}*/
		//hacky hacky hacky hack ahcahfavughv
		qDebug() << next_chains.empty() << active_chains.empty();
		std::unordered_set< uint32_t >& next_claimed = anm.first.second == -1U ? empty_dummy[0] : all_next_claimed[anm.first.second];
		

		auto do_range = [&](float begin, float end) {
			//std::cout << "do_range [" << begin << ", " << end << "): "; //DEBUG
			assert(begin <= end);
			auto const& ns = next_stitches[anm.first.second];
			uint32_t count = 0;
			for (auto const& s : ns) {
				if (begin <= s.t && s.t < end) {
					uint32_t si = &s - &ns[0];
					//std::cout << " " << si; std::cout.flush(); //DEBUG
					next_stitch_indices.emplace_back(si);
					next_stitch_locations.emplace_back(all_next_stitch_locations[anm.first.second][si]);
					next_stitch_linkones.emplace_back(all_next_stitch_linkones[anm.first.second][si]);
					auto ret = next_claimed.insert(si); //PARANOIA
					assert(ret.second);
					++count;
				}
			}
			//std::cout << std::endl; //DEBUG
			return count;
		};
		for (auto const& be : match.next) {
			uint32_t count;
			if (be.begin <= be.end) {
				count = do_range(be.begin, be.end);
			}
			else {
				assert(&be == &match.next.back()); //only last one should be split/merge
				count = do_range(be.begin, 1.0f);
				count += do_range(0.0f, be.end);
			}
			assert(count == be.stitches); //should have found the right number of stitches
		}
		qDebug() << "wowee";
		{ //PARANOIA: because of the care we took in the look-up above, should have (at most one) decrease in t-coord in the next_stitches:
			uint32_t decreases = 0;
			for (uint32_t i = 0; i < next_stitch_indices.size(); ++i) {
				float t0 = next_stitches[anm.first.second][next_stitch_indices[i == 0 ? next_stitch_indices.size() - 1 : i - 1]].t;
				float t1 = next_stitches[anm.first.second][next_stitch_indices[i]].t;
				if (t0 < t1) {
					//expected
				}
				else {
					assert(t0 < t1 + 1.0f); //wrapped
					decreases += 1;
				}
			}
			assert(decreases <= 1);
		}

		qDebug() << "check6";
		std::vector< uint32_t > active_stitch_indices;
		std::vector< glm::vec3 > active_stitch_locations;
		std::vector< bool > active_stitch_linkones;
		std::unordered_set< uint32_t >& active_claimed = all_active_claimed[anm.first.first];
		for (auto const& be : match.active) {
			for (auto si : be.stitches) {
				active_stitch_indices.emplace_back(si);
				active_stitch_locations.emplace_back(all_active_stitch_locations[anm.first.first][si]);
				active_stitch_linkones.emplace_back(all_active_stitch_linkones[anm.first.first][si]);
				auto ret = active_claimed.insert(si); //PARANOIA
				assert(ret.second);
			}
		}

		{ //PARANOIA: because of the care we take in book-keeping, should have (at most one) decrease in t-coord in the active_stitches:
			uint32_t decreases = 0;
			for (uint32_t i = 0; i < active_stitch_indices.size(); ++i) {
				float t0 = active_stitches[anm.first.first][active_stitch_indices[i == 0 ? active_stitch_indices.size() - 1 : i - 1]].t;
				float t1 = active_stitches[anm.first.first][active_stitch_indices[i]].t;
				if (t0 < t1) {
					//expected
				}
				else {
					assert(t0 < t1 + 1.0f); //wrapped
					decreases += 1;
				}
			}
			assert(decreases <= 1);
		}


		if (match.active.empty()) {
			qDebug() << "Ignoring match with empty active chain.";
			continue;
		}
		else if (match.next.empty()) {
			qDebug() << "Ignoring match with empty next chain.";
			continue;
		}
		assert(!match.active.empty());
		assert(!match.next.empty());


		//figure out how to link to next stitches:
		{ //DEBUG:
			uint32_t new_ones = 0;
			for (auto o : next_stitch_linkones) {
				if (o) ++new_ones;
			}
			uint32_t active_ones = 0;
			for (auto o : active_stitch_linkones) {
				if (o) ++active_ones;
			}
			qDebug() << " About to connect " << active_stitch_locations.size() << " active stitches (" << active_ones << " linkones) to " << next_stitch_locations.size() << " new stitches (" << new_ones << " linkones).";
		}

		//actually build links:
		{ //least-clever linking solution: FlagLinkOne's link 1-1, others link to keep arrays mostly in sync

			std::vector< std::pair< uint32_t, uint32_t > > best_links;
			float best_cost = std::numeric_limits< float >::infinity();
#define USE_OPTIMAL
#ifdef USE_OPTIMAL
			optimalLink(2.0f * stitchHeight / modelUnitLength,
				true, //allow roll
				active_stitch_locations, active_stitch_linkones,
				next_stitch_locations, next_stitch_linkones,
				&best_links);
			(void)best_cost; //unused!
#endif
			for (auto const& p : best_links) {
				Link link;
				link.from_chain = anm.first.first;
				assert(p.first < active_stitch_indices.size());
				link.from_stitch = active_stitch_indices[p.first];

				link.to_chain = anm.first.second;
				assert(p.second < next_stitch_indices.size());
				link.to_stitch = next_stitch_indices[p.second];
				links.emplace_back(link);
				}
			}
	}


		//PARANOIA: every stitch should have been claimed
	
		for (uint32_t ai = 0; ai < active_chains.size(); ++ai) {
			qDebug() << all_active_claimed[ai].size() << " vs " << active_stitches[ai].size();
			assert(all_active_claimed[ai].size() == active_stitches[ai].size());
		}
		for (uint32_t ni = 0; ni < next_chains.size(); ++ni) {
			assert(all_next_claimed[ni].size() == next_stitches[ni].size());
		}

		//mark next stitches in discard range as 'discard':
		assert(next_discard_after.size() == next_stitches.size());
		uint32_t marked = 0;
		uint32_t total = 0;
		for (auto& stitches : next_stitches) {
			auto const& discard_after = next_discard_after[&stitches - &next_stitches[0]];

			auto di = discard_after.begin();
			for (auto& s : stitches) {
				assert(di != discard_after.end());
				while (di + 1 != discard_after.end() && (di + 1)->first <= s.t) ++di;
				assert(di->first <= s.t);
				if (di->second) {
					s.flag = Stitch::FlagDiscard;
					++marked;
				}
				++total;
			}
		}
	qDebug() << "Marked " << marked << " of " << total << " newly created stitches as 'discard'.";

}



/*
*	Function: compute the shortest path from a source vertex to a target vertex
*	Return: void, but the path between the vertices is passed through the path_ argument
*	Called: called during the buildNextActiveChains() function
*/
void KnitGrapher::embeddedPathSimple(
	AutoKnitMesh const& model,
	EmbeddedVertex const& source,
	EmbeddedVertex const& target,
	std::vector< EmbeddedVertex >* path_ //out: path; path[0] will be source and path.back() will be target
) 
{

	assert(source != target);

	assert(path_);
	auto& path = *path_;
	path.clear();


	//idea: place distance storage along each edge and on each corner.
	std::vector< EmbeddedVertex > loc_ev;
	std::vector< glm::vec3 > loc_pos;

	//TODO: eventually, incrementally build locs starting from source / target verts

	//add locs for each vertex in the mesh:
	std::vector< uint32_t > vertex_locs;
	vertex_locs.reserve(model.vertices.size());
	for (uint32_t vi = 0; vi < model.vertices.size(); ++vi) {
		vertex_locs.emplace_back(loc_ev.size());
		loc_ev.emplace_back(EmbeddedVertex::on_vertex(vi));
		loc_pos.emplace_back(loc_ev.back().interpolate(model.vertices));
	}

	//make an edge-to-triangle look-up structure:
	std::unordered_multimap< glm::uvec2, uint32_t > edge_triangles;
	std::unordered_map< glm::uvec3, uint32_t > simplex_triangle;
	std::unordered_set< glm::uvec2 > edges;

	//std::vector<glm::uvec3> modelTriangles = getTriangles(model);

	for (auto const& tri : model.triangles) {
		uint32_t ti = &tri - &model.triangles[0];
		auto do_edge = [&](uint32_t a, uint32_t b) {
			if (a > b) std::swap(a, b);
			edge_triangles.insert(std::make_pair(glm::uvec2(a, b), ti));
			edges.insert(glm::uvec2(a, b));
		};
		do_edge(tri.x, tri.y);
		do_edge(tri.y, tri.z);
		do_edge(tri.z, tri.x);

		glm::uvec3 simplex = tri;
		if (simplex.x > simplex.y) std::swap(simplex.x, simplex.y);
		if (simplex.y > simplex.z) std::swap(simplex.y, simplex.z);
		if (simplex.x > simplex.y) std::swap(simplex.x, simplex.y);
		auto ret = simplex_triangle.insert(std::make_pair(simplex, ti));
		assert(ret.second);
	}

	//add (several?) locs along each edge:
	std::unordered_map< glm::uvec2, std::pair< uint32_t, uint32_t > > edge_locs;
	edge_locs.reserve(edges.size());
	float const max_spacing = getMaxPathSampleSpacing();
	for (auto const& e : edges) {
		uint32_t count = std::max(0, int32_t(std::floor(glm::length(model.vertices[e.y] - model.vertices[e.x]) / max_spacing)));
		uint32_t begin = loc_ev.size();
		uint32_t end = begin + count;
		edge_locs.insert(std::make_pair(e, std::make_pair(begin, end)));
		for (uint32_t i = 0; i < count; ++i) {
			loc_ev.emplace_back(EmbeddedVertex::on_edge(e.x, e.y, (i + 0.5f) / float(count)));
			loc_pos.emplace_back(loc_ev.back().interpolate(model.vertices));
		}
		assert(loc_ev.size() == end);
	}

	//build adjacency lists for each triangle:
	std::vector< std::vector< uint32_t > > loc_tris(loc_ev.size());
	std::vector< std::vector< uint32_t > > tri_adj(model.triangles.size());
	for (auto const& tri : model.triangles) {
		uint32_t ti = &tri - &model.triangles[0];
		auto do_edge = [&](uint32_t a, uint32_t b) {
			if (a > b) std::swap(b, a);
			auto f = edge_locs.find(glm::uvec2(a, b));
			assert(f != edge_locs.end());
			for (uint32_t i = f->second.first; i < f->second.second; ++i) {
				loc_tris[i].emplace_back(ti);
				tri_adj[ti].emplace_back(i);
			}
		};

		auto do_vertex = [&](uint32_t a) {
			uint32_t i = vertex_locs[a];
			loc_tris[i].emplace_back(ti);
			tri_adj[ti].emplace_back(i);
		};

		do_edge(tri.x, tri.y);
		do_edge(tri.y, tri.z);
		do_edge(tri.z, tri.x);
		do_vertex(tri.x);
		do_vertex(tri.y);
		do_vertex(tri.z);
	}


	//add source and target to the locs lists:
	auto add_embedded = [&](EmbeddedVertex const& ev) {
		assert(ev.simplex.x != -1U);
		if (ev.simplex.y == -1U) {
			//at a vertex
			return vertex_locs[ev.simplex.x];
		}
		else if (ev.simplex.z == -1U) {
			//on an edge
			uint32_t idx = loc_ev.size();
			loc_ev.emplace_back(ev);
			loc_pos.emplace_back(loc_ev.back().interpolate(model.vertices));

			loc_tris.emplace_back();
			auto r = edge_triangles.equal_range(glm::uvec2(ev.simplex));
			assert(r.first != r.second);
			for (auto ri = r.first; ri != r.second; ++ri) {
				uint32_t ti = ri->second;
				loc_tris.back().emplace_back(ti);
				tri_adj[ti].emplace_back(idx);
			}

			return idx;
		}
		else {
			//on a triangle
			uint32_t idx = loc_ev.size();
			loc_ev.emplace_back(ev);
			loc_pos.emplace_back(loc_ev.back().interpolate(model.vertices));

			auto f = simplex_triangle.find(ev.simplex);
			assert(f != simplex_triangle.end());
			uint32_t ti = f->second;

			loc_tris.emplace_back();
			loc_tris.back().emplace_back(ti);
			tri_adj[ti].emplace_back(idx);

			return idx;
		}
	};

	uint32_t source_idx = add_embedded(source);
	uint32_t target_idx = add_embedded(target);


	//now do actual search:

	std::vector< float > loc_dis(loc_pos.size(), std::numeric_limits< float >::infinity());
	std::vector< uint32_t > loc_from(loc_pos.size(), -1U);

	glm::vec3 target_pos = target.interpolate(model.vertices);

	std::vector< std::pair< float, std::pair< uint32_t, float > > > todo;

	auto queue = [&](uint32_t at, float distance, uint32_t from) {
		assert(distance < loc_dis[at]);
		loc_dis[at] = distance;
		loc_from[at] = from;

		float heuristic = glm::length(target_pos - loc_pos[at]);
		todo.emplace_back(std::make_pair(-(heuristic + distance), std::make_pair(at, distance)));
		std::push_heap(todo.begin(), todo.end());
	};

	queue(source_idx, 0.0f, -1U);
	while (!todo.empty()) {
		std::pop_heap(todo.begin(), todo.end());
		uint32_t at = todo.back().second.first;
		float distance = todo.back().second.second;
		todo.pop_back();

		if (distance > loc_dis[at]) continue;
		if (at == target_idx) break; //bail out early -- don't need distances to everything.

		assert(distance == loc_dis[at]);
		for (auto t : loc_tris[at]) {
			for (auto n : tri_adj[t]) {
				if (n == at) continue;
				float d = distance + glm::length(loc_pos[n] - loc_pos[at]);
				if (d < loc_dis[n]) queue(n, d, at);
			}
		}
	}

	//read back path:
	if (loc_from[target_idx] == -1U) {
		throw std::runtime_error("embedded_path requested between disconnected vertices");
	}

	uint32_t at = target_idx;
	do {
		path.emplace_back(loc_ev[at]);
		at = loc_from[at];
	} while (at != -1U);
	assert(path.size() >= 2);
	std::reverse(path.begin(), path.end());

	assert(path[0] == source);
	assert(path.back() == target);

}

/*
*	Function: compute a path from source to target estimating the distance between the source and target 
*		and then trimming the mesh to include only the triangles that could be traversed by the path.
*	Return: void, but the path between the source and target is passed through the path_ argument
*	Called: called during the builNextActiveChain() function
*/
void KnitGrapher::embeddedPath(
	AutoKnitMesh const& model,
	EmbeddedVertex const& source,
	EmbeddedVertex const& target,
	std::vector< EmbeddedVertex >* path_ //out: path; path[0] will be source and path.back() will be target
) {

	assert(path_);
	auto& path = *path_;
	path.clear();

	//first do a vertex-to-vertex distance computation to bound the computation:

	std::vector< std::vector< uint32_t > > adj(model.vertices.size());

	std::unordered_set< glm::uvec2 > edges;
	for (auto const& tri : model.triangles) {
		auto do_edge = [&](uint32_t a, uint32_t b) {
			if (a > b) std::swap(a, b);
			edges.insert(glm::uvec2(a, b));
		};
		do_edge(tri.x, tri.y);
		do_edge(tri.y, tri.z);
		do_edge(tri.z, tri.x);
	}

	for (auto const& e : edges) {
		adj[e.x].emplace_back(e.y);
		adj[e.y].emplace_back(e.x);
	}

	uint32_t target_idx = target.simplex.x;

	std::vector< float > dis(model.vertices.size(), std::numeric_limits< float >::infinity());

	std::vector< std::pair< float, std::pair< uint32_t, float > > > todo;

	auto queue = [&](uint32_t at, float distance) {
		assert(distance < dis[at]);
		dis[at] = distance;

		float heuristic = glm::length(model.vertices[target_idx] - model.vertices[at]);
		todo.emplace_back(std::make_pair(-(heuristic + distance), std::make_pair(at, distance)));
		std::push_heap(todo.begin(), todo.end());
	};

	queue(source.simplex.x, glm::length(source.interpolate(model.vertices) - model.vertices[source.simplex.x]));

	while (!todo.empty()) {
		std::pop_heap(todo.begin(), todo.end());
		uint32_t at = todo.back().second.first;
		float distance = todo.back().second.second;
		todo.pop_back();

		if (distance > dis[at]) continue;
		assert(distance == dis[at]);

		if (at == target_idx) break; //bail out early -- don't need distances to everything.

		for (auto n : adj[at]) {
			float d = distance + glm::length(model.vertices[n] - model.vertices[at]);
			if (d < dis[n]) queue(n, d);
		}
	}

	//okay, so this is a conservative (long) estimate of path length:
	float dis2 = dis[target_idx] + glm::length(target.interpolate(model.vertices) - model.vertices[target_idx]);
	dis2 = dis2 * dis2;

	//come up with a model containing only triangles that might be used in the path:

	AutoKnitMesh trimmed;
	trimmed.vertices.reserve(model.vertices.size());
	trimmed.triangles.reserve(model.triangles.size());

	std::vector< uint32_t > to_trimmed(model.vertices.size(), -1U);
	std::vector< uint32_t > from_trimmed;
	from_trimmed.reserve(model.vertices.size());
	auto vertex_to_trimmed = [&to_trimmed, &from_trimmed, &trimmed, &model](uint32_t v) {
		if (to_trimmed[v] == -1U) {
			to_trimmed[v] = trimmed.vertices.size();
			from_trimmed.emplace_back(v);
			trimmed.vertices.emplace_back(model.vertices[v]);
		}
		return to_trimmed[v];
	};

	{ //keep triangles that are close enough to source and target that the path could possible pass through them:
		glm::vec3 src = source.interpolate(model.vertices);
		glm::vec3 tgt = target.interpolate(model.vertices);
		for (auto const& tri : model.triangles) {
			glm::vec3 min = glm::min(model.vertices[tri.x], glm::min(model.vertices[tri.y], model.vertices[tri.z]));
			glm::vec3 max = glm::max(model.vertices[tri.x], glm::max(model.vertices[tri.y], model.vertices[tri.z]));

			float len2_src = glm::length2(glm::max(min, glm::min(max, src)) - src);
			float len2_tgt = glm::length2(glm::max(min, glm::min(max, tgt)) - tgt);
			if (len2_src + len2_tgt < dis2) {
				trimmed.triangles.emplace_back(glm::uvec3(
					vertex_to_trimmed(tri.x), vertex_to_trimmed(tri.y), vertex_to_trimmed(tri.z)
				));
			}
		}
	}

	EmbeddedVertex trimmed_source = source;

	trimmed_source.simplex.x = vertex_to_trimmed(trimmed_source.simplex.x);
	if (trimmed_source.simplex.y != -1U) trimmed_source.simplex.y = vertex_to_trimmed(trimmed_source.simplex.y);
	if (trimmed_source.simplex.z != -1U) trimmed_source.simplex.z = vertex_to_trimmed(trimmed_source.simplex.z);
	trimmed_source = EmbeddedVertex::canonicalize(trimmed_source.simplex, trimmed_source.weights);

	EmbeddedVertex trimmed_target = target;
	trimmed_target.simplex.x = vertex_to_trimmed(trimmed_target.simplex.x);
	if (trimmed_target.simplex.y != -1U) trimmed_target.simplex.y = vertex_to_trimmed(trimmed_target.simplex.y);
	if (trimmed_target.simplex.z != -1U) trimmed_target.simplex.z = vertex_to_trimmed(trimmed_target.simplex.z);
	trimmed_target = EmbeddedVertex::canonicalize(trimmed_target.simplex, trimmed_target.weights);

	assert(from_trimmed.size() == trimmed.vertices.size());


	embeddedPathSimple(
		trimmed,
		trimmed_source,
		trimmed_target,
		&path);

	for (auto& v : path) {
		v.simplex.x = from_trimmed[v.simplex.x];
		if (v.simplex.y != -1U) v.simplex.y = from_trimmed[v.simplex.y];
		if (v.simplex.z != -1U) v.simplex.z = from_trimmed[v.simplex.z];
		v = EmbeddedVertex::canonicalize(v.simplex, v.weights);
	}

}



/*
*	Function: build the new active chains that will serve as the jumping off point for the next find active->peel slice->link->build iteration
*	Return: void, but modifies the next active_chains variable to include the info of the next iteration
*		also changes the RowColGraph
*	Called: called when the user clicks the 'step' button in the GUI, is the 4th step after linkChains()
*/
void KnitGrapher::buildNextActiveChains(
	AutoKnitMesh const& slice,
	std::vector< EmbeddedVertex > const& slice_on_model, //in: vertices of slice (on model)
	std::vector< std::vector< uint32_t > > const& active_chains,  //in: current active chains (on slice)
	std::vector< std::vector< Stitch > > const& active_stitches, //in: current active stitches
	std::vector< std::vector< uint32_t > > const& next_chains, //in: next chains (on slice)
	std::vector< std::vector< Stitch > > const& next_stitches, //in: next stitches
	std::vector< bool > const& next_used_boundary, //in: did next chain use boundary?
	std::vector< Link > const& links_in, //in: links between active and next
	std::vector< std::vector< EmbeddedVertex > >* next_active_chains_, //out: next active chains (on model)
	std::vector< std::vector< Stitch > >* next_active_stitches_, //out: next active stitches
	RowColGraph* graph_ //in/out (optional): graph to update
) 
{

	
	for (auto const& chain : active_chains) {
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i - 1] != chain[i]);
		}
	}

	assert(active_stitches.size() == active_chains.size());

	for (auto const& chain : next_chains) {
		for (auto v : chain) {
			assert(v < slice.vertices.size());
		}
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i - 1] != chain[i]);
		}
	}

	assert(next_stitches.size() == next_chains.size());
	assert(next_used_boundary.size() == next_chains.size());


	for (auto const& l : links_in) {
		assert(l.from_chain < active_chains.size());
		assert(l.from_stitch < active_stitches[l.from_chain].size());
		assert(l.to_chain < next_chains.size());
		assert(l.to_stitch < next_stitches[l.to_chain].size());
	}

	assert(next_active_chains_);
	auto& next_active_chains = *next_active_chains_;
	next_active_chains.clear();

	assert(next_active_stitches_);
	auto& next_active_stitches = *next_active_stitches_;
	next_active_stitches.clear();

	//PARANOIA:
	if (graph_) {
		for (auto const& stitches : active_stitches) {
			for (auto const& s : stitches) {
				assert(s.vertex != -1U);
				assert(s.vertex < graph_->vertices.size());
			}
		}
		for (auto const& stitches : next_stitches) {
			for (auto const& s : stitches) {
				assert(s.vertex == -1U);
			}
		}

	}



	//any active chain with no links out is considered inactive and discarded:
	std::vector< bool > discard_active(active_chains.size(), true);
	for (auto const& l : links_in) {
		discard_active[l.from_chain] = false;
	}

	//filter to links that target non-discarded stitches only:
	std::vector< Link > links;
	for (auto const& l : links_in) {
		if (next_stitches[l.to_chain][l.to_stitch].flag == Stitch::FlagDiscard) continue;
		links.emplace_back(l);
	}

	//build a lookup structure for links:
	struct ChainStitch {
		ChainStitch(uint32_t chain_, uint32_t stitch_) : chain(chain_), stitch(stitch_) { }
		uint32_t chain;
		uint32_t stitch;
		bool operator<(ChainStitch const& o) const {
			if (chain != o.chain) return chain < o.chain;
			else return stitch < o.stitch;
		}
		bool operator==(ChainStitch const& o) const {
			return chain == o.chain && stitch == o.stitch;
		}
		bool operator!=(ChainStitch const& o) const {
			return !(*this == o);
		}
	};

	std::map< ChainStitch, std::vector< ChainStitch > > active_next;
	std::map< ChainStitch, std::vector< ChainStitch > > next_active;

	//NOTE: link_chains guarantees that links are in "direction of chain" order, so old sorting code removed:
	for (auto const& l : links) {
		active_next[ChainStitch(l.from_chain, l.from_stitch)]
			.emplace_back(l.to_chain, l.to_stitch);
		next_active[ChainStitch(l.to_chain, l.to_stitch)]
			.emplace_back(l.from_chain, l.from_stitch);
	}



	//record whether the segments adjacent to every next stitch is are marked as "discard" or "keep":
	std::vector< std::vector< std::pair< bool, bool > > > keep_adj(next_chains.size());
	for (uint32_t nc = 0; nc < next_chains.size(); ++nc) {
		auto const& chain = next_chains[nc];
		bool is_loop = (chain[0] == chain.back());
		auto const& stitches = next_stitches[nc];
		auto& ka = keep_adj[nc];
		ka.assign(stitches.size(), std::make_pair(true, true));
		for (uint32_t ns = 0; ns < stitches.size(); ++ns) {
			bool discard_before = false;
			bool discard_after = false;
			//This first ("easy") case is needed because of the case:
			// n0 --- x --- x --- n1
			//  \    /       \   /
			//    a0  -------  a1
			// which should get marked as
			// n0 xxx  xxx  xxx n1
			//   \             /
			//    a0  ------ a1
			// (otherwise, could not prune links_in -> links and probably skip this case.)
			if (stitches[ns].flag == Stitch::FlagDiscard) {
				discard_before = true;
				discard_after = true;
			}
			//stitches with a link to an active stitch whose next stitch doesn't have a link are marked discard-adj:
			discard_before = discard_before || [&]() -> bool {
				//okay, which active stitch is linked to this?
				auto fa = next_active.find(ChainStitch(nc, ns));
				if (fa == next_active.end()) return false; //nothing? no reason to discard before (though weird, I guess)
				//something? take earlier something:
				ChainStitch a = fa->second[0];
				auto const& a_chain = active_chains[a.chain];
				bool a_is_loop = (a_chain[0] == a_chain.back());
				auto const& a_stitches = active_stitches[a.chain];
				assert(a.stitch < a_stitches.size());
				{ //check for the case where a has an earlier non-discard link:
					auto fn = active_next.find(a);
					assert(fn != active_next.end());
					assert(fn->second.size() == 1 || fn->second.size() == 2);
					if (fn->second.size() == 2 && fn->second.back() == ChainStitch(nc, ns)) {
						// p -- n
						//  \  /
						//   a
						return false;
					}
					else {
						assert(fn->second[0] == ChainStitch(nc, ns));
					}
				}
				//check previous stitch:
				if (a.stitch == 0 && !a_is_loop) return false; //no previous stitch
				{
					ChainStitch pa(a.chain, (a.stitch > 0 ? a.stitch - 1 : a_stitches.size() - 1));
					auto fn = active_next.find(pa);
					if (fn == active_next.end()) {
						//        n  
						//   x    | /
						//  pa -- a
						return true; //aha! unlinked stitch; mark segment as not-keep.
					}
				}
				return false; //didn't have unlinked previous stitch
			}();
			discard_after = discard_after || [&]() -> bool {
				//okay, which active stitch is linked to this?
				auto fa = next_active.find(ChainStitch(nc, ns));
				if (fa == next_active.end()) return false; //nothing? no reason to discard after (though weird, I guess)
				//something? take later something:
				ChainStitch a = fa->second.back();
				auto const& a_chain = active_chains[a.chain];
				bool a_is_loop = (a_chain[0] == a_chain.back());
				auto const& a_stitches = active_stitches[a.chain];
				assert(a.stitch < a_stitches.size());
				{ //check for the case where a has a later non-discard link:
					auto fn = active_next.find(a);
					assert(fn != active_next.end());
					assert(fn->second.size() == 1 || fn->second.size() == 2);
					if (fn->second.size() == 2 && fn->second[0] == ChainStitch(nc, ns)) {
						// n -- nn
						//  \  /
						//   a
						return false;
					}
					else {
						assert((fn->second.size() == 1 && fn->second[0] == ChainStitch(nc, ns))
							|| (fn->second.size() == 2 && fn->second.back() == ChainStitch(nc, ns)));
					}
				}
				//check next stitch:
				if (a.stitch + 1 == a_stitches.size() && !a_is_loop) return false; //no next stitch
				{
					ChainStitch na(a.chain, (a.stitch + 1 < a_stitches.size() ? a.stitch + 1 : 0));
					auto fn = active_next.find(na);
					if (fn == active_next.end()) {
						//   n  
						// \ |     x
						//   a --- na
						return true; //aha! unlinked stitch; mark segment as not-keep.
					}
				}
				return false; //didn't have an unlinked next stitch
			}();

			if (discard_before) {
				if (ns > 0) ka[ns - 1].second = false;
				else if (is_loop) ka.back().second = false;
				ka[ns].first = false;
			}

			if (discard_after) {
				ka[ns].second = false;
				if (ns + 1 < stitches.size()) ka[ns + 1].first = false;
				else if (is_loop) ka[0].first = false;
			}
		}
	}



	//need lengths to figure out where stitches are on chains:
	//(duplicated, inelegantly, from link-chains)
	auto make_lengths = [&slice](std::vector< std::vector< uint32_t > > const& chains) {
		std::vector< std::vector< float > > all_lengths;
		all_lengths.reserve(chains.size());
		for (auto const& chain : chains) {
			all_lengths.emplace_back();
			std::vector< float >& lengths = all_lengths.back();
			lengths.reserve(chain.size());
			float total_length = 0.0f;
			lengths.emplace_back(total_length);
			for (uint32_t vi = 1; vi < chain.size(); ++vi) {
				glm::vec3 const& a = slice.vertices[chain[vi - 1]];
				glm::vec3 const& b = slice.vertices[chain[vi]];
				total_length += glm::length(b - a);
				lengths.emplace_back(total_length);
			}
		}
		return all_lengths;
	};
	std::vector< std::vector< float > > active_lengths = make_lengths(active_chains);
	std::vector< std::vector< float > > next_lengths = make_lengths(next_chains);

	std::vector< std::vector< uint32_t > > next_vertices;
	next_vertices.reserve(next_chains.size());

	if (!graph_) {
		//blank vertex info if no graph:
		for (auto const& stitches : next_stitches) {
			next_vertices.emplace_back(stitches.size(), -1U);
		}
	}
	else {
		assert(graph_);
		//for non-discard stitches on next chains, write down vertices:
		for (uint32_t nc = 0; nc < next_chains.size(); ++nc) {
			auto const& chain = next_chains[nc];
			bool is_loop = (chain[0] == chain.back());
			auto const& stitches = next_stitches[nc];
			auto const& ka = keep_adj[nc];

			auto const& lengths = next_lengths[nc];

			next_vertices.emplace_back();
			auto& vertices = next_vertices.back();
			vertices.reserve(stitches.size());

			auto li = lengths.begin();
			for (uint32_t ns = 0; ns < stitches.size(); ++ns) {
				assert(stitches[ns].vertex == -1U);
				if (stitches[ns].flag == Stitch::FlagDiscard) {
					vertices.emplace_back(-1U);
					continue;
				}
				else {
					float l = lengths.back() * stitches[ns].t;

					while (li != lengths.end() && *li <= l) ++li;
					assert(li != lengths.begin());
					assert(li != lengths.end());
					float m = (l - *(li - 1)) / (*li - *(li - 1));
					uint32_t i = li - lengths.begin();

					vertices.emplace_back(graph_->vertices.size());
					graph_->vertices.emplace_back();
					graph_->vertices.back().at = EmbeddedVertex::mix(
						slice_on_model[chain[i - 1]], slice_on_model[chain[i]], m
					);
				}
			}
			assert(ka.size() == stitches.size());
			assert(vertices.size() == stitches.size());
			if (!stitches.empty()) {
				uint32_t prev = (is_loop ? vertices.back() : -1U);
				for (uint32_t ns = 0; ns < stitches.size(); ++ns) {
					uint32_t cur = vertices[ns];
					if (ka[ns].first && prev != -1U) {
						assert(cur != -1U);
						//assert(prev != -1U); //<-- have to add to condition above because previous can be -1U in chains
						assert(cur < graph_->vertices.size() && prev < graph_->vertices.size());
						assert(graph_->vertices[cur].row_in == -1U);
						graph_->vertices[cur].row_in = prev;
						assert(graph_->vertices[prev].row_out == -1U);
						graph_->vertices[prev].row_out = cur;
					}
					prev = cur;
				}
			}
		}

		for (auto const& l : links) {
			uint32_t fv = active_stitches[l.from_chain][l.from_stitch].vertex;
			uint32_t tv = next_vertices[l.to_chain][l.to_stitch];
			assert(fv < graph_->vertices.size());
			assert(tv < graph_->vertices.size());
			graph_->vertices[fv].add_col_out(tv);
			graph_->vertices[tv].add_col_in(fv);
		}

	}

	//build a lookup structure for stitches:

	std::map< std::pair< OnChainStitch, OnChainStitch >, OnChainStitch > next_vertex;

	//edges where middle vertex is on a next chain:
	for (uint32_t nc = 0; nc < next_chains.size(); ++nc) {
		auto const& chain = next_chains[nc];
		bool is_loop = (chain[0] == chain.back());
		auto const& ak = keep_adj[nc];
		auto const& stitches = next_stitches[nc];

		assert(ak.size() == stitches.size());

		for (uint32_t ns = 0; ns < stitches.size(); ++ns) {
			if (stitches[ns].flag == Stitch::FlagDiscard) continue;
			OnChainStitch cur_ocs(OnChainStitch::OnNext, nc, ns);
			OnChainStitch prev_ocs;
			if (ak[ns].first) {
				//keep before segment, so prev is previous stitch (or to the begin of a non-loop):
				if (ns > 0 || is_loop) {
					prev_ocs.on = OnChainStitch::OnNext;
					prev_ocs.chain = nc;
					prev_ocs.stitch = (ns > 0 ? ns - 1 : stitches.size() - 1);
				}
				else {
					prev_ocs.on = OnChainStitch::OnNext;
					prev_ocs.chain = nc;
					prev_ocs.stitch = -1U;
					prev_ocs.type = OnChainStitch::TypeBegin;
				}
			}
			else {
				//don't keep before segment, so prev involves walking down a link:
				auto fa = next_active.find(ChainStitch(nc, ns));
				assert(fa != next_active.end()); //there should always be a link to walk down, right?
				prev_ocs.on = OnChainStitch::OnActive;
				prev_ocs.chain = fa->second[0].chain;
				prev_ocs.stitch = fa->second[0].stitch;
			}
			OnChainStitch next_ocs;
			if (ak[ns].second) {
				//keep after segment, so next is next stitch (or to the end of a non-loop):
				if (ns + 1 < stitches.size() || is_loop) {
					next_ocs.on = OnChainStitch::OnNext;
					next_ocs.chain = nc;
					next_ocs.stitch = (ns + 1 < stitches.size() ? ns + 1 : 0);
				}
				else {
					next_ocs.on = OnChainStitch::OnNext;
					next_ocs.chain = nc;
					next_ocs.stitch = -1U;
					next_ocs.type = OnChainStitch::TypeEnd;
				}
			}
			else {
				//don't keep next segment, so next involves walking down a link:
				auto fa = next_active.find(ChainStitch(nc, ns));
				assert(fa != next_active.end()); //there should always be a link to walk down, right?
				next_ocs.on = OnChainStitch::OnActive;
				next_ocs.chain = fa->second.back().chain;
				next_ocs.stitch = fa->second.back().stitch;
			}

			if (prev_ocs.on != OnChainStitch::OnNone && next_ocs.on != OnChainStitch::OnNone) {
				auto ret = next_vertex.insert(std::make_pair(
					std::make_pair(prev_ocs, cur_ocs),
					next_ocs
				));
				assert(ret.second);
			}
		}
	}

	//now edges from active chains:
	for (uint32_t ac = 0; ac < active_chains.size(); ++ac) {
		if (discard_active[ac]) {
			qDebug() << "Will discard active chain of " << active_stitches[ac].size() << " stitches because it had no outgoing links.";
			continue;
		}

		auto const& chain = active_chains[ac];
		bool is_loop = (chain[0] == chain.back());
		auto const& stitches = active_stitches[ac];

		for (uint32_t as = 0; as < stitches.size(); ++as) {
			OnChainStitch cur_ocs(OnChainStitch::OnActive, ac, as);

			OnChainStitch prev_ocs;
			prev_ocs.on = OnChainStitch::OnActive;
			prev_ocs.chain = ac;
			if (as > 0 || is_loop) {
				prev_ocs.stitch = (as > 0 ? as - 1 : stitches.size() - 1);
			}
			else {
				prev_ocs.type = OnChainStitch::TypeBegin;
			}
			OnChainStitch next_ocs;
			next_ocs.on = OnChainStitch::OnActive;
			next_ocs.chain = ac;
			if (as + 1 < stitches.size() || is_loop) {
				next_ocs.stitch = (as + 1 < stitches.size() ? as + 1 : 0);
			}
			else {
				next_ocs.type = OnChainStitch::TypeEnd;
			}

			auto fn = active_next.find(ChainStitch(ac, as));
			if (fn == active_next.end()) {
				//no link, so edge is previous to next:
				//      x    
				// p -> c -> n
				auto ret = next_vertex.insert(std::make_pair(
					std::make_pair(prev_ocs, cur_ocs),
					next_ocs
				));
				assert(ret.second);
			}
			else {
				//have a link. check for discarded segments:

				//previous segment is discarded, so link prev up:
				if (!keep_adj[fn->second[0].chain][fn->second[0].stitch].first) {
					auto n = fn->second[0];
					auto fa = next_active.find(n);
					assert(fa != next_active.end());
					if (fa->second[0] == ChainStitch(ac, as)) {
						//   xxx n0
						//       | /
						// p --- c
						auto ret = next_vertex.insert(std::make_pair(
							std::make_pair(prev_ocs, cur_ocs),
							OnChainStitch(OnChainStitch::OnNext, n.chain, n.stitch)
						));
						assert(ret.second);
					}
					else {
						assert(fa->second.size() == 2 && fa->second.back() == ChainStitch(ac, as));
						//don't link (this sort of case):
						//   xxx n0
						//     /  | 
						//   p -- c
					}
				}

				//next segment is discarded, so link down to next:
				if (!keep_adj[fn->second.back().chain][fn->second.back().stitch].second) {
					auto n = fn->second.back();
					auto fa = next_active.find(n);
					assert(fa != next_active.end());
					if (fa->second.back() == ChainStitch(ac, as)) {
						//   n0 xxx
						// \ |
						//   c --- n
						auto ret = next_vertex.insert(std::make_pair(
							std::make_pair(OnChainStitch(OnChainStitch::OnNext, n.chain, n.stitch), cur_ocs),
							next_ocs
						));
						assert(ret.second);
					}
					else {
						assert(fa->second.size() == 2 && fa->second[0] == ChainStitch(ac, as));
						//don't link (this sort of case):
						//    n0 xxx
						//   / |
						// c - n
					}
				}
			}
		}
	}


	//Walk through created edges array, creating chains therefrom:

	std::vector< std::vector< OnChainStitch > > loops;
	std::map< std::pair< OnChainStitch, OnChainStitch >, std::vector< OnChainStitch > > partials;

	while (!next_vertex.empty()) {
		std::vector< OnChainStitch > chain;
		chain.emplace_back(next_vertex.begin()->first.first);
		chain.emplace_back(next_vertex.begin()->first.second);
		chain.emplace_back(next_vertex.begin()->second);
		
		next_vertex.erase(next_vertex.begin());
		while (true) {
			auto f = next_vertex.find(std::make_pair(chain[chain.size() - 2], chain[chain.size() - 1]));
			if (f == next_vertex.end()) break;
			chain.emplace_back(f->second);
			next_vertex.erase(f);
		}

		{ //check if a partial chain comes after this one; if so, append it:
			auto f = partials.find(std::make_pair(chain[chain.size() - 2], chain[chain.size() - 1]));
			if (f != partials.end()) {
				chain.pop_back();
				chain.pop_back();
				chain.insert(chain.end(), f->second.begin(), f->second.end());
				partials.erase(f);
			}
		}

		//loops should look like abcdab
		//because abcd -> ab-c bc-d cd-a da-b
		if (chain[0] == chain[chain.size() - 2] && chain[1] == chain[chain.size() - 1]) {
			//great -- full loop.
			chain.pop_back();
			loops.emplace_back(chain);
		}
		else {
			//partial loop -- save for later
			auto ret = partials.insert(std::make_pair(std::make_pair(chain[0], chain[1]), chain));
			assert(ret.second);
		}
	}

	auto output = [&](std::vector< OnChainStitch > const& path) {
		assert(path.size() >= 2);

		{ //don't pass onward any chains that touch a boundary:
			uint32_t on_discard = 0;
			uint32_t on_non_discard = 0;
			for (auto const& ocs : path) {
				if (ocs.on == OnChainStitch::OnNext && next_used_boundary[ocs.chain]) {
					++on_discard;
				}
				else {
					++on_non_discard;
				}
			}
			if (on_discard) {
				assert(on_non_discard == 0); //should either be entirely discard or entirely not!
				return;
			}
		}

		//build embedded vertices for all stitches:
		std::vector< EmbeddedVertex > path_evs;
		std::vector< uint32_t > path_lefts; //also loop up i such that the stitch is in [ chain[i], chain[l+1] )
		path_evs.reserve(path.size());
		path_lefts.reserve(path.size());
		for (auto const& ocs : path) {
			std::vector< uint32_t > const& src_chain = (ocs.on == OnChainStitch::OnActive ? active_chains : next_chains)[ocs.chain];
			std::vector< float > const& src_lengths = (ocs.on == OnChainStitch::OnActive ? active_lengths : next_lengths)[ocs.chain];
			std::vector< Stitch > const& src_stitches = (ocs.on == OnChainStitch::OnActive ? active_stitches : next_stitches)[ocs.chain];
			assert(!src_chain.empty());
			assert(src_lengths.size() == src_chain.size());
			
			if (ocs.type == OnChainStitch::TypeBegin) {
				assert(src_chain[0] != src_chain.back());
				path_evs.emplace_back(EmbeddedVertex::on_vertex(src_chain[0]));
				path_lefts.emplace_back(0);
			}
			else if (ocs.type == OnChainStitch::TypeEnd) {
				assert(src_chain[0] != src_chain.back());
				path_evs.emplace_back(EmbeddedVertex::on_vertex(src_chain.back()));
				path_lefts.emplace_back(src_chain.size() - 2);
			}
			else {
				assert(ocs.type == OnChainStitch::TypeStitch);
				assert(ocs.stitch < src_stitches.size());
				float l = src_lengths.back() * src_stitches[ocs.stitch].t;
				auto li = std::upper_bound(src_lengths.begin(), src_lengths.end(), l);
				assert(li != src_lengths.end());
				assert(li != src_lengths.begin());
				float m = (l - *(li - 1)) / (*li - *(li - 1));
				uint32_t i = li - src_lengths.begin();
				assert(i > 0);
				path_evs.emplace_back(EmbeddedVertex::on_edge(src_chain[i - 1], src_chain[i], m));
				path_lefts.emplace_back(i - 1);
			}
		}

		std::vector< EmbeddedVertex > chain;
		std::vector< Stitch > stitches;
		std::vector< uint32_t > remove_stitches; //indices of stitches to remove in a moment.
		float length = 0.0f;

		auto append_ev = [&chain, &slice, &length](EmbeddedVertex const& ev, char const* why) {
			if (!chain.empty()) {
				assert(ev != chain.back());
				EmbeddedVertex::common_simplex(chain.back().simplex, ev.simplex); //make sure this works
				length += glm::length(
					chain.back().interpolate(slice.vertices) - ev.interpolate(slice.vertices)
				);
			}
			chain.emplace_back(ev);
		};

		for (uint32_t pi = 0; pi + 1 < path.size(); ++pi) {
			auto const& a = path[pi];
			auto const& a_ev = path_evs[pi];
			auto const& a_left = path_lefts[pi];
			std::vector< uint32_t > const& a_chain = (a.on == OnChainStitch::OnActive ? active_chains : next_chains).at(a.chain);
			assert(a_left + 1 < a_chain.size());
			auto const& b = path[pi + 1];
			auto const& b_ev = path_evs[pi + 1];
			auto const& b_left = path_lefts[pi + 1];
			std::vector< uint32_t > const& b_chain = (b.on == OnChainStitch::OnActive ? active_chains : next_chains).at(b.chain);
			assert(b_left + 1 < b_chain.size());

			
			if (pi == 0) append_ev(a_ev, "first a");
			else assert(!chain.empty() && chain.back() == a_ev);
			if (a.type == OnChainStitch::TypeBegin) {
				assert(pi == 0);
				assert(a_left == 0);
			}
			else {
				assert(a.type == OnChainStitch::TypeStitch);
				std::vector< Stitch > const& a_stitches = (a.on == OnChainStitch::OnActive ? active_stitches : next_stitches).at(a.chain);
				assert(a.stitch < a_stitches.size());
				Stitch::Flag a_flag = a_stitches[a.stitch].flag;
				uint32_t a_vertex;
				if (a.on == OnChainStitch::OnActive) {
					a_vertex = active_stitches.at(a.chain).at(a.stitch).vertex;
				}
				else {
					a_vertex = next_vertices.at(a.chain).at(a.stitch);
				}
				if (pi == 0) stitches.emplace_back(length, a_flag, a_vertex);
				else assert(!stitches.empty() && stitches.back().t == length && stitches.back().flag == a_flag && stitches.back().vertex == a_vertex);

				if (a.on != b.on && a.on == OnChainStitch::OnActive) {
					remove_stitches.emplace_back(stitches.size() - 1);
				}
			}

			if (a.on == b.on) {
				assert(a.chain == b.chain);
				bool is_loop = (a_chain[0] == a_chain.back());
				if (a_left != b_left) {
					assert(b_left + 1 < a_chain.size()); //RIIIIIGHT?
					uint32_t v = a_left;
					do {
						//advance v
						assert(v + 1 < a_chain.size());
						v += 1;
						if (v + 1 == a_chain.size()) {
							assert(is_loop);
							v = 0;
						}
						append_ev(EmbeddedVertex::on_vertex(a_chain[v]), "lefts");
					} while (v != b_left);
				}
			}
			else {
				//find an embedded path between a and b:
				std::vector< EmbeddedVertex > ab;
				embeddedPath(slice, a_ev, b_ev, &ab);
				assert(ab[0] == a_ev);
				assert(ab.back() == b_ev);
				for (uint32_t i = 1; i + 1 < ab.size(); ++i) {
					append_ev(ab[i], "embedded path");
				}
			}

			append_ev(b_ev, "b");
			if (b.type == OnChainStitch::TypeEnd) {
				assert(pi + 2 == path.size());
				assert(b_left == b_chain.size() - 2);
			}
			else {
				assert(b.type == OnChainStitch::TypeStitch);
				std::vector< Stitch > const& b_stitches = (b.on == OnChainStitch::OnActive ? active_stitches : next_stitches).at(b.chain);
				assert(b.stitch < b_stitches.size());
				Stitch::Flag b_flag = b_stitches[b.stitch].flag;
				uint32_t b_vertex;
				if (b.on == OnChainStitch::OnActive) {
					b_vertex = active_stitches.at(b.chain).at(b.stitch).vertex;
				}
				else {
					b_vertex = next_vertices.at(b.chain).at(b.stitch);
				}

				stitches.emplace_back(length, b_flag, b_vertex);

				if (a.on != b.on && b.on == OnChainStitch::OnActive) {
					remove_stitches.emplace_back(stitches.size() - 1);
				}
			}
		}

		//should turn loops into loops:
		assert((chain[0] == chain.back()) == (path[0] == path.back()));

		//remove last stitch from loops:
		if (path[0] == path.back()) {
			assert(!stitches.empty());
			assert(stitches[0].flag == stitches.back().flag);
			assert(stitches[0].vertex == stitches.back().vertex);
			assert(stitches[0].t == 0.0f && stitches.back().t == length);
			//if last was tagged for removal, well, tag first instead:
			if (!remove_stitches.empty() && remove_stitches.back() == stitches.size() - 1) {
				remove_stitches.back() = 0;
			}
			stitches.pop_back();
		}

		{ //remove any stitches marked for discard:
			for (auto s : remove_stitches) {
				assert(s < stitches.size());
				stitches[s].flag = Stitch::FlagDiscard;
			}
			auto out = stitches.begin();
			for (auto in = stitches.begin(); in != stitches.end(); ++in) {
				assert(out <= in);
				if (in->flag != Stitch::FlagDiscard) {
					*(out++) = *in;
				}
			}
			if (out != stitches.end()) {
				qDebug() << "Removed " << stitches.end() - out << " stitches under short-row ends.";
				stitches.erase(out, stitches.end());
			}
		}

		//convert stitch t values from lengths to [0,1) range:
		for (auto& s : stitches) {
			s.t /= length;
		}

		//convert chain from being embedded on slice to being embedded on model:
		for (auto& ev : chain) {
			glm::uvec3 simplex = slice_on_model[ev.simplex.x].simplex;
			if (ev.simplex.y != -1U) simplex = EmbeddedVertex::common_simplex(simplex, slice_on_model[ev.simplex.y].simplex);
			if (ev.simplex.z != -1U) simplex = EmbeddedVertex::common_simplex(simplex, slice_on_model[ev.simplex.z].simplex);

			glm::vec3 weights = ev.weights.x * slice_on_model[ev.simplex.x].weights_on(simplex);
			if (ev.simplex.y != -1U) weights += ev.weights.y * slice_on_model[ev.simplex.y].weights_on(simplex);
			if (ev.simplex.z != -1U) weights += ev.weights.z * slice_on_model[ev.simplex.z].weights_on(simplex);

			ev = EmbeddedVertex::canonicalize(simplex, weights);
		}

		next_active_chains.emplace_back(chain);
		next_active_stitches.emplace_back(stitches);

	};

	qDebug() << "Found " << loops.size() << " loops and " << partials.size() << " chains.";
	for (auto const& loop : loops) {
		output(loop);
	}
	for (auto const& pp : partials) {
		output(pp.second);
	}

	//HACK: sometimes duplicate vertices after splatting back to model somehow:
	uint32_t trimmed = 0;
	for (auto& chain : next_active_chains) {
		for (uint32_t i = 1; i < chain.size(); /* later */) {
			if (chain[i - 1] == chain[i]) {
				chain.erase(chain.begin() + i);
				++trimmed;
			}
			else {
				++i;
			}
		}
	}

	if (trimmed) {
		qDebug() << "Trimmed " << trimmed << " identical-after-moving-to-model vertices from next active chains.";
	}


	//PARANOIA:
	assert(next_active_stitches.size() == next_active_chains.size());
	for (auto const& chain : next_active_chains) {
		for (uint32_t i = 1; i < chain.size(); ++i) {
			assert(chain[i - 1] != chain[i]);
		}
	}

	//PARANOIA:
	if (graph_) {
		for (auto const& stitches : next_active_stitches) {
			for (auto const& s : stitches) {
				assert(s.vertex != -1U);
				assert(s.vertex < graph_->vertices.size());
			}
		}
	}
	

}


/*
*	Function: cleaning function that is called after an iteration of the 4 main steps (find, slice, link, build) is completed
*	Return: void, but clears a lot of class member variables
*	Called: called before the find step of the algorithm, is called when the user clicks the 'step' button on the GUI at least 4 times
*/
void KnitGrapher::clearPeeling() {
	graph.clear();
	active_chains.clear();
	active_stitches.clear();
	slice.clear();
	sliceOnModel.clear();
	sliceActiveChains.clear();
	sliceNextChains.clear();
	nextUsedBoundary.clear();
	sliceTimes.clear();
	nextStitches.clear();
	links.clear();
	nextActiveChains.clear();
	nextActiveStitches.clear();
}


/*
*	Function: a SLOT function that is warned that the user clicked the 'step' button, does a different action based on which step is current
*	Return: void, starts the different autknit algorithms
*	Called: is a SLOT function, called when the user clicks the 'step' button
*/
void KnitGrapher::stepButtonClicked()
{
	if (stepCount % 4 == 0) {
		auto old_next_active_chains = nextActiveChains;
		auto old_next_active_stitches = nextActiveStitches;
		auto old_rowcol_graph = graph;
		clearPeeling();
		
		graph = old_rowcol_graph;
		if (stepCount == 0) {
			qDebug() << "[Step 0] - peel begin, calling findFirstActiveChains()---------------------------------------------------";

			//generate the 3 variables used in findFirstActiveChains()
			findFirstActiveChains(&active_chains, &active_stitches, &graph);
		}
		else {
			qDebug() << "[Step " << stepCount << "] - repeat, moving active and next chains...--------------------------------------";
			active_chains = old_next_active_chains;
			active_stitches = old_next_active_stitches;
		}
		if (!active_chains.empty())
			emit firstActiveChainsCreated(&active_chains, &active_stitches, &graph);
		else {
			stepCount--;
			emit firstActiveChainsCreated(&active_chains, &active_stitches, &graph);
		}
	}
	else if (stepCount % 4 == 1) {
		qDebug() << "[Step " << stepCount << " ] - slice, calling peelSlice()-------------------------------------------------------";
		
		peelSlice(active_chains, &slice, &sliceOnModel, &sliceActiveChains, &sliceNextChains, &nextUsedBoundary);
		
		qDebug() << "sliceonmodel size" << sliceOnModel.size();
		sliceTimes.clear();
		sliceTimes.reserve(sliceOnModel.size());
		for (auto& ev : sliceOnModel) {
			sliceTimes.emplace_back(ev.interpolate(constrained_values));
		}
		
		ObjectMesh emittedMesh = slice.toObjMesh();
		emit peelSliceDone(&emittedMesh, &sliceActiveChains, &sliceNextChains);
	}
	else if (stepCount % 4 == 2) {
		qDebug() << "[Step "<< stepCount << " ] - link, calling linkChains()-------------------------------------------------------------------";
		
		qDebug() << "slice size" << sliceNextChains.size();
		linkChains(slice, sliceTimes, sliceActiveChains, active_stitches, sliceNextChains, nextUsedBoundary, &nextStitches, &links);
		emit linkChainsDone(&nextStitches, &links);
	}
	else if (stepCount % 4 == 3) {
		qDebug() << "[Step " << stepCount << " ] - build, calling buildNextActiveChains()------------------------------------------------------------";
		buildNextActiveChains(slice, sliceOnModel, sliceActiveChains, active_stitches, sliceNextChains, nextStitches, nextUsedBoundary, links, &nextActiveChains, &nextActiveStitches, &graph);
		emit nextActiveChainsDone(&nextActiveChains);
	}
	stepCount++;
}


/*
*	Function: construct the mesh that will be used in the autoknit algorithms 
*	Return: void, but creates the newMesh class variable
*	Called: called after the user selects the 'interpolate' button in the GUI
*/
void KnitGrapher::constructNewMesh(std::vector<Constraint*> constraints)
{
	qDebug() << "Constructing new mesh, sizes:" << stitchWidth << stitchHeight << modelUnitLength;
	//create a new mesh that conforms to the constraints and stitch size given by the user
	qDebug() << "Constraints:" << constraints.size();
	this->constraints = constraints;   
	
	remesh();
	interpolateValues();
	
	ObjectMesh emittedMesh = newMesh.toObjMesh();
	emit knitGraphInterpolated(emittedMesh, constrained_values);
}

/*
*	Function: the following 3 functions are SLOT functions that change the class parameters based on users input in the GUI
*/
void KnitGrapher::setStitchWidth(float width)
{
	qDebug() << "Setting stitch width to " << width;
	stitchWidth = width;
}
void KnitGrapher::setStitchHeight(float height)
{
	stitchHeight = height;
}
void KnitGrapher::setModelUnitLength(float length)
{
	modelUnitLength = length;
}