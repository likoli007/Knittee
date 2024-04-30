#pragma once
#include "KnitGrapher.h"
#include "EmbeddedVertex.h"
#include "EmbeddedPlanarMap.h"
#include <vector>
#include <QSet>
#include <QPair>
#include <set>
#include <deque>
#include <Eigen/SparseCholesky>

std::vector<KnitGrapher::Edge> KnitGrapher::globedges;

KnitGrapher::KnitGrapher(QObject* parent) : QObject(parent)
{
	stitchWidth = 5;
	stitchHeight = 5;
	modelUnitLength = 1;

}

void KnitGrapher::setOriginalMesh(ObjectMesh mesh)
{
	originalMesh = mesh;
}

QPair<GLuint, GLuint> qMinMax(GLuint a, GLuint b)
{
	if (a < b)
		return QPair<GLuint, GLuint>(a, b);
	else
		return QPair<GLuint, GLuint>(b, a);
}



float KnitGrapher::getMaxEdgeLength()
{
	return 0.5f * std::min(stitchWidth, 2.0f * stitchHeight) / modelUnitLength;
}


void KnitGrapher::generateTriangles() {
	//populate the oldTriangles vector with uvec3 by traversing my oldMesh.indices
	for (int i = 0; i < originalMesh.indices.size(); i += 3) {
		glm::uvec3 triangle(originalMesh.indices[i], originalMesh.indices[i + 1], originalMesh.indices[i+2]);
		oldTriangles.push_back(triangle);
	}
}

//custom comparator used in visit(), since QPair does no have <= == operators
auto customComparator = [](const QPair<float, GLuint>& lhs, const QPair<float, GLuint>& rhs) {
	if (lhs.first > rhs.first)
		return true;
	else if (lhs.first == rhs.first)
		return lhs.second > rhs.second;
	else
		return false;
};

// Dijkstra's algorithm 'visit' function
void visit(std::vector<QPair<float, GLuint>>& todo,
	std::vector<QPair<float, GLuint>>& visited,
	GLuint vertex,
	float distance,
	GLuint from) {
	if (distance < visited[vertex].first) {
		visited[vertex] = QPair<float, GLuint>(distance, from);
		todo.emplace_back(distance, vertex);
		std::push_heap(todo.begin(), todo.end(), customComparator);
	}
}

GLuint KnitGrapher::lookup(GLuint a, GLuint b, QHash<QPoint, GLuint>& marked_verts) {
	auto f = marked_verts.find((a < b ? QPoint(a, b) : QPoint(b, a)));
	if (f != marked_verts.end()) return f.value();
	else return -1;
}

// quad not even implemented in my mesh loading, but good to have this here for later use
void KnitGrapher::quad(std::vector<glm::uvec3>& new_tris, std::vector<QVector3D>& verts, GLuint a, GLuint b, GLuint c, GLuint d) {
	float ac = QVector3D(verts[c] - verts[a]).lengthSquared();
	float bd = QVector3D(verts[d] - verts[b]).lengthSquared();
	if (ac < bd) {
		new_tris.emplace_back(glm::uvec3(a, b, c));
		new_tris.emplace_back(glm::uvec3(c, d, a));
	}
	else {
		new_tris.emplace_back(glm::uvec3(a, b, d));
		new_tris.emplace_back(glm::uvec3(b, c, d));
	}
}


void KnitGrapher::divide(QSet<QPoint>& marked,
	std::vector<QVector3D>& verts,
	std::vector<std::vector<GLuint>>& paths,
	std::vector<glm::uvec3>& tris)
{
	//assert(!marked.empty());
	if (marked.empty()) {
		qDebug() << "Marked is empty!";
		return;
	}

	QHash<QPoint, GLuint> marked_verts;
	marked_verts.reserve(marked.size());

	std::vector<QPoint> edges(marked.begin(), marked.end());

	//lambda sort, had to change to QPoint 
	std::sort(edges.begin(), edges.end(), [](const QPoint& a, const QPoint& b) {
		if (a.x() != b.x()) return a.x() < b.x();
		else return a.y() < b.y();
		});

	for (auto const& e : edges) {
		marked_verts.insert(e, verts.size());
		// marked_verts.insert(std::make_pair(e, verts.size()));
		verts.emplace_back((verts[e.x()] + verts[e.y()]) / 2.0f);   //why 2.0f?
	}

	/*
	auto lookup_func = [&marked_verts](int a, int b) {
		return lookup(a, b, marked_verts);
	};*/

	for (auto& path : paths) {
		std::vector<GLuint> new_path;
		new_path.emplace_back(path[0]);
		for (int i = 1; i < path.size(); ++i) {
			int v = lookup(path[i - 1], path[i], marked_verts);
			if (v != -1U) new_path.emplace_back(v);
			new_path.emplace_back(path[i]);
		}
		path = std::move(new_path);
	}

	std::vector<glm::uvec3> new_tris;

	for (auto const& tri : tris) {
		GLuint a = tri.x;
		GLuint b = tri.y;
		GLuint c = tri.z;
		GLuint ab = lookup(a, b, marked_verts);
		GLuint bc = lookup(b, c, marked_verts);
		GLuint ca = lookup(c, a, marked_verts);

		if (ab != -1U && bc != -1U && ca != -1U) {
			//1 -> 4 subdiv!
			new_tris.emplace_back(glm::uvec3(a, ab, ca));
			new_tris.emplace_back(glm::uvec3(b, bc, ab));
			new_tris.emplace_back(glm::uvec3(c, ca, bc));
			new_tris.emplace_back(glm::uvec3(ab, bc, ca));
		}
		else if (ab != -1U && bc != -1U && ca == -1U) {
			//1 -> 3 subdiv!
			//NOTE: should consider recursively subdividing to avoid this case
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
			//assert(ab == -1U && bc == -1U && ca == -1U);
			//no subdiv!
			new_tris.emplace_back(glm::uvec3(a, b, c));
		}
	}
	tris = std::move(new_tris);
	newTriangles = new_tris;
}

bool is_ccw(QVector2D const& a, QVector2D const& b, QVector2D const& c) {
	return QVector2D::dotProduct(QVector2D(-(b.y() - a.y()), (b.x() - a.x())), c - a) > 0.0f;
};

inline uint qHash(const QVector2D& key)
{
	return qHash(QPair<float, float>(key.x(), key.y()));
}

float& get_dis(GLuint a, GLuint b, QHash<QPoint, float>& min_dis) {
	if (a > b) std::swap(a, b);
	auto it = min_dis.insert(QPoint(a, b), qInf());
	return it.value();
};
void KnitGrapher::unfold(GLuint depth, GLuint root, QVector2D const& flat_root,
	GLuint ai, QVector2D const& flat_a,
	GLuint bi, QVector2D const& flat_b,
	QVector2D const& limit_a, QVector2D const& limit_b, QHash< QPoint, GLuint > const& opposite,
	std::vector<QVector3D> const& newVertices, QHash<QPoint, float>& min_dis) {

	//std::cout << "r: " << root << ": (" << flat_root.x << ", " << flat_root.y << ")" << std::endl; //DEBUG
				//std::cout << "a: " << ai << ": (" << flat_a.x << ", " << flat_a.y << ")" << std::endl; //DEBUG
				//std::cout << "b: " << bi << ": (" << flat_b.x << ", " << flat_b.y << ")" << std::endl; //DEBUG
	//assert(is_ccw(flat_root, flat_a, flat_b));
	//should go 'a - limit_a - limit_b - b':
	//assert(flat_a == limit_a || is_ccw(flat_root, flat_a, limit_a));
	//assert(is_ccw(flat_root, limit_a, limit_b));
	//assert(flat_b == limit_b || is_ccw(flat_root, limit_b, flat_b));




	uint32_t ci;
	QVector2D flat_c;
	{ //if there is a triangle over the ai->bi edge, find other vertex and flatten it:
		auto f = opposite.find(QPoint(bi, ai));
		if (f == opposite.end()) return;
		ci = f.value();
		//figure out c's position along ab and distance from ab:
		QVector3D const& a = newVertices[ai];
		QVector3D const& b = newVertices[bi];
		QVector3D const& c = newVertices[ci];

		QVector3D ab = (b - a).normalized();
		float along = QVector3D::dotProduct(c - a, ab);
		float perp = -QVector3D::crossProduct(c - a, ab).length();

		QVector2D flat_ab = (flat_b - flat_a).normalized();
		QVector2D flat_perp_ab = QVector2D(-flat_ab.y(), flat_ab.x());

		flat_c = flat_a + flat_ab * along + flat_perp_ab * perp;
	}

	//std::cout << "c: " << ci << ": (" << flat_c.x << ", " << flat_c.y << ")" << std::endl; //DEBUG

	//flat_a and flat_b should always be outside limit, it seems like we need to test anyway (thanks, numerics)

	bool ccw_rac = is_ccw(flat_root, limit_a, flat_c) && is_ccw(flat_root, flat_a, flat_c);
	bool ccw_rcb = is_ccw(flat_root, flat_c, limit_b) && is_ccw(flat_root, flat_c, flat_b);

	if (ccw_rac && ccw_rcb) {
		float& dis = get_dis(root, ci, min_dis);
		dis = std::min(dis, (flat_root - flat_c).length());

		//PARANOIA:
		float dis3 = (newVertices[root] - newVertices[ci]).length();


		/*
		if (dis3 > dis * (1.0f + 1e-6f) + 1e-6f) {
			qDebug() << "dis3: " << dis3 << " vs flat dis " << dis << " seems bad!" << std::endl;
			qDebug() << "  ra3: " << glm::length(verts[root] - verts[ai]) << " vs ra: " << glm::length(flat_root - flat_a) << std::endl;
			qDebug() << "  rb3: " << glm::length(verts[root] - verts[bi]) << " vs rb: " << glm::length(flat_root - flat_b) << std::endl;
			qDebug() << "  ab3: " << glm::length(verts[ai] - verts[bi]) << " vs ab: " << glm::length(flat_a - flat_b) << std::endl;
			qDebug() << "  ac3: " << glm::length(verts[ai] - verts[ci]) << " vs ac: " << glm::length(flat_a - flat_c) << std::endl;
			qDebug() << "  bc3: " << glm::length(verts[bi] - verts[ci]) << " vs bc: " << glm::length(flat_b - flat_c) << std::endl;
			assert(dis3 < dis + 1e-6);
		}*/

		if (depth > 1) {
			//assert(is_ccw(flat_root, flat_a, flat_c));
			unfold(depth - 1, root, flat_root, ai, flat_a, ci, flat_c, limit_a, flat_c, opposite, newVertices, min_dis);
			//assert(is_ccw(flat_root, flat_c, flat_b));
			unfold(depth - 1, root, flat_root, ci, flat_c, bi, flat_b, flat_c, limit_b, opposite, newVertices, min_dis);
		}
	}
	else if (ccw_rac && !ccw_rcb) {
		if (depth > 1) {
			//assert(!is_ccw(flat_root, flat_c, limit_b)); //DEBUG
			//assert(is_ccw(flat_root, limit_b, flat_c)); //DEBUG -- fails sometimes [thanks, numerics]
			//assert(is_ccw(flat_root, flat_a, flat_c));
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


void KnitGrapher::remesh()
{
	qDebug() << "Start of remesh() function...";

	//check if constraints are empty, no way to remesh without them
	if (constraints.size() == 0) {
		qDebug() << "No Constraints found! exiting remesh()...";
		return;
	}

	//extract edges from the model
	qDebug() << "extracting edges from all triangles...";
	std::vector< std::vector< QPair< GLuint, float > > > adjacents(originalMesh.vertices.size());
	QSet <QPair<GLuint, GLuint>> edges;
	{ 
		for (auto const& tri : oldTriangles) {
			edges.insert(qMinMax(tri.x, tri.y));
			edges.insert(qMinMax(tri.y, tri.z));
			edges.insert(qMinMax(tri.z, tri.x));
		}
		for (QPair<GLuint, GLuint> edge : edges)
		{
			float length = (originalMesh.vertices[edge.second] - originalMesh.vertices[edge.first]).length();

			QPair <GLuint, float> edge1(edge.second, length);
			QPair <GLuint, float> edge2(edge.first, length);

			adjacents[edge.first].push_back(edge1);
			adjacents[edge.second].push_back(edge2);
		}
	}


	//find chain paths on original model
	qDebug() << "finding chain paths on original model... (linking constraint vertices together)";
	std::vector<std::vector<GLuint>> paths;

	//altough done here, technically my constraints chain has nice 'single jump' vertex chains, so not really necessary
	//could be coopted in Visualizer for constraint by dijkstra adding...
	for (Constraint* constraint : constraints)
	{
		if (constraint->vertices.empty()) {
			qDebug() << "Constraint has no vertices, skipping...";
			continue;
		}
		std::vector<GLuint> path;
		for (GLuint goal : constraint->vertices) {
			if (path.empty()) {
				path.push_back(goal);
				continue;
			}
			std::vector<QPair<float, GLuint>> todo;
			std::vector<QPair<float, GLuint>> visited(originalMesh.vertices.size(), QPair<float, GLuint>(qInf(), -1U));
			visit(todo, visited, goal, 0.0f, -1U);

			while (!todo.empty()) {
				std::pop_heap(todo.begin(), todo.end(), customComparator);
				QPair<float, GLuint> at = todo.back();
				todo.pop_back();
				if (at.first > visited[at.second].first) continue;
				if (at.second == path.back()) break;
				for (QPair<float, GLuint> const& a : adjacents[at.second]) {
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

	qDebug() << "Done!";
	qDebug() << "Paths:" << paths.size();

	
	for (std::vector<GLuint> path : paths)
	{
		qDebug() << "Path:";
		for (int vertex : path)
		{
			qDebug() << vertex;
		}
	}

	// create a higher resolution mesh according to the parameters set by the user
	// in algorithm is part of the 'Parameters' class, i have those values inside this KnitGrapher class
	float maxEdgeLength = getMaxEdgeLength();
	float maxEdgeLengthSquared = maxEdgeLength * maxEdgeLength;


	qDebug() << "Max Edge Length:" << maxEdgeLength;

	std::vector<QVector3D> newVerts = originalMesh.vertices;
	std::vector<glm::uvec3> newTris = oldTriangles;

	// there is degenarate triangle checking in the original algorithm, will do the same here
	if (degenerateCheck(newTris)) {
		qDebug() << "Degenerate triangle found, error...";
		return;
	}

	qDebug() << "Starting divide() block, oldTris: " << newTris.size();
	// Lambda functions in original code were made into their own functions...
	while (true) {
		QSet<QPoint> marked;
		auto mark = [&marked](int a, int b) {
			if (b < a) std::swap(a, b);
			marked.insert(QPoint(a, b));
		};
		auto is_marked = [&marked](int a, int b) {
			if (b < a) std::swap(a, b);
			return marked.find(QPoint(a, b)) != marked.end();
		};
		(void)is_marked;
		(void)minEdgeRatioSquared;

		for (auto const& tri : newTris) {
			float len_ab2 = QVector3D(newVerts[tri.y] - newVerts[tri.x]).lengthSquared();
			float len_bc2 = QVector3D(newVerts[tri.z] - newVerts[tri.y]).lengthSquared();
			float len_ca2 = QVector3D(newVerts[tri.x] - newVerts[tri.z]).lengthSquared();

			//qDebug() << "Lengths: " << len_ab2 << len_bc2 << len_ca2;

			if (len_ab2 > maxEdgeLengthSquared) mark(tri.x, tri.y);
			if (len_bc2 > maxEdgeLengthSquared) mark(tri.y, tri.z);
			if (len_ca2 > maxEdgeLengthSquared) mark(tri.z, tri.x);
		}
		if (marked.empty()) {
			qDebug() << "No marked edges found, breaking...";
			break;
		}
		divide(marked, newVerts, paths, newTris);
	}
	qDebug() << "divide() sleu finished! newTris: " << newTris.size();

	//after this while cycle, have newTriangles variable equal to divides newTris...
	//once again the degenerate triangle check in original algortihm
	//if (degenerateCheck(newTris)) {
		//qDebug() << "Degenerate triangle found after remesh! skipping...";
		//return;
	//}


	//extract edges from subdivided model:
	qDebug() << "extracting edges from subdivided model...";
	adjacents.assign(newVerts.size(), std::vector< QPair< GLuint, float > >());
	{ //extract edges from subdivided model:
		QSet<QPair<GLuint, GLuint>> newEdges;
		for (auto const& tri : newTris) {
			edges.insert(qMinMax(tri.x, tri.y));
			edges.insert(qMinMax(tri.y, tri.z));
			edges.insert(qMinMax(tri.z, tri.x));
		}

		for (auto const& e : edges) {
			float len = QVector3D(newVerts[e.second] - newVerts[e.first]).length();
			//float len = glm::length(verts[e.second] - verts[e.first]);
			adjacents[e.first].emplace_back(e.second, len);
			adjacents[e.second].emplace_back(e.first, len);
		}
	}


	qDebug() << "Building opposites...";
	QHash< QPoint, GLuint > opposite; //vertex opposite each [oriented] triangle edge
	opposite.reserve(newTris.size());
	for (auto const& tri : newTris) {
		auto ret_xy = opposite.insert(QPoint(tri.x, tri.y), tri.z);
		//assert(ret_xy.second);
		auto ret_yz = opposite.insert(QPoint(tri.y, tri.z), tri.x);
		//assert(ret_yz.second);
		auto ret_zx = opposite.insert(QPoint(tri.z, tri.x), tri.y);
		//assert(ret_zx.second);
	}

	qDebug() << "build (+ add to adj) extra shortcut edges by unwrapping triangle neighborhoods:";
	{ //build (+ add to adj) extra "shortcut" edges by unwrapping triangle neighborhoods:
		QHash< QPoint, float > min_dis;
		for (auto const& tri : newTris) {
			QVector2D flat_x, flat_y, flat_z; //original verts
			QVector3D const& x = newVerts[tri.x];
			QVector3D const& y = newVerts[tri.y];
			QVector3D const& z = newVerts[tri.z];
			flat_x = QVector2D(0.0f, 0.0f);
			flat_y = QVector2D((y - x).length(), 0.0f);

			QVector3D xy = (y - x).normalized();
			QVector3D perp_xy = QVector3D::crossProduct(QVector3D::crossProduct(y - x, z - x), y - x).normalized();
			float along = QVector3D::dotProduct(z - x, xy);
			float perp = QVector3D::dotProduct(z - x, perp_xy);

			flat_z = QVector2D(along, perp);

			//look through edge [ai,bi] from point [root], where edge [ai,bi] is ccw oriented.
			const constexpr GLuint D = 3; //depth to unfold triangles to for more adjacency information; makes slightly nicer geodesics at the expense of increased compute time.
			//THIS COULD BE A SETTABLE PARAMETER!
			//TODO: make this a parameter

			if (D > 0) {
				unfold(D, tri.x, flat_x, tri.y, flat_y, tri.z, flat_z, flat_y, flat_z, opposite, newVerts, min_dis);
				unfold(D, tri.y, flat_y, tri.z, flat_z, tri.x, flat_x, flat_z, flat_x, opposite, newVerts, min_dis);
				unfold(D, tri.z, flat_z, tri.x, flat_x, tri.y, flat_y, flat_x, flat_y, opposite, newVerts, min_dis);
			}
		}
		for (uint32_t x = 0; x < newVerts.size(); ++x) {
			for (auto const& yd : adjacents[x]) {
				float& dis = get_dis(x, yd.first, min_dis);
				dis = std::min(dis, yd.second);
			}
		}

		//clear adj + re-create from min_dis:
		uint32_t old_adj = 0;
		for (auto const& a : adjacents) {
			old_adj += a.size();
		}

		adjacents.assign(newVerts.size(), std::vector< QPair< GLuint, float > >());

		QHash<QPoint, float>::const_iterator xyd;
		for (xyd = min_dis.constBegin(); xyd != min_dis.constEnd(); ++xyd) {
			adjacents[xyd.key().x()].emplace_back(xyd.key().y(), xyd.value());
			adjacents[xyd.key().y()].emplace_back(xyd.key().x(), xyd.value());
			//generated from below statements...

			//adjacents[xyd.first.x].emplace_back(xyd.first.y, xyd.second);
			//adjacents[xyd.first.y].emplace_back(xyd.first.x, xyd.second);
		}

		uint32_t new_adj = 0;
		for (auto const& a : adjacents) {
			new_adj += a.size();
		}

		//std::cout << "Went from " << old_adj << " to " << new_adj << " by unfolding triangles." << std::endl;

		//for consistency:
		for (auto& a : adjacents) {
			std::sort(a.begin(), a.end());
		}
	}
	qDebug() << "build and unfolding complete!";

	//start creating embedded vertices
	std::vector< std::vector< EmbeddedVertex > > embedded_chains;
	qDebug() << "strating embedded_chains building...";
	// was a big for cycle, but that was because of radius checking, which is not part of my program
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

	//assert(embedded_chains.size() == constraints.size());
	qDebug() << "embedded chains size" << embedded_chains.size() << " constraints size (with empty dummy):" << constraints.size();

	
	EmbeddedPlanarMap< float, SameValue< float >, ReplaceValue< float > > epm;
	qDebug() << "EmbeddedPlanarMap created!";

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

		epm.split_triangles(newVerts, newTris, &split_evs, &split_tris, &epm_to_split);

		std::vector< QVector3D > split_verts;

		qDebug() << "populating split_verts...";
		split_verts.reserve(split_evs.size());
		for (auto const& ev : split_evs) {
			split_verts.emplace_back(ev.interpolate(newVerts));
		}

		qDebug() << "record constrained edges in terms of split_verts";
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
		std::vector< QVector3D > compressed_verts;
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

		qDebug() << "Went from " << newTris.size() << " to (via split) " << split_tris.size() << " to (via discard) " << compressed_tris.size() << " triangles."; //DEBUG

		newMesh.vertices = compressed_verts;
		newMesh.indices = toIntArray(compressed_tris);
		newTriangles = compressed_tris;

		qDebug() << compressed_values.size() << " constrained values and " << compressed_verts.size() << " vertices";


		constrained_values = compressed_values;
}

void KnitGrapher::findFirstActiveChains(std::vector< std::vector< EmbeddedVertex > >* active_chains_,
	std::vector< std::vector< Stitch > >* active_stitches_,
	RowColGraph* graph_) {
	
	assert(active_chains_);
	auto& active_chains = *active_chains_;
	active_chains.clear();

	assert(active_stitches_);
	auto& active_stitches = *active_stitches_;
	active_stitches.clear();

	newTriangles = getTriangles(newMesh);
	//PARANOIA: triangles must reference valid time values:
	for (glm::uvec3 const& tri : newTriangles) {
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
		for (glm::uvec3 const& tri : newTriangles) {
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

		//std::cout << "Considering chain with value range [" << chain_min << ", " << chain_max << "] and neighbor value range [" << adj_min << ", " << adj_max << "]." << std::endl;

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
			QVector3D a = divided_chain[ci - 1].interpolate(newMesh.vertices);
			QVector3D b = divided_chain[ci].interpolate(newMesh.vertices);
			total_length += (b - a).length();
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
	std::cout << "Found " << active_chains.size() << " first active chains." << std::endl;

	if (graph_) {
		for (uint32_t ci = 0; ci < active_chains.size(); ++ci) {
			auto const& chain = active_chains[ci];

			std::vector< float > lengths;
			lengths.reserve(chain.size());
			lengths.emplace_back(0.0f);
			for (uint32_t i = 1; i < chain.size(); ++i) {
				QVector3D a = chain[i - 1].interpolate(newMesh.vertices);
				QVector3D b = chain[i].interpolate(newMesh.vertices);
				lengths.emplace_back(lengths.back() + (b - a).length());
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
QVector3D mix(const QVector3D& wa, const QVector3D& wb, float m) {
	return wa + m * (wb - wa);
}

void KnitGrapher::sampleChain(float spacing,
	std::vector< EmbeddedVertex > const& chain, //in: chain to be sampled
	std::vector< EmbeddedVertex >* sampled_chain_){
	
	auto& sampled_chain = *sampled_chain_;
	sampled_chain.clear();

	//assert(sampled_flags_);
	//auto &sampled_flags = *sampled_flags_;
	//sampled_flags.clear();

	for (uint32_t ci = 0; ci + 1 < chain.size(); ++ci) {
		sampled_chain.emplace_back(chain[ci]);

		QVector3D a = chain[ci].interpolate(newMesh.vertices);
		QVector3D b = chain[ci + 1].interpolate(newMesh.vertices);
		glm::uvec3 common = EmbeddedVertex::common_simplex(chain[ci].simplex, chain[ci + 1].simplex);
		QVector3D wa = chain[ci].weights_on(common);
		QVector3D wb = chain[ci + 1].weights_on(common);

		float length = (b - a).length();
		int32_t insert = std::floor(length / spacing);
		for (int32_t i = 0; i < insert; ++i) {
			float m = float(i + 1) / float(insert + 1);
			sampled_chain.emplace_back(common, mix(wa, wb, m));
		}
	}
	sampled_chain.emplace_back(chain.back());
}



void KnitGrapher::interpolateValues() {
	//assert(constraints.size() == model.vertices.size());
	//assert(values_);
	if (constrained_values.size() != newMesh.vertices.size()) {
		qDebug() << "Constraints size does not match model vertices size" << constrained_values.size() << newMesh.vertices.size();
		return;
	}
	qDebug() << "Interpolating values...";

	qDebug() << newMesh.vertices.size() << " vertices and " << newMesh.indices.size() << " triangles.";
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



	for (const auto& tri : newTriangles) {
		const QVector3D& a = newMesh.vertices[tri.x];
		const QVector3D& b = newMesh.vertices[tri.y];
		const QVector3D& c = newMesh.vertices[tri.z];

		float weight_ab = QVector3D::dotProduct(c - a, b - a) / QVector3D::crossProduct(c - a, b - a).length();
		float weight_bc = QVector3D::dotProduct(a - b, c - b) / QVector3D::crossProduct(a - b, c - b).length();
		float weight_ca = QVector3D::dotProduct(b - c, a - c) / QVector3D::crossProduct(b - c, a - c).length();

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
	coefficients.reserve(newMesh.indices.size() * 3); //more than is needed
	Eigen::VectorXd rhs(total_dofs);

	for (uint32_t i = 0; i < dofs.size(); ++i) {
		if (dofs[i] == -1U) continue;
		//sum adj[x] + one * 1 - c * x = 0.0f
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
	//A = A * A.transpose();
	A.makeCompressed(); //redundant?

	//Eigen::SparseLU< Eigen::SparseMatrix< double > > solver;
	//Eigen::SparseQR< Eigen::SparseMatrix< double >, Eigen::COLAMDOrdering< int > > solver;
	Eigen::SimplicialLDLT< Eigen::SparseMatrix< double > > solver;
	//Eigen::ConjugateGradient< Eigen::SparseMatrix< double > > solver;
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
	//std::cout << solver.iterations() << " interations later..." << std::endl; //DEBUG
	//std::cout << solver.error() << " (estimated error)..." << std::endl; //DEBUG

	values = constrained_values;
	for (uint32_t i = 0; i < dofs.size(); ++i) {
		if (dofs[i] != -1U) values[i] = x[dofs[i]];
	}

	qDebug() << "Interpolation complete!";
	for (float val : values) {
		qDebug() << val;
	}
}

std::vector<GLuint> KnitGrapher::toIntArray(std::vector<glm::uvec3> triangles) {
	std::vector<GLuint> result;

	for (auto t : triangles) {
		result.push_back(t.x);
		result.push_back(t.y);
		result.push_back(t.z);
	}
	return result;
}

bool KnitGrapher::degenerateCheck(std::vector<glm::uvec3> tris) {
	for (auto const& tri : tris) {
		QVector3D const& x = originalMesh.vertices[tri.x];
		QVector3D const& y = originalMesh.vertices[tri.y];
		QVector3D const& z = originalMesh.vertices[tri.z];

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


std::vector<glm::uvec3> KnitGrapher::getTriangles(ObjectMesh const& mesh) {
	std::vector<glm::uvec3> result;

	for (int i = 0; i < mesh.indices.size(); i += 3) {
		glm::uvec3 tri;
		tri.x = mesh.indices[i];
		tri.y = mesh.indices[i + 1];
		tri.z = mesh.indices[i + 2];
		result.push_back(tri);
	}
	return result;

}

void KnitGrapher::extractLevelChains(
	ObjectMesh const& model, //in: model on which to embed vertices
	std::vector< float > const& values, //in: values at vertices
	float const level, //in: level at which to extract chains
	std::vector< std::vector< EmbeddedVertex > >* chains_ //chains of edges at given level
) {
	assert(chains_);
	auto& chains = *chains_;
	chains.clear();

	//embed points along all edges that start below level and end at or above it:
	std::vector<glm::uvec3> modelTriangles = getTriangles(model);
	std::vector< QVector3D > const& verts = model.vertices;
	std::vector< glm::uvec3 > const& tris = modelTriangles;

	std::unordered_map< glm::uvec2, EmbeddedVertex > embedded_pts;
	std::unordered_map< glm::uvec2, QVector3D > pts;
	auto add = [&](uint32_t a, uint32_t b) {
		assert(values[a] < level && values[b] >= level);
		float vmix = (level - values[a]) / (values[b] - values[a]);
		pts[glm::uvec2(a, b)] = mix(verts[a], verts[b], vmix);
		embedded_pts[glm::uvec2(a, b)] = EmbeddedVertex::on_edge(a, b, vmix);
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
		//NOTE: we treat level as "level + epsilon"
		if (values[a] >= level) continue; //all above border
		assert(values[a] < level);
		//NOTE: if values increase along +y, chains should be oriented in the +x direction
		//assuming ccw oriented triangles, this means:

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

void KnitGrapher::trimModel(std::vector< std::vector< EmbeddedVertex > >& left_of,
	std::vector< std::vector< EmbeddedVertex > >& right_of,
	ObjectMesh* clipped_,
	std::vector< EmbeddedVertex >* clipped_vertices_,
	std::vector< std::vector< uint32_t > >* left_of_vertices_, //out (optional): indices of vertices corresponding to left_of chains [may be some rounding]
	std::vector< std::vector< uint32_t > >* right_of_vertices_ //out (optional): indices of vertices corresponding to right_of chains [may be some rounding]
) {

	qDebug() << "trimModel() asserts and paranoias...";

	assert(clipped_);
	auto& clipped = *clipped_;
	clipped.clear();

	assert(clipped_vertices_);
	auto& clipped_vertices = *clipped_vertices_;
	clipped_vertices.clear();


	newTriangles = getTriangles(newMesh);
	{ //PARANOIA: make sure all chains are loops or edge-to-edge:
		std::unordered_set< glm::uvec2 > edges;
		for (auto const& tri : newTriangles) {
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

	qDebug() << "embedding chains using planar map...";

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

	qDebug() << "left_of chains have " << total_chain_edges << " edges.";

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

	qDebug() << "read out left_of_vertices and right_of_vertices from edge information";
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
						//if (s.den > 2) std::cout << s.a << "/" << s.b << " -> [" << ee.first << "-" << ee.second << "] (non-flipped)" << std::endl; //DEBUG
						edge_subedges[glm::uvec2(s.a, s.b)].emplace_back(ee.first, ee.second, s.num, s.den);
					}
					else {
						//if (s.den > 2) std::cout << s.b << "/" << s.a << " -> [" << ee.second << "-" << ee.first << "] (flipped)" << std::endl; //DEBUG
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
				//    a.num / a.den < b.num / b.den
				// => a.num * b.den < b.num * a.den
				return uint64_t(a.num) * uint64_t(b.den) < uint64_t(b.num) * uint64_t(a.den);
				});
			/*if (subedges.size() > 1) {
				//DEBUG:
				for (auto &s : subedges) {
					std::cout << " [" << s.a << "-" << s.b << "]@(" << s.num << "/" << s.den << ")";
				}
				std::cout << std::endl;
			}*/
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
				QVector3D a = epm.vertices[epm_chain[i - 1]].interpolate(newMesh.vertices);
				QVector3D b = epm.vertices[epm_chain[i]].interpolate(newMesh.vertices);
				lengths.emplace_back(lengths.back() + (b - a).length());
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
		} //while (again)

		if (verts_removed) {
			qDebug() << "Removed " << verts_removed << " vertices (that's " << length_removed << " units; " << length_removed / initial_length * 100.0 << "% of the initial length of " << initial_length << " units).";
		}

	};

	qDebug() << "strating loop cleanup";

	for (auto& epm_chain : left_of_epm) {
		cleanup_chain(epm_chain, 1);
	}

	for (auto& epm_chain : right_of_epm) {
		cleanup_chain(epm_chain, -(1 << 8));
	}

	qDebug() << "building split mesh";

	//build split mesh:
	std::vector< EmbeddedVertex > split_verts;
	std::vector< glm::uvec3 > split_tris;
	std::vector< uint32_t > epm_to_split;
	epm.split_triangles(newMesh.vertices, newTriangles, &split_verts, &split_tris, &epm_to_split);


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
				auto f = edge_to_tri.find(glm::uvec2(b, a));
				if (f == edge_to_tri.end()) return;
				int32_t nv = value;
				auto f2 = edge_values.find(glm::uvec2(b, a));
				if (f2 != edge_values.end()) nv += f2->second;

				if (values[f->second] == Unvisited) {
					values[f->second] = nv;
					component.emplace_back(f->second);
				}
				else {
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

	std::vector<glm::uvec3> clippedTriangles = getTriangles(clipped);
	for (auto const& tri : split_tris) {
		if (!keep[&tri - &split_tris[0]]) continue;
		clippedTriangles.emplace_back(
			use_vertex(tri.x),
			use_vertex(tri.y),
			use_vertex(tri.z)
		);
	}
	clipped.vertices.reserve(clipped_vertices.size());
	for (auto const& v : clipped_vertices) {
		clipped.vertices.emplace_back(v.interpolate(newMesh.vertices));
	}

	clipped.indices = toIntArray(clippedTriangles);

	qDebug() << "Trimmed model from " << newTriangles.size() << " triangles on " << newMesh.vertices.size() << " vertices to " << clippedTriangles.size() << " triangles on " << clipped.vertices.size() << " vertices.";
	

	qDebug() << "transform vertex indices for left_of and right_of vertices->clipped model";
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
			//std::cout << "left_of[" << (&chain - &left_of_epm[0]) << "]:"; //DEBUG
			left_of_vertices_->emplace_back(transform_chain(chain));
		}
	}
	if (right_of_vertices_) {
		right_of_vertices_->clear();
		right_of_vertices_->reserve(right_of_epm.size());
		for (auto const& chain : right_of_epm) {
			//std::cout << "right_of[" << (&chain - &right_of_epm[0]) << "]:"; //DEBUG
			right_of_vertices_->emplace_back(transform_chain(chain));
		}
	}
}

void KnitGrapher::peelSlice(std::vector< std::vector< EmbeddedVertex > > & active_chains,
	ObjectMesh* slice_,
	std::vector< EmbeddedVertex >* slice_on_model_,
	std::vector< std::vector< uint32_t > >* slice_active_chains_,
	std::vector< std::vector< uint32_t > >* slice_next_chains_,
	std::vector< bool >* used_boundary_) {

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
		//DEBUG:
		uint32_t loops = 0;
		uint32_t lines = 0;
		for (auto const& chain : active_chains) {
			if (chain[0] == chain.back()) ++loops;
			else ++lines;
		}
		qDebug() << "---- peel slice on [" << loops << " loops and " << lines << " lines] ----";
	}

	ObjectMesh clipped;
	std::vector< EmbeddedVertex > clipped_on_model;
	std::vector< std::vector< EmbeddedVertex > > dummy;

	qDebug() << "trimming model...";
	trimModel(active_chains, dummy, &clipped, &clipped_on_model, nullptr, nullptr);
	qDebug() << "trimming finished!";

	//This version of the code just uses the 3D distance to the curve.
	//might have problems with models that get really close to themselves.

	std::vector< float > values(clipped.vertices.size(), std::numeric_limits< float >::infinity());

	auto do_seg = [&values, &clipped](QVector3D const& a, QVector3D const& b) {
		if (a == b) return;
		QVector3D ab = b - a;
		float limit = QVector3D::dotProduct(ab, ab);
		float inv_limit = 1.0f / limit;
		for (auto const& v : clipped.vertices) {
			float amt = QVector3D::dotProduct(v - a, ab);
			amt = std::max(0.0f, std::min(limit, amt));
			QVector3D pt = (amt * inv_limit) * (b - a) + a;
			float dis2 = (v - pt).lengthSquared();
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

	qDebug() << "values computed!, extracting chains...";
	std::vector<glm::uvec3> clippedTriangles = getTriangles(clipped);
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
			for (auto const& tri : clippedTriangles) {
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
				float amt = ev.weights.y();
				if (next.count(e)) {
					e = glm::uvec2(ev.simplex.y, ev.simplex.x);
					amt = ev.weights.x();
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
				//path.emplace_back(EmbeddedVertex::on_edge(begin_e, begin_ce.along));

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
				QVector3D weights = v.weights.x() * clipped_on_model[v.simplex.x].weights_on(simplex);
				if (v.simplex.y != -1U) weights += v.weights.y() * clipped_on_model[v.simplex.y].weights_on(simplex);
				if (v.simplex.z != -1U) weights += v.weights.z() * clipped_on_model[v.simplex.z].weights_on(simplex);
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
	qDebug() << "starting paranoia check";
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


	qDebug() << "doing second trim (otce nas ktory si na nebesiach...)";
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
void KnitGrapher::stepButtonClicked()
{
	qDebug() << "KnitGrapher received step!";




	if (stepCount == 0) {
		qDebug() << "[Step 0] - peel begin, calling findFirstActiveChains()";
		
		//generate the 3 variables used in findFirstActiveChains()



		findFirstActiveChains(&active_chains, &active_stitches, &graph);

		qDebug() << "emitting the firstActiveChainsCreated signal...";
		emit firstActiveChainsCreated(&active_chains, &active_stitches, &graph);

	}
	if (stepCount == 1) {
		qDebug() << "[Step 1] - slice, calling peelSlice()";
		//peel_slice(parameters, constrained_model, active_chains, &slice, &slice_on_model, &slice_active_chains, &slice_next_chains, &slice_next_used_boundary);
		peelSlice(active_chains, &slice, &sliceOnModel, &sliceActiveChains, &sliceNextChains, &nextUsedBoundary);
		//emit peelSliceDone(&slice, &sliceOnModel, &sliceActiveChains, &sliceNextChains, &nextUsedBoundary);
	}

	stepCount++;
	
	//emit activeChainsFound(data);
}

void KnitGrapher::printConstrainedValues()
{
	for (float val : constrained_values) {
		qDebug() << val;
	}
	// get size of constrained values
	qDebug() << "Constrained Values Size:" << constrained_values.size();
}

void KnitGrapher::constructKnitGraph(std::vector<Constraint*> constraints)
{
	qDebug() << "Constructing Knit Graph, sizes:" << stitchWidth << stitchHeight << modelUnitLength;
	// Step 1: create a new mesh that conforms to the constraints and stitch size given by the user
	qDebug() << "Constraints:" << constraints.size();

	this->constraints = constraints;

	stitchWidth = 3.66f;
	stitchHeight = 1.73f;
	modelUnitLength = 10.0f;
	//this->constraints[2]->timeValue = 0.5f;

	generateTriangles();
	remesh();
	//printConstrainedValues();
	interpolateValues();
	qDebug() << "interpolation done!, emitting...";
	emit knitGraphInterpolated(newMesh, constrained_values);
}
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