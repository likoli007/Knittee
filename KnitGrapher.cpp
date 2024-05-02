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


bool KnitGrapher::fillUnassigned(std::vector< uint32_t >& closest, std::vector< float > const& weights, bool is_loop) {
	qDebug() << "fillUnassigned() start";
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
	qDebug() << "fillUnassigned() finish!";
	return true;
}

void KnitGrapher::optimalLink(
	float target_distance, bool do_roll,
	std::vector< QVector3D > const& source,
	std::vector< bool > const& source_linkone,
	std::vector< QVector3D > const& target,
	std::vector< bool > const& target_linkone,
	std::vector< std::pair< uint32_t, uint32_t > >* links_) {

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
		float dis = (source[si] - target[ti]).length() - target_distance;
		return dis * dis;
	};

	//TODO: could refine this a fair bit & build a table with real link counts
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

void KnitGrapher::flatten(std::vector< uint32_t >& closest, std::vector< float > const& weights, bool is_loop) {
	assert(closest.size() == weights.size());
	if (closest.empty()) return;
	qDebug() << "flatten start!";
	//make sure that 'closest' looks like:
	//  a a a b b b c c c
	//   a a b b b b c c
	// that is, can be flattened to the knitting machine bed while preserving constituent cycles
	// One view of this: if you start at some stitch A, then the left side should be a mirror of the right side
	//  (with the possible exception that some symbols may be skipped on the left or right)

	//(a) condense closest into short list of symbols:
	qDebug() << "condensing...";
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
	//DEBUG: don't do these checks; exercise the code a bit more instead:
	//if (symbols.size() == 2) return; //two symbols without alternation
	//if (is_loop && symbols.size() == 3 && symbols[0] == symbols.back()) return; //two symbols without alternation (loop version)

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
		//std::cout << state.to_string(bits) << " from " << from.to_string(bits) << " cost " << cost << std::endl; //DEBUG

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

	qDebug() << "queueing starting states...";
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

	qDebug() << "reading back states...";
	std::vector< State::Packed > path;
	path.emplace_back(finished.state);
	path.emplace_back(finished.from);
	while (true) {
		auto f = visited.find(path.back());
		if (f == visited.end()) break;
		path.emplace_back(f->second.second);
	}
	std::reverse(path.begin(), path.end());

	qDebug() << "Starting keep loop...";
	std::vector< int8_t > keep(bit_symbols.size(), -1);
	for (uint32_t i = 1; i < path.size(); ++i) {
		State state = State::unpack(path[i - 1]);
		State next = State::unpack(path[i]);

		//std::cout << state.to_string(bits) << " -> " << next.to_string(bits) << ": "; std::cout.flush(); //DEBUG
		if (state.min != next.min && state.max != next.max) {
			//a(bc)d -> (abcd), keep 'a' (next.min), 'd' (state.max)
			assert(bit_symbols[next.min].first == bit_symbols[state.max].first);
			assert(next.current == bit_symbols[state.min].first);
			//std::cout << "keep " << int32_t(next.min) << " (\"" << symbols[next.min].first << "\")" << ", " << int32_t(state.max) << " (\"" << symbols[state.max].first << "\")" << std::endl; //DEBUG
			assert(keep[next.min] == -1);
			assert(keep[state.max] == -1);
			keep[next.min] = keep[state.max] = 1;
		}
		else if (state.min != next.min) {
			//a(bc)d -> (abc)d, keep/discard next.min
			if (bit_symbols[next.min].first == next.current) {
				//std::cout << "keep " << int32_t(next.min) << " (\"" << symbols[next.min].first << "\")" << std::endl; //DEBUG
				assert(keep[next.min] == -1);
				keep[next.min] = 1;
			}
			else {
				//std::cout << "discard " << int32_t(next.min) << " (\"" << symbols[next.min].first << "\")" << std::endl; //DEBUG
				assert(keep[next.min] == -1);
				keep[next.min] = 0;
			}
		}
		else {
			assert(state.max != next.max);
			//a(bc)d -> a(bcd), keep/discard state.max
			if (bit_symbols[state.max].first == next.current) {
				//std::cout << "keep " << int32_t(state.max) << " (\"" << symbols[state.max].first << "\")" << std::endl; //DEBUG
				assert(keep[state.max] == -1);
				keep[state.max] = 1;
			}
			else {
				//std::cout << "discard " << int32_t(state.max) << " (\"" << symbols[state.max].first << "\")" << std::endl; //DEBUG
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

	qDebug() << "doing relabelling";
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
			//if (relabel.back()) {
			//	std::cout << ' ' << 'x' << int32_t(c) << 'x'; //DEBUG
			//} else {
			//	std::cout << ' ' << ' ' << int32_t(c) << ' '; //DEBUG
			//}
		}
		//std::cout << std::endl; //DEBUG
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
			//std::cout << "Relabelling [" << first << ", " << last << "] using " << before << "/" << after << std::endl; //DEBUG

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
void KnitGrapher::linkChains(
	ObjectMesh const& slice, //in: slice on which the chains reside
	std::vector< float > const& slice_times, //in: time field (times @ vertices), for slice
	std::vector< std::vector< uint32_t > > const& active_chains, //in: current active chains (slice vertex #'s)
	std::vector< std::vector< Stitch > > const& active_stitches, //in: current active stitches, sorted by time
	std::vector< std::vector< uint32_t > > const& next_chains, //in: current next chains (slice vertex #'s)
	std::vector< bool > const& next_used_boundary, //in: did next chain use boundary? (forces no discard)
	//need this or slice_times (above) std::vector< std::vector< bool > > const &discard_segments,
	std::vector< std::vector< Stitch > >* next_stitches_, //out: next active stitches
	std::vector< Link >* links_ //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
) {

	qDebug() << "linkChains() starting!";
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

	qDebug() << "argument preparation phase complete!";

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
				QVector3D const& a = slice.vertices[chain[vi - 1]];
				QVector3D const& b = slice.vertices[chain[vi]];
				total_length += (b - a).length();
				lengths.emplace_back(total_length);
			}
			assert(lengths.size() == chain.size());
		}
		return all_lengths;
	};
	std::vector< std::vector< float > > active_lengths = make_lengths(active_chains);
	std::vector< std::vector< float > > next_lengths = make_lengths(next_chains);

	qDebug() << "make_lengths() within linkChains() complete!";

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
	qDebug() << "stitch removal complete!";

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

	qDebug() << "finding segments of active and next chains that are mutual nearest neighbors";

	std::vector< std::vector< uint32_t > > active_closest;
	std::vector< std::vector< uint32_t > > next_closest;

	std::vector<glm::uvec3> sliceTriangles = getTriangles(slice);

	{ //find closest pairs:

		std::vector< std::vector< uint32_t > > adj(slice.vertices.size());
		{ //adjacency matrix -- always handy:
			std::unordered_set< glm::uvec2 > edges;
			for (auto const& tri : sliceTriangles) {
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
						float nd = d + (slice.vertices[n] - slice.vertices[at]).length();
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
#if 0
		//DEBUG:
		for (auto const& ac : active_closest) {
			std::cout << "active_closest[" << (&ac - &active_closest[0]) << "]:";
			for (auto c : ac) {
				std::cout << " " << int32_t(c);
			}
			std::cout << "\n";
		}
		std::cout.flush();
		for (auto const& nc : next_closest) {
			std::cout << "next_closest[" << (&nc - &next_closest[0]) << "]:";
			for (auto c : nc) {
				std::cout << " " << int32_t(c);
			}
			std::cout << "\n";
		}
		std::cout.flush();
#endif
	}


	qDebug() << "discarding non-mutuals...";
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

		if (discarded) qDebug() << "Discarded " << discarded << " non-mutual segment matches.";

		return discarded > 0;
	};

	//start by removing any non-mutual links:
	discard_nonmutual();

	qDebug() << "filling discarded areas with adjacent information";
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
	qDebug() << "[CHECKPOINT] linkChains() fill-flatten loop complete!";

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

	qDebug() << "building parametric segments into matches";
	//build parametric segments into matches:
	for (auto& closest : active_closest) {
		uint32_t ai = &closest - &active_closest[0];
		auto const& lengths = active_lengths[ai];
		assert(lengths.size() == closest.size() + 1);

		for (uint32_t begin = 0; begin < closest.size(); /* later */) {
			uint32_t end = begin + 1;
			while (end < closest.size() && closest[end] == closest[begin]) ++end;
			assert(end < lengths.size());
			//std::cout << "[" << begin << ", " << end << ") becomes "; //DEBUG
			matches[std::make_pair(ai, closest[begin])].active.emplace_back(lengths[begin] / lengths.back(), lengths[end] / lengths.back(), closest[begin]);
			auto m = matches[std::make_pair(ai, closest[begin])].active;
			//std::cout << "[" << m.back().begin << ", " << m.back().end << ')'; //DEBUG
			//std::cout << " = [" << lengths[begin] << ", " << lengths[end] << ") / " << lengths.back() << std::endl; //DEBUG
			begin = end;
		}
	}

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

	qDebug() << "doing stitch assignments";

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
				//std::cout << "match[" << int32_t(anm.first.first) << "," << int32_t(anm.first.second) << "] gives segment [" << bse.begin << ',' << bse.end << ')' << std::endl; //DEBUG
			}
		}

		//assign stitches to segments:
		for (uint32_t ai = 0; ai < active_chains.size(); ++ai) {
			auto& segments = active_segments[ai];
			assert(!segments.empty());

			std::sort(segments.begin(), segments.end(), [](BeginEndStitches const* a, BeginEndStitches const* b) {
				return a->begin < b->begin;
				});

			//DEBUG:
			//std::cout << "active[" << ai << "] segments:";
			//for (auto seg : segments) {
			//	std::cout << ' ' << '[' << seg->begin << ',' << seg->end << ')';
			//}
			//std::cout << std::endl;

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

	qDebug() << "merging segments";
	//merge segments because it's probably more convenient for balancing:
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
	if (active_merges) qDebug() << "Merged " << active_merges << " active segments." ;

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
	if (next_merges) qDebug() << "Merged " << next_merges << " next segments." ;

	qDebug() << "balancing stitch assignments";
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
				qDebug() << "Balanced a merge: old: " << old_back << " " << old_front << " new: " << new_back << " " << new_front;
				//} else {
				//if (old_back == new_back && old_front == new_front) {
					//std::cout << "NOTE: merge was already balanced." << std::endl;
			}

		}
	}

	qDebug() << "looking for merges/splits";
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
	qDebug() << "allocating next stitches...";
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
		//std::cout << "Matching to"; //DEBUG
		for (auto const& be : match.next) {
			//	std::cout << " [" << be.begin << ", " << be.end << ")"; std::cout.flush(); //DEBUG
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
					match.next[a.bi].stitches += 1; //.emplace_back(next_stitches[anm.first.second].size() + new_stitches.size() - 1); //track stitch index
				}
			}
		}
		assert(new_stitches.size() == stitches);

		next_stitches[anm.first.second].insert(next_stitches[anm.first.second].end(), new_stitches.begin(), new_stitches.end());
	} //end stitch allocation
	qDebug() << "balancingnew stitch allocation";

	{ //balance new stitch allocations for merges:
		//sort ranges from matches back to nexts:
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
				qDebug() << "Balanced a split:\n" << "old: " << old_back << "\n" << "     " << old_front << "\n" << "new: " << new_back << "\n" << "     " << new_front << "\n";


				//} else {
				//if (old_back == new_back && old_front == new_front) {
				//	std::cout << "NOTE: split was already balanced." << std::endl;
			}


		}
	}

	qDebug() << "making stitch info...";
	for (auto& stitches : next_stitches) {
		std::stable_sort(stitches.begin(), stitches.end(), [](Stitch const& a, Stitch const& b) {
			return a.t < b.t;
			});
	}

	auto make_stitch_info = [&slice](
		std::vector< uint32_t > const& chain,
		std::vector< float > const& lengths,
		std::vector< Stitch > const& stitches,
		std::vector< QVector3D >* stitch_locations_,
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
					mix(
						slice.vertices[chain[i - 1]],
						slice.vertices[chain[i]],
						(l - *(li - 1)) / (*li - *(li - 1))
					)
				);
				stitch_linkones.emplace_back(si->flag == Stitch::FlagLinkOne);
			}
	};

	std::vector< std::vector< QVector3D > > all_active_stitch_locations(active_chains.size());
	std::vector< std::vector< bool > > all_active_stitch_linkones(active_chains.size());

	std::vector< std::vector< QVector3D > > all_next_stitch_locations(next_chains.size());
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

	qDebug() << "starting to connect stitches...";
	for (auto const& anm : matches) {
		Match const& match = anm.second;

		std::vector< uint32_t > next_stitch_indices;
		std::vector< QVector3D > next_stitch_locations;
		std::vector< bool > next_stitch_linkones;
		std::unordered_set< uint32_t >& next_claimed = all_next_claimed[anm.first.second];
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


		std::vector< uint32_t > active_stitch_indices;
		std::vector< QVector3D > active_stitch_locations;
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

		qDebug() << "building links...";
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
#ifdef USE_EVEN
			auto try_links_even = [&](uint32_t roll_active, uint32_t roll_new) {
				std::vector< std::pair< uint32_t, uint32_t > > possible_links;

				if (active_stitch_locations.size() <= next_stitch_locations.size()) {
					//evenly distribute increases among the non-linkone stitches:
					uint32_t total = 0;
					for (auto l : active_stitch_linkones) {
						if (!l) ++total;
					}
					uint32_t increases = next_stitch_locations.size() - active_stitch_locations.size();
					std::vector< bool > inc(total, false);
					for (uint32_t i = 0; i < increases; ++i) {
						assert(inc[i * total / increases] == false);
						inc[i * total / increases] = true;
					}
					uint32_t n = 0;
					uint32_t i = 0;
					for (uint32_t a = 0; a < active_stitch_locations.size(); ++a) {
						uint32_t ra = (a + roll_active) % active_stitch_locations.size();
						uint32_t rn = (n + roll_new) % next_stitch_locations.size();
						possible_links.emplace_back(ra, rn);
						++n;
						if (!active_stitch_linkones[ra]) {
							assert(i < inc.size());
							if (inc[i]) {
								rn = (n + roll_new) % next_stitch_locations.size();
								possible_links.emplace_back(ra, rn);
								++n;
							}
							++i;
						}
					}
					assert(i == total);
					assert(n == next_stitch_locations.size());
				}
				else if (active_stitch_locations.size() > next_stitch_locations.size()) {
					//evenly distribute decreases among the non-linkone stitches:
					uint32_t total = 0;
					for (auto l : next_stitch_linkones) {
						if (!l) ++total;
					}
					assert(total > 0);
					uint32_t decreases = active_stitch_locations.size() - next_stitch_locations.size();
					std::vector< bool > dec(total, false);
					for (uint32_t i = 0; i < decreases; ++i) {
						assert(dec[i * total / decreases] == false);
						dec[i * total / decreases] = true;
					}
					uint32_t a = 0;
					uint32_t i = 0;
					for (uint32_t n = 0; n < next_stitch_locations.size(); ++n) {
						uint32_t ra = (a + roll_active) % active_stitch_locations.size();
						uint32_t rn = (n + roll_new) % next_stitch_locations.size();
						possible_links.emplace_back(ra, rn);
						++a;
						if (!next_stitch_linkones[rn]) {
							assert(i < dec.size());
							if (dec[i]) {
								ra = (a + roll_active) % active_stitch_locations.size();
								possible_links.emplace_back(ra, rn);
								++a;
							}
							++i;
						}
					}
					assert(i == total);
					assert(a == active_stitch_locations.size());
				}
				float const row_height = 2.0f * parameters.stitch_height_mm / parameters.model_units_mm;
				float cost = 0.0f;
				for (auto const& p : possible_links) {
					float len = glm::length(next_stitch_locations[p.second] - active_stitch_locations[p.first]);
					cost += (len - row_height) * (len - row_height);
				}

				if (cost < best_cost) {
					best_cost = cost;
					best_links = possible_links;
				}
			};


			if (active_stitch_locations.size() >= next_stitch_locations.size()) {
				for (uint32_t roll_active = 0; roll_active < active_stitch_locations.size(); ++roll_active) {
					try_links_even(roll_active, 0);
				}
			}
			else {
				for (uint32_t roll_new = 0; roll_new < next_stitch_locations.size(); ++roll_new) {
					try_links_even(0, roll_new);
				}
			}
#endif //USE_EVEN

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
	qDebug() << "links built!";

	qDebug() << "doing final paranoia...";
	//PARANOIA: every stitch should have been claimed
	for (uint32_t ai = 0; ai < active_chains.size(); ++ai) {
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
	qDebug() << "linkChains() finished!";
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
		
		sliceTimes.clear();
		sliceTimes.reserve(sliceOnModel.size());
		for (auto& ev : sliceOnModel) {
			sliceTimes.emplace_back(ev.interpolate(constrained_values));
		}
		
		//emit peelSliceDone(&slice, &sliceOnModel, &sliceActiveChains, &sliceNextChains, &nextUsedBoundary);
	}
	if (stepCount == 2) {
		qDebug() << "[Step 2] - link, calling linkChains()";
		linkChains(slice, sliceTimes, sliceActiveChains, active_stitches, sliceNextChains, nextUsedBoundary, &nextStitches, &links);
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