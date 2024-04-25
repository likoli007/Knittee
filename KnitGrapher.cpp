#pragma once
#include "KnitGrapher.h"
#include "EmbeddedVertex.h"
#include "EmbeddedPlanarMap.h"
#include <vector>
#include <QSet>
#include <QPair>

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
			if (len_ab2 > maxEdgeLengthSquared) mark(tri.x, tri.y);
			if (len_bc2 > maxEdgeLengthSquared) mark(tri.y, tri.z);
			if (len_ca2 > maxEdgeLengthSquared) mark(tri.z, tri.x);
		}
		if (marked.empty()) {
			break;
		}
		divide(marked, newVerts, paths, newTris);
	}
	qDebug() << "divide() sleu finished! newTris: " << newTris.size();

	//after this while cycle, have newTriangles variable equal to divides newTris...
	//once again the degenerate triangle check in original algortihm
	if (degenerateCheck(newTris)) {
		qDebug() << "Degenerate triangle found after remesh! skipping...";
		return;
	}


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

	//EmbeddedPlanarMap< float, SameValue< float >, ReplaceValue< float > > epm;




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

void KnitGrapher::constructKnitGraph(std::vector<Constraint*> constraints)
{
	qDebug() << "Constructing Knit Graph, sizes:" << stitchWidth << stitchHeight << modelUnitLength;
	// Step 1: create a new mesh that conforms to the constraints and stitch size given by the user
	qDebug() << "Constraints:" << constraints.size();

	this->constraints = constraints;

	generateTriangles();
	remesh();

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