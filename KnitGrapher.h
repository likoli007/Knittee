#pragma once

#include "ObjectMesh.h"
#include "Visualizer.h"
//#include "myQVectors.h"

#include <glm/glm.hpp>
#include <QVector3D>
#include "Stitch.h"
#include "RowColGraph.h"
#include "EmbeddedVertex.h"
//#include "Constraint.h"
#include "EmbeddedPlanarMap.h"
#include "RowColGraph.h"



/*
* This class will be used to construct the KnitGraph from the parameters given by the user
* A lot of this code will be implemented by refering to the algorithms described by AUTOKNIT
*/
class KnitGrapher : public QObject
{
	Q_OBJECT

//for now defining helper classes inside KnitGrapher, may be moved in the future
//EPM value that can track edge splits:
struct Edge {
		uint32_t a;
		uint32_t b;
		enum Type : uint8_t {
			Initial, //a,b are vertex indices
			Reverse, //a,b are the same edge index
			Combine, //a,b are edge indices
			SplitFirst, //a,b are same edge index
			SplitSecond, //a,b are edge edge index
		} type;
		Edge(Type type_, uint32_t a_, uint32_t b_) : a(a_), b(b_), type(type_) { }
	};

	static std::vector< Edge > globedges;

	struct Value {
		int32_t sum;
		uint32_t edge;
		struct Reverse {
			static inline void reverse(Value* v) {
				v->sum = -v->sum;
				globedges.emplace_back(Edge::Reverse, v->edge, v->edge);
				v->edge = globedges.size() - 1;
			}
		};
		struct Combine {
			static inline void combine(Value* v, Value const& b) {
				v->sum += b.sum;
				globedges.emplace_back(Edge::Combine, v->edge, b.edge);
				v->edge = globedges.size() - 1;
			}
		};
		struct Split {
			static inline void split(Value const& v, Value* first, Value* second) {
				first->sum = v.sum;
				second->sum = v.sum;

				globedges.emplace_back(Edge::SplitFirst, v.edge, v.edge);
				first->edge = globedges.size() - 1;
				globedges.emplace_back(Edge::SplitSecond, v.edge, v.edge);
				second->edge = globedges.size() - 1;
			}
		};
	};














private:
	ObjectMesh originalMesh, newMesh;
	float stitchWidth;
	float stitchHeight;		//both float values in mm
	float modelUnitLength;	//the length of a unit in the model in mm

	const float minEdgeRatio = 0.3f; //'smallest allowed smallest-to-largest edge ratio in a triangle' 
	const float minEdgeRatioSquared = minEdgeRatio * minEdgeRatio;
	//maybe let them be set by the user?

	std::vector<Constraint*> constraints;
	std::vector<float> constrained_values;

	std::vector<glm::uvec3> oldTriangles;
	std::vector<glm::uvec3> newTriangles;


	//////////////////////class member variables used in peeling///////////////////////////////
	std::vector< std::vector< EmbeddedVertex > > active_chains;
	std::vector< std::vector< Stitch > > active_stitches;
	RowColGraph graph;

	ObjectMesh slice;
	std::vector< EmbeddedVertex > sliceOnModel;
	std::vector< std::vector< uint32_t > > sliceActiveChains;
	std::vector< std::vector< uint32_t > > sliceNextChains;
	std::vector< bool > nextUsedBoundary;



	float getChainSampleSpacing() const {
		return 0.25f * stitchWidth / modelUnitLength;
	
	}
	std::vector<glm::uvec3> getTriangles(ObjectMesh const& model);

	void extractLevelChains(
		ObjectMesh const& model, //in: model on which to embed vertices
		std::vector< float > const& values, //in: values at vertices
		float const level, //in: level at which to extract chains
		std::vector< std::vector< EmbeddedVertex > >* chains_ //chains of edges at given level
	);

	//trim the constrained model according to the parameters
	void trimModel(std::vector< std::vector< EmbeddedVertex > > & left_of,
		std::vector< std::vector< EmbeddedVertex > > & right_of, 
		ObjectMesh* clipped_,
		std::vector< EmbeddedVertex >* clipped_vertices_,
		std::vector< std::vector< uint32_t > >* left_of_vertices_, //out (optional): indices of vertices corresponding to left_of chains [may be some rounding]
		std::vector< std::vector< uint32_t > >* right_of_vertices_ //out (optional): indices of vertices corresponding to right_of chains [may be some rounding]
	);



	void sampleChain(
		float spacing,
		std::vector< EmbeddedVertex > const& chain,
		std::vector< EmbeddedVertex >* sampled_chain_//in: chain to be sampled
		);


	void peelSlice(std::vector< std::vector< EmbeddedVertex > > const& active_chains,
		ObjectMesh* slice_,
		std::vector< EmbeddedVertex >* slice_on_model_,
		std::vector< std::vector< uint32_t > >* slice_active_chains_,
		std::vector< std::vector< uint32_t > >* slice_next_chains_,
		std::vector< bool >* used_boundary_);

	int stepCount = 0;

	void quad(std::vector<glm::uvec3>& new_tris, std::vector<QVector3D>& verts, GLuint a, GLuint b, GLuint c, GLuint d);
	void divide(QSet<QPoint>& marked,
		std::vector<QVector3D>& verts,
		std::vector<std::vector<GLuint>>& paths,
		std::vector<glm::uvec3>& tris);
	void remesh();
	float getMaxEdgeLength();
	bool degenerateCheck(std::vector<glm::uvec3> tris);

	std::vector<GLuint> toIntArray(std::vector<glm::uvec3>);

	void generateTriangles();
	void interpolateValues();
	void printConstrainedValues();

	void findFirstActiveChains(std::vector< std::vector< EmbeddedVertex > >* active_chains_,
		std::vector< std::vector< Stitch > >* active_stitches_,
		RowColGraph* graph_);

	//so. much. passing. by. reference. need to make some variables class members.... TODO!!!
	void unfold(GLuint depth, GLuint root, QVector2D const& flat_root,
		GLuint ai, QVector2D const& flat_a,
		GLuint bi, QVector2D const& flat_b,
		QVector2D const& limit_a, QVector2D const& limit_b, QHash< QPoint, GLuint > const& opposite,
		std::vector<QVector3D> const& newVertices, QHash< QPoint, float >& min_dis);

	GLuint lookup(GLuint a, GLuint b, QHash<QPoint, GLuint>& marked_verts);
public:
	KnitGrapher(QObject* parent = nullptr);
	void constructKnitGraph(std::vector<Constraint*> constraints);
public slots:
	void setStitchWidth(float width);
	void setStitchHeight(float height);
	void setModelUnitLength(float length);
	void setOriginalMesh(ObjectMesh mesh);
	void stepButtonClicked();
signals:
	void knitGraphInterpolated(ObjectMesh mesh, std::vector<float> values);
	void firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
		std::vector< std::vector< Stitch > >* active_stitches,
		RowColGraph* graph);
};

