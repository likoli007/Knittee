#pragma once
#define GLM_ENABLE_EXPERIMENTAL	

#include <vector>
#include <QSet>
#include <QPair>
#include <set>
#include <deque>
#include <Eigen/SparseCholesky>
#include <glm/gtx/norm.hpp>
#include <glm/gtx/hash.hpp>
#include <glm/glm.hpp>
#include <QVector3D>
#include "TracedStitch.h"



#include "ObjectMesh.h"
#include "Visualizer.h"
//#include "Stitch.h"
#include "RowColGraph.h"
#include "EmbeddedVertex.h"
#include "EmbeddedPlanarMap.h"
#include "RowColGraph.h"
#include "Link.h"
#include "EmbeddedVertex.h"
#include "EmbeddedPlanarMap.h"

#include "HelperText.h"





/*
* Helper class used in the autoknit algorithms, does not require its own .h file...
*/
struct OnChainStitch {
	enum On : uint8_t { OnNone, OnActive, OnNext } on;
	uint32_t chain;
	uint32_t stitch;
	enum Type : uint8_t { TypeBegin, TypeStitch, TypeEnd } type = TypeStitch;
	OnChainStitch(On on_ = OnNone, uint32_t chain_ = -1U, uint32_t stitch_ = -1U) : on(on_), chain(chain_), stitch(stitch_) { }

	bool operator<(OnChainStitch const& o) const {
		if (on != o.on) return on < o.on;
		else if (chain != o.chain) return chain < o.chain;
		else return stitch < o.stitch;
	}
	bool operator==(OnChainStitch const& o) const {
		return on == o.on && chain == o.chain && stitch == o.stitch && type == o.type;
	}
	bool operator!=(OnChainStitch const& o) const {
		return !(*this == o);
	}
};

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
	/////////////////////////////////parameter values/////////////////////////////////////////////
	float stitchWidth;
	float stitchHeight;		//both float values in mm
	float modelUnitLength;	//the length of a unit in the model in mm
	int stepCount;

	///////////class member variables and functions used mainly in remesh and interpolation///////
	AutoKnitMesh originalMesh, newMesh;
	std::vector<Constraint*> constraints;
	std::vector<float> constrained_values;

	void quad(std::vector<glm::uvec3>& new_tris, std::vector<glm::vec3> const& verts, uint32_t a, uint32_t b, uint32_t c, uint32_t d);
	void divide(std::unordered_set< glm::uvec2 >const& marked,
		std::vector< glm::vec3 >& verts,
		std::vector<std::vector<uint32_t>>& paths,
		std::vector<glm::uvec3>& tris);
	void remesh();
	float getMaxEdgeLength() const;
	bool degenerateCheck(std::vector<glm::uvec3> tris);

	//so. much. passing. by. reference...
	void unfold(uint32_t depth, uint32_t root, glm::vec2 const& flat_root,
		uint32_t ai, glm::vec2 const& flat_a,
		uint32_t bi, glm::vec2 const& flat_b,
		glm::vec2 const& limit_a, glm::vec2 const& limit_b, std::unordered_map< glm::uvec2, uint32_t >  const& opposite,
		std::vector<glm::vec3> const& newVertices, std::unordered_map< glm::uvec2, float >& min_dis);

	GLuint lookup(uint32_t a, uint32_t b, std::unordered_map< glm::uvec2, uint32_t >& marked_verts) const;
	bool isCcw(glm::vec2 const& a, glm::vec2 const& b, glm::vec2 const& c);
	float& get_dis(uint32_t a, uint32_t b, std::unordered_map< glm::uvec2, float >& min_dis);
	void interpolateValues();

	/////////////class member variables and functions used mainly in find step////////////////////
	void findFirstActiveChains(std::vector< std::vector< EmbeddedVertex > >* active_chains_,
		std::vector< std::vector< Stitch > >* active_stitches_,
		RowColGraph* graph_);

	///////////////class member variables and functions used mainly in peeling////////////////////
	std::vector< std::vector< EmbeddedVertex > > active_chains;
	std::vector< std::vector< Stitch > > active_stitches;
	RowColGraph graph;
	AutoKnitMesh slice;
	std::vector< EmbeddedVertex > sliceOnModel;
	std::vector< std::vector< uint32_t > > sliceActiveChains;
	std::vector< std::vector< uint32_t > > sliceNextChains;
	std::vector< bool > nextUsedBoundary;
	std::vector <float> sliceTimes;
	std::vector< std::vector< EmbeddedVertex > > nextActiveChains;
	std::vector< std::vector< Stitch > > nextActiveStitches;

	void extractLevelChains(
		AutoKnitMesh const& model, //in: model on which to embed vertices
		std::vector< float > const& values, //in: values at vertices
		float const level, //in: level at which to extract chains
		std::vector< std::vector< EmbeddedVertex > >* chains_ //chains of edges at given level
	);

	void trimModel(std::vector< std::vector< EmbeddedVertex > >& left_of,
		std::vector< std::vector< EmbeddedVertex > >& right_of,
		AutoKnitMesh* clipped_,
		std::vector< EmbeddedVertex >* clipped_vertices_,
		std::vector< std::vector< uint32_t > >* left_of_vertices_, //out (optional): indices of vertices corresponding to left_of chains [may be some rounding]
		std::vector< std::vector< uint32_t > >* right_of_vertices_ //out (optional): indices of vertices corresponding to right_of chains [may be some rounding]
	);
	void sampleChain(
		float spacing,
		std::vector< EmbeddedVertex > const& chain,
		std::vector< EmbeddedVertex >* sampled_chain_//in: chain to be sampled
	);
	void peelSlice(std::vector< std::vector< EmbeddedVertex > >& active_chains,
		AutoKnitMesh* slice_,
		std::vector< EmbeddedVertex >* slice_on_model_,
		std::vector< std::vector< uint32_t > >* slice_active_chains_,
		std::vector< std::vector< uint32_t > >* slice_next_chains_,
		std::vector< bool >* used_boundary_);

	///////////////class member variables and functions used mainly in linking/////////////////////////////////////
	std::vector<std::vector< Stitch>> nextStitches;
	std::vector<Link> links;

	bool fillUnassigned(std::vector< uint32_t >& closest, std::vector< float > const& weights, bool is_loop);
	void flatten(std::vector< uint32_t >& closest, std::vector< float > const& weights, bool is_loop);
	void optimalLink(
		float target_distance, bool do_roll,
		std::vector< glm::vec3 > const& source,
		std::vector< bool > const& source_linkone,
		std::vector< glm::vec3 > const& target,
		std::vector< bool > const& target_linkone,
		std::vector< std::pair< uint32_t, uint32_t > >* links_);
	void linkChains(
		AutoKnitMesh const& slice, //in: slice on which the chains reside
		std::vector< float > const& slice_times, //in: time field (times @ vertices), for slice
		std::vector< std::vector< uint32_t > > const& active_chains, //in: current active chains (slice vertex #'s)
		std::vector< std::vector< Stitch > > const& active_stitches, //in: current active stitches, sorted by time
		std::vector< std::vector< uint32_t > > const& next_chains, //in: current next chains (slice vertex #'s)
		std::vector< bool > const& next_used_boundary, //in: did next chain use boundary? (forces no discard)
		//need this or slice_times (above) std::vector< std::vector< bool > > const &discard_segments,
		std::vector< std::vector< Stitch > >* next_stitches, //out: next active stitches
		std::vector< Link >* links //out: active_chains[from_chain][from_vertex] -> linked_next_chains[to_chain][to_vertex] links
	);

	///////////////////////////class member variables and functions used mainly in building////////////////
	void embeddedPathSimple(
		AutoKnitMesh const& model,
		EmbeddedVertex const& source,
		EmbeddedVertex const& target,
		std::vector< EmbeddedVertex >* path_ //out: path; path[0] will be source and path.back() will be target
	);

	void embeddedPath(
		AutoKnitMesh const& model,
		EmbeddedVertex const& source,
		EmbeddedVertex const& target,
		std::vector< EmbeddedVertex >* path_ //out: path; path[0] will be source and path.back() will be target
	);


	void buildNextActiveChains(
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
	);
	///////////////////////class variables and functions used for graph tracing/////////////////////////////
	std::vector< TracedStitch > tracedMesh;
	void traceGraph(
		RowColGraph const& graph, //in: row-column graph
		std::vector< TracedStitch >* traced_, //out:traced list of stitches
		AutoKnitMesh* DEBUG_model_ //in (optional): model
	);


	////////////////////////////////helper functions///////////////////////////////////////////////////////
	void clearPeeling();

	float getChainSampleSpacing() const {
		return 0.25f * stitchWidth / modelUnitLength;
	
	}

	float getMaxPathSampleSpacing() const {
		return 0.02f * std::min(stitchWidth, 2.0f * stitchHeight) / modelUnitLength;
	}

	std::vector<glm::uvec3> getTriangles(ObjectMesh const& model);
	std::vector<GLuint> toIntArray(std::vector<glm::uvec3>);
	



//signals and slots, public functions, all 'interface' with the rest of the program	
public:
	KnitGrapher(QObject* parent = nullptr);
	void constructNewMesh(std::vector<Constraint*> constraints);
	void reset();
public slots:
	void setStitchWidth(float width);
	void setStitchHeight(float height);
	void setModelUnitLength(float length);
	void setOriginalMesh(ObjectMesh mesh);
	void stepButtonClicked();
	void traceButtonClicked();
signals:
	void knitGraphInterpolated(ObjectMesh mesh, std::vector<float> values);
	void firstActiveChainsCreated(std::vector< std::vector< EmbeddedVertex > >* active_chains,
		std::vector< std::vector< Stitch > >* active_stitches,
		RowColGraph* graph);
	void peelSliceDone(ObjectMesh* slice_, std::vector< std::vector< uint32_t > >* slice_active_chains_, std::vector< std::vector< uint32_t > >* slice_next_chains_);
	void linkChainsDone(std::vector< std::vector< Stitch > >*, std::vector< Link >*);
	void nextActiveChainsDone(std::vector< std::vector< EmbeddedVertex > >*);
	void knitGraphCreated();
	void knitGraphTraced(std::vector< TracedStitch >* traced_);
	void helpBoxCommunication(QString message);
	void helpBoxAppend(QString message);
};

