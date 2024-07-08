#ifndef TYPES_H
#define TYPES_H

// CGAL
// #ifndef GLOG_NO_ABBREVIATED_SEVERITIES
// #define GLOG_NO_ABBREVIATED_SEVERITIES
// #endif
// #define CGAL_PMP_USE_CERES_SOLVER
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Triangulation_data_structure_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_segment_primitive.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Surface_mesh_shortest_path.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
// Boost
#include <boost/property_map/property_map.hpp>

// Eigen
#include <Eigen/Core>
#include <Eigen/SparseCore>
//#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseCholesky>
#include <Eigen/Geometry> 
// std
#include <map>
#include <set>
#include <queue>
#include <vector>
#include <unordered_set>
#include <string>

// OpenMP
//#include <omp.h>

namespace OrthoLab {

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

// Simple geometric types
typedef Kernel::FT                  FT;
typedef Kernel::Point_2             Point_2;
typedef Kernel::Point_3             Point;
typedef Kernel::Segment_2           Segment_2;
typedef Kernel::Vector_3            Vector;
typedef Kernel::Triangle_3          Triangle;
typedef Kernel::Segment_3           Segment;
typedef CGAL::Bbox_3                Bbox;
typedef Kernel::Ray_3               Ray;
typedef Kernel::Line_3              Line;
typedef Kernel::Plane_3             Plane;
typedef Kernel::Line_2              Line_2;
typedef Kernel::Sphere_3            Sphere;
typedef Kernel::Aff_transformation_3  Transformation;

// Polyhedron
typedef CGAL::Polyhedron_3<Kernel>                              Polyhedron;
typedef Polyhedron::Vertex_handle                               Vertex_handle;
typedef Polyhedron::Edge_iterator                               Edge_iterator;
typedef Polyhedron::Halfedge_handle                             Halfedge_handle;
typedef Polyhedron::Facet_handle                                Facet_handle;
typedef Polyhedron::Facet_iterator                              Facet_iterator;
typedef Polyhedron::Halfedge_iterator                           Halfedge_iterator;
typedef Polyhedron::Halfedge_around_facet_circulator            Halfedge_facet_circulator;
typedef Polyhedron::Halfedge_around_vertex_circulator           Halfedge_vertex_circulator;

// Polyhedron related structures
typedef std::map<Facet_handle, bool>                     FaceBoolMap;
typedef std::map<Vertex_handle, bool>                    VertBoolMap;
typedef std::map<Vertex_handle, int>                     VertIntMap;
typedef std::map<Halfedge_iterator, int>                 HEdgeIntMap;
typedef std::vector<Halfedge_handle>                     HEdgeList;
typedef std::set<Halfedge_handle>                        HEdgeSet;
typedef std::set<Facet_handle>                           FacetSet;
typedef std::set<Vertex_handle>                          VertexSet;
typedef std::vector<Facet_handle>                        FacetList;
typedef std::vector<Vertex_handle>                       VertexList;
typedef std::pair<Facet_handle, bool>                    FacetPair;
typedef std::queue<FacetPair>                            FacetQueue;
typedef std::queue<Vertex_handle>                        VertexQueue;
typedef std::map<Vertex_handle, double>                  VertDoubleMap;
typedef boost::associative_property_map<VertDoubleMap>   VertDouble_property_map;
typedef std::map<Vertex_handle, Vector>                  VertVecMap;
typedef boost::associative_property_map< VertVecMap>     VertVec_property_map;
typedef std::map<Facet_handle, Vector>                   FaceVecMap;
typedef boost::associative_property_map< FaceVecMap >    FaceVec_property_map;
typedef std::map<int, double>                            IntDoubleMap;


// UV mapping
typedef CGAL::Unique_hash_map<Vertex_handle, Point_2>           UV_uhm;
typedef boost::associative_property_map<UV_uhm>                 UV_pmap;
typedef CGAL::Unique_hash_map<Vertex_handle, Point>             UV_uhm_3;
typedef boost::associative_property_map<UV_uhm_3>               UV_pmap_3;
typedef CGAL::Aff_transformation_3<Kernel>                      Aff_transformation_3;

// Data Structure
typedef std::pair<Point, Vector>            Point_with_normal;
typedef std::vector<Point_with_normal>      PwnList;
typedef std::vector<Point>                  PointList;
typedef std::vector<float>                  DataList;
typedef std::vector<int>                    IntList;
typedef std::vector<double>                 DoubleList;
typedef std::vector<IntList>                PolyList;
typedef std::vector<Triangle>               TriangleList;
typedef std::vector<bool>                   BoolList;
typedef std::vector<Segment>                SegmentList;
typedef std::vector<Vector>                 VectorList;
typedef std::vector<Point_2>                Point2dList;
typedef std::vector<std::string>            StringList;
typedef std::vector<Segment_2>              Segment2dList;


// AABBTree for triangles
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron>           AABBPrimitive;
typedef CGAL::AABB_traits<Kernel, AABBPrimitive>                       AABBTraits;
typedef CGAL::AABB_tree<AABBTraits>                                    AABBTree;
typedef AABBTree::Point_and_primitive_id                               Point_and_primitive_id;

// AABBTree for segments
typedef SegmentList::iterator                                   SegmentIterator;
typedef CGAL::AABB_segment_primitive<Kernel, SegmentIterator>   SegmentPrimitive;
typedef CGAL::AABB_traits<Kernel, SegmentPrimitive>             SegmentTraits;
typedef CGAL::AABB_tree<SegmentTraits>                          SegmentTree;

// Remeshing
typedef FT (*ValueFunction)(Point);
typedef CGAL::Surface_mesh_default_triangulation_3              Tr;       
typedef CGAL::Complex_2_in_triangulation_3<Tr>                  C2t3;
typedef CGAL::Implicit_surface_3<Kernel, ValueFunction>         Surface_3;

// KDTree
typedef CGAL::Search_traits_3<Kernel>                           KDTraits;
typedef CGAL::Fuzzy_sphere<KDTraits>                            Fuzzy_circle;
typedef CGAL::Kd_tree<KDTraits>                                 KDTree;

// Eigen
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>          EMatrix;
typedef Eigen::VectorXd                                                EVector;
typedef Eigen::SparseMatrix<double>                                    ESMatrix;
typedef std::vector<Eigen::Triplet<double> >                           ESTripleList;
typedef Eigen::SimplicialLDLT<ESMatrix>                                ESolver;
typedef Eigen::Matrix3d     EMatrix3d;
typedef Eigen::Vector3d     EVector3d;
typedef Eigen::Quaterniond  EQuaterniond;

// Geodesic distance
typedef typename CGAL::Surface_mesh_shortest_path_traits<Kernel, Polyhedron>                    Geodesic_traits;
typedef typename Geodesic_traits::Barycentric_coordinates                                       Geodesic_coordinates;
typedef typename boost::property_map<Polyhedron, boost::vertex_external_index_t>::const_type    Vertex_index_map;
typedef typename boost::property_map<Polyhedron, CGAL::halfedge_external_index_t>::const_type   Halfedge_index_map;
typedef typename boost::property_map<Polyhedron, CGAL::face_external_index_t>::const_type       Face_index_map;
typedef typename CGAL::Surface_mesh_shortest_path<  Geodesic_traits, 
                                                    Vertex_index_map,
                                                    Halfedge_index_map,
                                                    Face_index_map >                            Geodesic_tree;
typedef typename Geodesic_tree::Face_location                                                   Face_location;


#ifndef ortho_max
#define ortho_max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef ortho_min
#define ortho_min(a,b)            (((a) < (b)) ? (a) : (b))
#endif


#ifndef M_PI
#define M_PI 3.14159265358979323846264338327950288
#endif


enum SceneObjectType
{
    SceneObjectTypeStart = 0,
    SceneUpperTeeth = 0,
    SceneLowerTeeth,
    SceneMeshEdges,
    SceneUpperInsideRay,
    SceneUpperInsideFunction,
    SceneUpperInsideFaces,
    SceneLowerInsideRay,
    SceneLowerInsideFunction,
    SceneLowerInsideFaces,
    SceneUpperLocalConvexity,
    SceneUpperLaplacianConvexity,
    SceneLowerLocalConvexity,
    SceneLowerLaplacianConvexity,
    SceneUpperConvexityIsoEdges,
    SceneLowerConvexityIsoEdges,
    SceneSelection,
    SceneSegmentation,
    SceneValidatedSegmentation,
    SceneGeodesicDistance,
    ScenePadOutlines,
    SceneToothBound,
    SceneToothCutPoints,
    SceneToothCutPlane,
    ScenePairedFittingPlanes,
    ScenePairedFittingLines,
    ScenePairedFittingPoints,
    ScenePairedPlanes,
    ScenePairedWires,
    ScenePairedWireShapes,
    ScenePairedWireExtrudeShapes,
    ScenePadShapes,
    SceneUpperTeethDist,
    SceneLowerTeethDist,
    SceneObjectTypeCount
};

enum SceneElementFlag
{
    FlagUpper = (1 << 0),
    FlagLower = (1 << 1),
    FlagMesh = (1 << 2),
    FlagInside = (1 << 3),
    FlagConvex = (1 << 4),
    FlagConvexIso = (1 << 5),
    FlagSegment = (1 << 6),
    FlagValid = (1 << 7),
    FlagPadGeodesic = (1 << 8),
    FlagPadOutline = (1 << 9),
    FlagToothCut = (1 << 10),
    FlagPairCut = (1 << 11),
    FlagPairWire = (1 << 12),
    FlagPadShape = (1 << 13),
    FlagTeethDist = (1 << 14),
    // Can go up to (1 << 31)
    FlagAll = 0xFFFFFFFF
};

struct MyPlane
{
	Plane m_plane = Plane(0., 0., 0., 1.);
    Point m_center = CGAL::ORIGIN;
    Point m_right = CGAL::ORIGIN;
};

const std::string UPPER_TEETH = "upper";
const std::string LOWER_TEETH = "lower";
const std::string PATH_SEPERATION = "/";
const std::string DB_EXT = ".db";
} // namespace OrthoLab

#endif
