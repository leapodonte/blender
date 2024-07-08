#ifndef MY_TEETH_H
#define MY_TEETH_H

#include "types.h"
#include "scene_object.h" 
#include "manu_asm_guide.h"

class FileReader;
class FileWriter;


namespace OrthoLab {
class Tooth;
class Wire;

class Teeth : public SceneObject
{
public:
    typedef typename std::vector< std::unique_ptr<Tooth>> ToothList;
    typedef typename std::vector<std::unique_ptr<Wire>>  WireList;

private:
    Polyhedron      m_teeth_poly;
    VertDoubleMap   m_teeth_distance_map;
    Polyhedron      m_teeth_poly_closed;
    Polyhedron      m_teeth_poly_expanded;
    VertIntMap      m_vert_ind_map;
    VertDoubleMap   m_vert_ray_map;
    VertDoubleMap   m_vert_smooth_ray_map;
    VertDoubleMap   m_vert_inside_func_map;
    
    Polyhedron      m_teeth_inside;
    VertIntMap      m_inside_vert_ind_map;
    VertDoubleMap   m_inside_vert_convex_map;
    VertDoubleMap   m_inside_vert_smooth_convex_map;
    VertDoubleMap   m_inside_vert_convex_func_map;
    AABBTree        m_inside_tree;
    KDTree          m_inside_boundary_tree;

    ToothList       m_tooth;
    WireList        m_wires;
    PointList       m_flattened_wire;

    ManuAsmGuide    m_manu_asm_guide;

public:

    //----------------------------- Initialization -----------------------------//

    Teeth(const std::string& id, const Project& prj);

    ~Teeth();

    void save(FileWriter& file);
    void load(FileReader& file);
    
	void save_to_db();
	void load_from_db();
    std::string get_base_path() const; // for db save and load

    void close_teeth_mesh();
    void expandTeethMesh(double dist);
    void save_teeth_model(const std::string &filename);

    void flattenWires(double dist);

    void reset();

    void reset_inside();

    void reset_convexity();

    bool load_teeth_mesh(std::string filename);

    bool init_teeth_mesh(const float* vertices, size_t vertices_size, const int* triangles, size_t triangles_size);

    Polyhedron& get_teeth_poly();
    Polyhedron& get_closed_teeth_poly();
    VertDoubleMap& get_teeth_distance_map();
    VertDoubleMap& get_smoothed_ray_map();
    VertDoubleMap& get_vertex_inside_function_map();
    Polyhedron& get_inside_teeth_poly();
    VertDoubleMap& get_smoothed_convex_map();
    VertDoubleMap& get_vertex_convex_function_map();
    ManuAsmGuide& get_manu_asm_guide(){return m_manu_asm_guide;};
    int size_of_facets();
    int size_of_vertices();
    int size_of_inside_facets();
    int size_of_inside_vertices();

    void read_polygon_vertices(PointList& points);

    void initialize_vertex_indices();

    double get_zmax();
 
    // void exportTeethMesh(QJsonObject &meshObject);

    void compute_teeth_distance(Polyhedron& other_teeth);

    //----------------------------- Inside Function -----------------------------//

    bool compute_laplacian_based_inside(AABBTree& teeth_tree, Point& center, int smooth_range, double lambda);

    void compute_intersection_ray(AABBTree& teeth_tree, Point& center);

    double compute_distance_to_inside_teeth(Point& point);

    void extract_inside_from_function(double value);

    void initialize_inside_vertex_index();

    void initialize_inside_boundary_kdtree();

    void fill_holes_in_vertmap(Polyhedron& poly, VertDoubleMap& vert_map, double isovalue, VertDoubleMap& new_vert_map);

    //----------------------------- Convexity Function -----------------------------//

    bool compute_laplacian_based_convexity(double radius_ratio, double bounding_radius, int smooth_range, double lambda);

    void compute_convexity_from_sampling(double radius_ratio, double bounding_radius);

    //----------------------------- Segmentation -----------------------------//

    void compute_selected_segmentation(PointList& selected_points, double isovalue, Polyhedron& segment_poly);

    bool validate_selected_segmentation(Polyhedron& segment_poly, bool flag_lower);

    //----------------------------- Tooth Function -----------------------------//

    Polyhedron& get_tooth_poly(int index);
    VertDoubleMap& get_geodesic_distance_map(int index);
    Polyhedron& get_pad_outline_poly(int index);

    HEdgeSet& get_upper_bound_set(int index);
    PointList& get_cut_points(int index);
    Plane& get_cut_plane(int index);

    Wire& get_wire(int index);

    Tooth& get_tooth(int index);

    void compute_geodesic_pad_outlines(int index, int discretize_step, double isovalue, bool flag_smooth);

    void compute_circular_pad_outlines(int index, int discretize_step, double isovalue, double min_area);

    void recompute_geodesic_pad_outlines(int index, double isovalue, bool flag_smooth);

    void compute_tooth_cut_plane(int index);

    void compute_pad_shapes(int index, bool flag_left, bool flag_right, double pad_height, double wire_height, double wire_radius);

    //----------------------------- Wire Function -----------------------------//

    void compute_paired_wires(int index, double distance, double tolerance, int max_bend, double radius, double offset);

    void compute_coplanar_wires(double distance);

    static void copy_boundary_edges_to_segments(HEdgeSet& bound_hedge_set, PointList& points, SegmentList& segments);

    static void project_points_on_line(PointList& candidate_points, Line& fitting_line, Point& left_point, Point& right_point);

    void compute_smooth_wire(double tolerance, int max_bend, PointList& polylines, Plane& wire_plane);

    void compute_smooth_intersection(std::vector<Point2dList>& multi_proj_polyline_pts, std::vector<Line_2>& multi_fitting_lines, Point2dList& polylines_2d);

    bool split_fitting_line_with_criterion(std::vector<Point2dList>& multi_proj_polyline_pts, std::vector<Line_2>& multi_fitting_lines, DoubleList& multi_distortions, double tol, int max_bend);

    double compute_fitting_line_2d(Point2dList& proj_seg_points, Line_2& proj_line);

    int add_split_point_for_polyline(Point2dList& proj_points, Line_2& fitting_line);

    void compute_pair_plane(double distance, Point& plane_center, Vector& line_direction, Point& pair_center, Plane& fitting_plane, Plane& pair_plane);

    // void compute_plane_mesh(double scale, Plane& plane, Point& center, Point& right, Polyhedron& plane_mesh);

    static bool compute_intersection_plane_line(Plane& plane, Line& line, Point& point);

    //void compute_wire_shapes(int index, int discretize_step, double width);
    void compute_wire_shapes(int index, int discretize_step, double radius);

    void compute_enlarge_wire_shapes(int index, int discretize_step, double radius, double enlarge_ratio);

    void clip_mesh_by_plane(Polyhedron& mesh, Plane& plane, Point& center);

    void clip_volume_by_plane(Polyhedron& mesh, Plane& plane, Point& center);

    static bool extract_wire_shape_for_halftooth(  double width, Plane& cut_plane, Geodesic_tree& geodesic_tree, AABBTree& geodesic_aabb_tree, Polyhedron& halftooth,  // input
                                            Polyhedron& halfwire, Point& upper_bound, Point& lower_bound);
    static bool extract_wire_shape_for_halftooth(  double width, Geodesic_tree& geodesic_tree, AABBTree& geodesic_aabb_tree, Polyhedron& halftooth,  // input
                                            Polyhedron& halfwire);

    void find_vnormal_from_adjacent_mesh(Polyhedron& halfwire, VertexList& vert_handles, VectorList& vert_normals, double threshold = 1e-6);

    //----------------------------- surface expansion pad Function -----------------------------//
    void compute_minimum_surface_pad(int index, double surface_pad);

    //----------------------------- Inside / Convexity Help Function -----------------------------//

    static void compute_smoothed_field(int smooth_range, VertDoubleMap& init_map, VertDoubleMap& value_map);

    double locate_and_evaluate_function(Point& query, Facet_handle face, VertDoubleMap& vert_func_map);

    void barycentric_coordinates(Point& p, Point& pa, Point& pb, Point& pc, double& a, double& b, double& c);

    /******************* test **********************/
	void save_validated_segmentation(const std::string& filename);
	void compute_manu_guide(bool flag_lower, const int start_index, const int end_index, const bool flattened);
	void compute_bonding_surface(const std::string& pad, const std::string& padoutline, const std::string& padout);
    void load_validated_segmentation(const std::string& filename);
	void load_guide(const std::string& filename);
    void get_pad_statistics(StringList& names, DoubleList& left, DoubleList& right);
}; // end of class Teeth

} // namespace OrthoLab

#endif
