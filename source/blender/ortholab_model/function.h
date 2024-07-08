#ifndef MY_FUNCTION_H
#define MY_FUNCTION_H

#include "types.h"
#include "teeth.h"

class FileReader;
class FileWriter;

namespace OrthoLab { 

class Project;
class Function
{

public:
    bool        m_upper_init;
    bool        m_lower_init;
    Teeth       m_upper_teeth;
    Teeth       m_lower_teeth;

    Point       m_center;
    double      m_radius;
    AABBTree    m_teeth_tree;

    Polyhedron  m_segment;
    int         m_segment_status; // 0 if non available segment, 1 if upper, -1 if lower

    const Project&  m_prj;


    //----------------------------- Initialization -----------------------------//

    Function(const Project&  m_prj);

    ~Function();

	void save(FileWriter& file);

    void load(FileReader& file);

    Teeth& get_upper_teeth(){ return m_upper_teeth;};
    Teeth& get_lower_teeth(){ return m_lower_teeth;};
    Polyhedron* get_segment(){ return &m_segment;};
    
    void reset();

    void clear_segmentation();

    bool load_upper_teeth_mesh(std::string filename);

    bool load_lower_teeth_mesh(std::string filename);

    bool load_upper_lower_teeth(const std::string &filename1,
        const std::string &filename2);

    void import_segment(const std::string &filename, bool is_upper);

    bool init_upper_teeth_mesh(const float* vertices, size_t vertices_size,
     const int* triangles, size_t triangles_size);

    bool init_lower_teeth_mesh(const float* vertices, size_t vertices_size,
     const int* triangles, size_t triangles_size);

    void compute_bounding_box();

    void update_teeth_tree();

    const double get_radius();
    const Point& get_center();
    int number_of_upper_teeth_faces();
    int number_of_upper_teeth_vertices();
    int number_of_lower_teeth_faces();
    int number_of_lower_teeth_vertices();

    void close_mesh(const std::string &pad, bool shell_solid);
    void Expand_mesh(const std::string &pad,double thickness);
    //----------------------------- Selection -----------------------------//

    bool update_selection(Point camera_pos, Vector camera_dir, PointList& selects, DataList& selected_points);

    bool add_selected_point(const double x, const double y, const double z, 
    PointList& selects, DataList& selected_points);
    //----------------------------- Algorithm -----------------------------//

    bool compute_selected_segmentation(PointList& selected_points, double isovalue);

    int validate_selected_segmentation();

    void compute_geodesic_pad_outlines(int index, bool flag_lower, int discretize_step, double isovalue);

    void compute_circular_pad_outlines(int index, bool flag_lower, int discretize_step, double isovalue, double min_area);

    void recompute_geodesic_pad_outlines(int index, bool flag_lower, double isovalue);

    void compute_tooth_cut_plane(int index, bool flag_lower);

    void compute_paired_wires(int index, bool flag_lower, double distance, double tol, int max_bend, double radius, double offset);

    void compute_coplanar_wires(bool flag_lower, double distance);

    void compute_wire_shapes(int index, bool flag_lower, int geodesic_step, double width, double enlarge_ratio);

    void compute_pad_shapes(int index, bool flag_lower, double pad_height, double wire_height);

    void compute_minimum_surface_pad(int index, bool flag_lower, double pad_surface);

    void compute_teeth_distances();

    //----------------------------- Validation -----------------------------//
    void save_pad_statistics();

    //----------------------------- Visualization -----------------------------//
    static void generate_mesh_visu(Polyhedron &teeth_poly, DataList &teeth_facets,
                                   DataList &teeth_normals);

    void update_teeth_mesh(DataList& teeth_facets, DataList& teeth_normals, bool flag_lower);

    void update_inside_function(DataList& inside_ray_facets,
                                DataList& inside_ray_normals,
                                DataList& inside_ray_colors,
                                DataList& inside_func_colors,
                                DataList& inside_facets,
                                DataList& inside_normals,
                                bool flag_lower);

    void update_convex_function(DataList& local_convex_facets,
                                DataList& local_convex_normals,
                                DataList& local_convex_colors,
                                DataList& convex_func_colors,
                                bool flag_lower);

    void update_convexity_isovalue(DataList& convex_isoedges, double value, bool flag_lower);

    void update_segmentation(DataList& seg_facets, DataList& seg_normals);

    void update_validated_segmentation(int index, bool flag_lower, DataList& seg_facets, 
    DataList& seg_normals);

    void update_tooth_geodesic_function(int index, bool flag_lower, DataList& poly_colors);

    void update_pad_outline(int index, bool flag_lower, DataList& poly_facets, DataList& poly_normals);

    void update_tooth_cut_plane(int index, bool flag_lower, DataList& upper_bound, 
    DataList& upper_points, DataList& cut_plane_points, DataList& cut_plane_normals);

    void update_pair_cut_plane( int index, bool flag_lower, 
                                DataList& fitting_plane_points, 
                                DataList& fitting_plane_normals,
                                DataList& fitting_lines,
                                DataList& fitting_points,
                                DataList& cut_plane_points,
                                DataList& cut_plane_normals,
                                DataList& cut_polylines);

    static void compute_plane_corners(Plane& plane, Point& center, Point& right, 
    PointList& corners, Vector& normal);

    void update_pair_wire_shape( int index, bool flag_lower, 
                                DataList& wire_shape_points, 
                                DataList& wire_shape_normals,
                                DataList& wire_extrude_shape_points, 
                                DataList& wire_extrude_shape_normals);

    void update_pad_shape(  int index, bool flag_lower, 
                            DataList& pad_shape_points, 
                            DataList& pad_shape_normals);

    void update_teeth_distance_colors(bool flag_lower, double threshold, DataList& poly_colors);

}; // end of class Function

} // namespace OrthoLab

#endif
