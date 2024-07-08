#ifndef MY_TOOTH_H
#define MY_TOOTH_H

#include "types.h"
#include "pad.h"
#include "scene_object.h"

namespace OrthoLab {
class Teeth;
class Tooth : public SceneObject
{
private:
    std::string     m_name;
    Polyhedron      m_tooth_poly;
    Polyhedron      m_tooth_poly_t;
    Point           m_center;
    VertIntMap      m_vert_ind_map;
    VertDoubleMap   m_vert_geodesic_map;
    VertDoubleMap   m_vert_smooth_geodesic_map;

    HEdgeSet        m_upper_bound;
    VertexSet       m_upper_points;
    PointList       m_cut_points;
    Plane           m_cut_plane;

    Polyhedron      m_pad_outlines;
    Polyhedron      m_pad_outlines_t;
    Polyhedron      m_pad_shape;
    Polyhedron      m_pad_shape_t;
    bool            m_pad_split_flag;
    Polyhedron      m_pad_left;
    Polyhedron      m_pad_left_t;
    Polyhedron      m_pad_right;
    Polyhedron      m_pad_right_t;
    Polyhedron      m_pad_outlines_left;
    Polyhedron      m_pad_outlines_right;
    Polyhedron      m_pad_outlines_left_t;
    Polyhedron      m_pad_outlines_right_t;

    double          m_minimun_surface;
    Pad             m_pad;
    Vector          m_pad_extrusion_direction;

    Polyhedron      m_bonding_surface;
    const Teeth&    m_teeth;

    Transformation m_transform_left;
    Transformation m_transform_right;
    // 0 transform not set, 1 transform set
    int            m_transform_status = 0; 
public:

    Tooth(const std::string& id, const Teeth& prj);

    ~Tooth();


	void save_to_db();

	void load_from_db();

    std::string get_base_path() const;

    Pad& get_pad(){ return m_pad;}

    void reset();

    void reset_cut_plane();

    void reset_pad();

    void set_tooth_mesh(Polyhedron& seg_poly, Polyhedron& teeth_poly);

    void set_name(std::string name);

    double compute_length_for_hole(Halfedge_handle he, Polyhedron& mesh);

    void initialize_vertex_indices();

    std::string get_name();
    bool get_init_flag();
    const Point& get_center();
    Polyhedron& get_tooth_poly(bool transformed=false);
    VertDoubleMap& get_geodesic_distance_map();
    Polyhedron& get_pad_outline_poly(bool transformed=false);
    HEdgeSet& get_upper_bound_set();
    PointList& get_cut_points();
    Plane& get_cut_plane();
    bool get_split_flag();
    Polyhedron& get_pad_shape_poly(bool transformed=false);
    Polyhedron& get_left_pad_shape_poly(bool transformed=false);
    Polyhedron& get_right_pad_shape_poly(bool transformed = false);
    Polyhedron& get_left_pad_outline_poly(bool transformed = false);
    Polyhedron& get_right_pad_outline_poly(bool transformed = false);
    Polyhedron& get_bonding_surface();
    
    // int type
    // 0 => origin
    // 1 left pad extrusion direction transformed
    // 2 right pad extrusion direction transformed
    Vector get_pad_extrusion_direction(const int type = 0) 
    { 
        if (type == 1)
        {
            return m_transform_left.transform(m_pad_extrusion_direction);
        }
        else if (type == 2)
        {
            return m_transform_right.transform(m_pad_extrusion_direction);
        }
        //  if (type == 0)
        return m_pad_extrusion_direction;
    };

    Transformation get_left_transformation()
    {return m_transform_left;}
    Transformation get_right_transformation()
    {return m_transform_right;}
    void set_left_transformation(const Transformation t);
    void set_right_transformation(const Transformation t)
    {m_transform_right = t;}

    bool compute_geodedic_pad_outlines(int discretize_step, double isovalue, Polyhedron& teeth_poly, bool flag_smooth);

    bool compute_circular_pad_outlines(int discretize_step, double isovalue, double min_area);

    void compute_geodesic_distance(int discretize_step);

    bool update_geodesic_pad_outlines(double isovalue, Polyhedron& teeth_poly, bool flag_smooth);

    void smooth_geodesic_pad_outlines(Polyhedron& teeth_poly);

    Point compute_centroid(Polyhedron& poly);

    bool update_circular_pad_outlines(double isovalue, double radius);

    bool compute_tooth_cut_plane(KDTree& teeth_boundary_tree, AABBTree& teeth_tree);

    void compute_geodesic_center(Point& geodesic_center);

    void compute_minimum_surface_pad(const double surface);

    void compute_boundary_halfedges(KDTree& boundary_tree, SegmentList& upper_bound_segs);

    void compute_boundary_midpoint(Line& fitting_line, Point& upper_bound_center);

    void compute_pad_shape(Polyhedron& wire_shape, Polyhedron& enlarge_wire_shape, double pad_height, double wire_height, Point& teeth_center, double wire_radius);

    void compute_splited_pad_shape( Polyhedron& left_wire_shape, Polyhedron& left_enlarge_wire_shape, 
                                    Polyhedron& right_wire_shape, Polyhedron& right_enlarge_wire_shape, 
                                    double pad_height, double wire_height, Point& teeth_center, double wire_radius);
    
    void compute_bonding_surface(const std::string& pad, const std::string& padoutline, const std::string& padout);

    void save_pad(const std::string &filename);

    void get_pad_statistics(StringList& names, DoubleList& left, DoubleList& right);

}; // end of class Tooth

} // namespace OrthoLab

#endif
