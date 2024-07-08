#ifndef MY_WIRE_H
#define MY_WIRE_H

#include "types.h"
#include "scene_object.h"

namespace OrthoLab {


class Teeth;
class Wire : public SceneObject
{
private:
    bool        m_init;
    MyPlane     m_fitting_plane;
    Segment     m_fitting_segment;
    Point       m_fitting_point;
    MyPlane     m_cut_plane;
    SegmentList m_polyline;
    PointList   m_intersection_list;
    Polyhedron  m_wire_left_shape;
    Polyhedron  m_wire_left_shape1; // increased 10% for pad construction
    Polyhedron  m_wire_left_shape2; // increased to 2mm for guide convering windows substraction
    Polyhedron  m_wire_middle_shape;
    Polyhedron  m_wire_right_shape;
    Polyhedron  m_wire_right_shape1; // increased 10% for pad construction
    Polyhedron  m_wire_right_shape2; // increased to 2mm for guide convering windows substraction
    
    Polyhedron  m_wire_shape;
    Polyhedron  m_wire_shape_t;
    Polyhedron  m_wire_shape1; // increased 10% for pad construction
    Polyhedron  m_wire_shape1_t; // increased 10% for pad construction
    Polyhedron  m_wire_shape2; // increased to 2mm for guide convering windows substraction
    Polyhedron  m_wire_shape2_t; // increased to 2mm for guide convering windows substraction
    
    Polyhedron  m_wire_extrude_shape;
    Polyhedron  m_wire_extrude_shape_t;
    // 0 transform not set, 1 transform set
    int         m_transform_status = 0;
    Transformation m_transformation; //new Transformation(); // strange CD will happen if change to simple obj
    const Teeth&    m_teeth;
public:

    Wire(const std::string& id, const Teeth& teeth);

    ~Wire();

    void reset();

    std::string get_base_path();
    void save_to_db();
	void load_from_db();

    void set_fitting_plane(Plane& fitting_plane, Point& center, Point& right);

    void set_fitting_segment(Segment& fitting_segment);
    void set_fitting_point(Point& fitting_point);

    void set_cut_plane(Plane& cut_plane, Point& center, Point& right);

    void set_polylines(PointList& intersection);

    void set_flag(bool flag);

    void set_wire_left_shape(Polyhedron& wire_left_shape,int type=0);

    void set_wire_middle_shape(Polyhedron& wire_middle_shape,int type=0);

    void set_wire_right_shape(Polyhedron& wire_right_shape,int type=0);

    void set_wire_shape(Polyhedron& wire_shape,int type=0);

    void set_wire_extrude_shape(Polyhedron& wire_extrude_shape);

    // used for flatten wires
    Transformation& getTransformation() { return m_transformation; };
    void setTransformation(Transformation &transformation)
    {
        //if (m_transformation == nullptr)
        //    m_transformation = new Transformation();
        m_transform_status = 1;
        m_transformation = transformation;
        m_wire_extrude_shape_t.clear();
        m_wire_shape_t.clear();
        m_wire_shape1_t.clear();
        m_wire_shape2_t.clear();
    };

    MyPlane& get_fitting_plane();
    Segment& get_fitting_segment();
    Point& get_fitting_point();
    MyPlane& get_cut_plane();
    SegmentList& get_polylines();
    bool get_flag();

    Polyhedron& get_wire_left_shape(int type=0);
    Polyhedron& get_wire_right_shape(int type=0);
    Polyhedron &get_wire_shape(int type = 0, bool transformed = false);

    Polyhedron& get_wire_extrude_shape(bool transformed = false);

    void update_pair_cut_plane(DataList &fitting_plane_points,
                               DataList &fitting_plane_normals,
                               DataList &fitting_lines,
                               DataList &fitting_points,
                               DataList &cut_plane_points,
                               DataList &cut_plane_normals,
                               DataList &cut_polylines);
}; // end of class Wire

} // namespace OrthoLab

#endif