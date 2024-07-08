#include "wire.h"

#include "types.h"

#include "function.h"
#include "project.h"
#include "save.h"
#include "meshing.h"
#include "teeth.h"
#include "mesh_io.h"
#include "file_services.h"


namespace OrthoLab {

Wire::Wire(const std::string &id, const Teeth &teeth) : 
    m_teeth(teeth),
    SceneObject(id, teeth)
{
    m_init = false;
}

Wire::~Wire()
{
    reset();
    // if (m_transformation)
    //     delete m_transformation;
    // m_transformation = nullptr;
}

void Wire::reset()
{
    m_polyline.clear();
    m_wire_left_shape.clear();
    m_wire_middle_shape.clear();
    m_wire_right_shape.clear();
    m_wire_shape.clear();
    m_wire_extrude_shape.clear();
    m_init = false;
}

std::string Wire::get_base_path()
{
    return m_teeth.get_base_path() +
           PATH_SEPERATION + "w_" + get_id();           
}

void Wire::save_to_db()
{
    std::string wire_key = get_base_path();
    std::string wire_path = ensure_directory_exists(wire_key);
    wire_path += PATH_SEPERATION;
    // save trasformation
    std::string wire_key_transform = wire_key + "transform";
    get_project().get_store().store_transformation(wire_key_transform, m_transformation);
    
    std::cout << "save_to_db_transform: " << wire_key << "transform";
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            std::cout << (m_transformation.m(i, j));
    std::cout << std::endl;
    Transformation transformation2;
    get_project().get_store().retrieve_transformation(wire_key_transform, transformation2);

    std::cout << "load_from_transform: " << wire_key << "transform";
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
            std::cout << (transformation2.m(i, j));
    std::cout << std::endl;
    
    std::string wire_key_intersection = wire_key + "m_intersection_list";
    get_project().get_store().store_pointlist(
        wire_key_intersection, m_intersection_list);
    // std::cout << "save_to_db_transform: " << wire_key << "transform";
    // for (int i = 0; i < 3; i++)
    //     for (int j = 0; j < 4; j++)
    //         std::cout << (m_transformation.m(i, j));
    // std::cout << std::endl;
    mesh_export(get_wire_extrude_shape(), wire_path + "wire_extrude_shape_"+get_id()+".stl");
    mesh_export(get_wire_extrude_shape(true), wire_path + "wire_extrude_shape_t"+get_id()+".stl");
}

void Wire::load_from_db()
{
    std::string wire_key = get_base_path();
    std::string wire_path = ensure_directory_exists(wire_key);
    wire_path += PATH_SEPERATION;

    get_project().get_store().retrieve_transformation(
                         wire_key + "transform", m_transformation);

	PointList intersection_list;
	get_project().get_store().retrieve_pointlist(
		wire_key + "m_intersection_list", intersection_list);

	set_polylines(intersection_list);
    std::string shape_file_path = wire_path + "wire_extrude_shape_" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_wire_extrude_shape);

    shape_file_path = wire_path + "wire_extrude_shape_t" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_wire_extrude_shape_t);

    // std::cout << "load_from_db_transform: "<<wire_key<<"transform";
    // for (int i = 0; i < 3; i++)
    //         for (int j = 0; j < 4; j++)
    //            std::cout << (m_transformation.m(i, j));
    // std::cout << std::endl;

}

void Wire::set_fitting_plane(Plane& fitting_plane, Point& center, Point& right)
{
    m_fitting_plane.m_plane = fitting_plane;
    m_fitting_plane.m_center = center;
    m_fitting_plane.m_right = right;
}

void Wire::set_fitting_segment(Segment& fitting_segment) { m_fitting_segment = fitting_segment; }
void Wire::set_fitting_point(Point& fitting_point) { m_fitting_point = fitting_point; }

void Wire::set_cut_plane(Plane& cut_plane, Point& center, Point& right)
{
    m_cut_plane.m_plane = cut_plane;
    m_cut_plane.m_center = center;
    m_cut_plane.m_right = right;
}

void Wire::set_polylines(PointList& intersection)
{
    if (intersection.size() <= 1)
        return;
    m_polyline.clear();
    m_intersection_list = intersection;
    for (int i = 0; i < intersection.size() - 1; i++)
        m_polyline.push_back(Segment(intersection[i], intersection[i + 1]));
}

void Wire::set_flag(bool flag) { m_init = flag; }

void Wire::set_wire_left_shape(Polyhedron& wire_left_shape,int type)
{
    if (type == 0)
        copy_polyehdral_surface(wire_left_shape, m_wire_left_shape);
    else if (type == 1)
        copy_polyehdral_surface(wire_left_shape, m_wire_left_shape1);
    else if (type == 2)
        copy_polyehdral_surface(wire_left_shape, m_wire_left_shape2);
}

void Wire::set_wire_middle_shape(Polyhedron& wire_middle_shape,int type)
{
    copy_polyehdral_surface(wire_middle_shape, m_wire_middle_shape);
}

void Wire::set_wire_right_shape(Polyhedron& wire_right_shape,int type)
{ 
    if (type == 0)
        copy_polyehdral_surface(wire_right_shape, m_wire_right_shape);
    else if (type == 1)
        copy_polyehdral_surface(wire_right_shape, m_wire_right_shape1);
    else if (type == 2)
        copy_polyehdral_surface(wire_right_shape, m_wire_right_shape2);
}

void Wire::set_wire_shape(Polyhedron& wire_shape,int type)
{
    if (type == 0)
        copy_polyehdral_surface(wire_shape, m_wire_shape);
    else if (type == 1)
        copy_polyehdral_surface(wire_shape, m_wire_shape1);
    else if (type == 2)
        copy_polyehdral_surface(wire_shape, m_wire_shape2);
}

void Wire::set_wire_extrude_shape(Polyhedron& wire_extrude_shape)
{
    copy_polyehdral_surface(wire_extrude_shape, m_wire_extrude_shape);
}

MyPlane& Wire::get_fitting_plane() { return m_fitting_plane; }
Segment& Wire::get_fitting_segment() { return m_fitting_segment; }
Point& Wire::get_fitting_point() { return m_fitting_point; }
MyPlane& Wire::get_cut_plane() { return m_cut_plane; }
SegmentList& Wire::get_polylines() { return m_polyline; }
bool Wire::get_flag() { return m_init; }

Polyhedron &Wire::get_wire_left_shape(int type)
{
    if (type == 1)
        return m_wire_left_shape1;
    else if (type == 2)
        return m_wire_left_shape2;
    return m_wire_left_shape;
}
Polyhedron &Wire::get_wire_right_shape(int type)
{
    if (type == 1)
        return m_wire_right_shape1;
    else if (type == 2)
        return m_wire_right_shape2;
    return m_wire_right_shape;
}
Polyhedron &Wire::get_wire_shape(int type, bool transformed)
{
    Polyhedron* shape = &m_wire_shape;
    if (type == 1)
        shape = &m_wire_shape1;
    else if (type == 2)
        shape = &m_wire_shape2;
    if(transformed && m_transform_status == 1)
    {
        if(m_wire_shape_t.empty())
        {
            copy_polyehdral_surface(m_wire_shape, m_wire_shape_t);
            transform_mesh(m_wire_shape_t, m_transformation);

            copy_polyehdral_surface(m_wire_shape1, m_wire_shape1_t);
            transform_mesh(m_wire_shape1_t, m_transformation);

            copy_polyehdral_surface(m_wire_shape2, m_wire_shape2_t);
            transform_mesh(m_wire_shape2_t, m_transformation);
        }
        shape = &m_wire_shape_t;
        if (type == 1)
            shape = &m_wire_shape1_t;
        else if (type == 2)
            shape = &m_wire_shape2_t;
    }
    return *shape;
}

Polyhedron& Wire::get_wire_extrude_shape(bool transformed) 
{ 
    if(transformed && m_transform_status == 1)
    {
        copy_polyehdral_surface(m_wire_extrude_shape, m_wire_extrude_shape_t);
        transform_mesh(m_wire_extrude_shape_t, m_transformation);
        return m_wire_extrude_shape_t;
    }
    return m_wire_extrude_shape; 
}

void Wire::update_pair_cut_plane(DataList &fitting_plane_points,
                                 DataList &fitting_plane_normals,
                                 DataList &fitting_lines,
                                 DataList &fitting_points,
                                 DataList &cut_plane_points,
                                 DataList &cut_plane_normals,
                                 DataList &cut_polylines)
{
    Wire& wire = *this;// (flag_lower) ? m_lower_teeth.get_wire(index) : m_upper_teeth.get_wire(index);
    if(!wire.get_flag())
        return;
    
    MyPlane& wire_fitting_plane = wire.get_fitting_plane();
    Segment& wire_fitting_seg = wire.get_fitting_segment();
    Point& wire_fitting_point = wire.get_fitting_point();
    MyPlane& wire_cut_plane = wire.get_cut_plane();
    SegmentList& wire_polylines = wire.get_polylines();

    fitting_points.push_back((float)wire_fitting_point.x());
    fitting_points.push_back((float)wire_fitting_point.y());
    fitting_points.push_back((float)wire_fitting_point.z());

    fitting_lines.push_back((float)wire_fitting_seg.source().x());
    fitting_lines.push_back((float)wire_fitting_seg.source().y());
    fitting_lines.push_back((float)wire_fitting_seg.source().z());
    fitting_lines.push_back((float)wire_fitting_seg.target().x());
    fitting_lines.push_back((float)wire_fitting_seg.target().y());
    fitting_lines.push_back((float)wire_fitting_seg.target().z());

    for (int i = 0; i < wire_polylines.size(); i++)
    {
        Segment seg = wire_polylines[i];
        Point source = seg.source();
        Point target = seg.target();
        cut_polylines.push_back((float)source.x());
        cut_polylines.push_back((float)source.y());
        cut_polylines.push_back((float)source.z());
        cut_polylines.push_back((float)target.x());
        cut_polylines.push_back((float)target.y());
        cut_polylines.push_back((float)target.z());
    }

    Vector normal;
    PointList corners;
    Function::compute_plane_corners(wire_fitting_plane.m_plane, wire_fitting_plane.m_center, wire_fitting_plane.m_right, corners, normal);
    for (int i = 0; i < corners.size(); i++)
    {
        Point p = corners[i];
        fitting_plane_points.push_back((float)p.x());
        fitting_plane_points.push_back((float)p.y());
        fitting_plane_points.push_back((float)p.z());
        fitting_plane_normals.push_back((float)normal.x());
        fitting_plane_normals.push_back((float)normal.y());
        fitting_plane_normals.push_back((float)normal.z());
    }

    corners.clear();
    Function::compute_plane_corners(wire_cut_plane.m_plane, wire_cut_plane.m_center, wire_cut_plane.m_right, corners, normal);
    for (int i = 0; i < corners.size(); i++)
    {
        Point p = corners[i];
        cut_plane_points.push_back((float)p.x());
        cut_plane_points.push_back((float)p.y());
        cut_plane_points.push_back((float)p.z());
        cut_plane_normals.push_back((float)normal.x());
        cut_plane_normals.push_back((float)normal.y());
        cut_plane_normals.push_back((float)normal.z());
    }
}

} // namespace OrthoLab
