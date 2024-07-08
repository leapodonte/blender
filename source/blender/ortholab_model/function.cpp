#include "function.h"

#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include "wire.h"
#include "meshing.h"
#include "mesh_io.h"
#include "ramp.h"
#include "save.h"
#include "project.h" 
#include "tooth.h" 

#define USE_OPENMP

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include "spdlog/spdlog.h"
#include "timer.h"

namespace OrthoLab {

//----------------------------- Initialization -----------------------------//

Function::Function(const Project &prj) : m_prj(prj),
                                         m_upper_init(false),
                                         m_lower_init(false),
                                         m_radius(1.),
                                         m_center(Point(0., 0., 0.)),
                                         m_upper_teeth(UPPER_TEETH, prj),
                                         m_lower_teeth(LOWER_TEETH, prj)
{
    
}

Function::~Function()
{
    reset();
}

void Function::reset()
{
    m_upper_init = false;
    m_lower_init = false;
    m_center = Point(0., 0., 0.);
    m_radius = 1.;
    m_upper_teeth.reset();
    m_lower_teeth.reset();
    m_teeth_tree.clear();
    clear_segmentation();
}

void Function::clear_segmentation()
{
    m_segment.clear();
    m_segment_status = 0;
}

void write_point(FileWriter& file, const Point& p)
{
    file.write(p.x());
    file.write(p.y());
    file.write(p.z());
}

void read_point(FileReader& file, Point& p)
{
    file.read((double&)p.x());
    file.read((double&)p.y());
    file.read((double&)p.z());
}

void write_point(FileWriter& file, const char* key, const Point& p)
{
    file.write_key_value_cb(key, [&](FileWriter& file) { write_point(file, p); });
}

void read_point(FileReader& file, const char* key, Point& p)
{
    if (!file.empty && file.current_key == key)
        read_point(file, p);
}

void Function::save(FileWriter& file)
{
    m_center = Point(321, 432, 543);
    file.write_key_value("m_upper_init", m_upper_init);
    file.write_key_value("m_lower_init", m_lower_init);
    file.write_key_value_cb("m_upper_teeth", std::function([&](FileWriter& file) { m_upper_teeth.save(file); }));
    file.write_key_value_cb("m_lower_teeth", std::function([&](FileWriter& file) { m_lower_teeth.save(file); }));
    write_point(file, "m_center", m_center);
    file.write_key_value("m_radius", m_radius);
    file.write_key_value("m_segment_status", m_segment_status);

    // TODO: AABBTree    m_teeth_tree;

    // TODO: Polyhedron  m_segment;
}

void Function::load(FileReader& file)
{
    while (!file.empty)
    {
        file.next_key();

        file.read_key_value("m_upper_init", m_upper_init);
        file.read_key_value("m_lower_init", m_lower_init);
        file.read_key_value_cb("m_upper_teeth", std::function([&](FileReader& file) { m_upper_teeth.load(file); }));
        file.read_key_value_cb("m_lower_teeth", std::function([&](FileReader& file) { m_lower_teeth.load(file); }));
        read_point(file, "m_center", m_center);
        file.read_key_value("m_radius", m_radius);
        file.read_key_value("m_segment_status", m_segment_status);
    }

    // TODO: AABBTree    m_teeth_tree;

    // TODO: Polyhedron  m_segment;
}

bool Function::load_upper_teeth_mesh(std::string filename)
{
    spdlog::trace("Loading upper teeth!");
    m_upper_teeth.reset();
    m_upper_init = m_upper_teeth.load_teeth_mesh(filename);

    if (m_upper_init)
    {
        compute_bounding_box();
        update_teeth_tree();
    }

    return m_upper_init;
}

bool Function::load_lower_teeth_mesh(std::string filename)
{
    spdlog::trace("Loading lower teeth!");
    m_lower_teeth.reset();
    m_lower_init = m_lower_teeth.load_teeth_mesh(filename);

    if (m_lower_init)
    {
        compute_bounding_box();
        update_teeth_tree();
    }

    return m_lower_init;
}

bool Function::load_upper_lower_teeth(const std::string &filename1,
                                      const std::string &filename2)
{
    spdlog::trace("Loading upper and lower teeth!");
    m_upper_teeth.reset();
    m_lower_teeth.reset();

#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            {
                m_upper_init = m_upper_teeth.load_teeth_mesh(filename1);
            }
#pragma omp section
            m_lower_init = m_lower_teeth.load_teeth_mesh(filename2);
        }
    }
    if (m_upper_init && m_lower_init)
    {
        compute_bounding_box();
        update_teeth_tree();
    }
    return true;
}

void Function::import_segment(const std::string &filename, bool is_upper)
{
    m_segment_status = 0;
    if (!mesh_import(filename, m_segment)) return; 
    m_segment_status = -1;
    if(is_upper)
        m_segment_status = 1;
}

bool Function::init_upper_teeth_mesh(const float* vertices, size_t vertices_size, const int* triangles, size_t triangles_size)
{
    spdlog::trace("Init upper teeth!");
    //std::cout << "Init upper teeth!" << std::endl;
    m_upper_teeth.reset();
    m_upper_init = m_upper_teeth.init_teeth_mesh(vertices, vertices_size, triangles, triangles_size);

    if (m_upper_init)
    {
        compute_bounding_box();
        update_teeth_tree();
    }

    return m_upper_init;
}

bool Function::init_lower_teeth_mesh(const float* vertices, size_t vertices_size, const int* triangles, size_t triangles_size)
{
    spdlog::trace("Init lower teeth!");
    //std::cout << "Init lower teeth!" << std::endl;
    m_lower_teeth.reset();
    m_lower_init = m_lower_teeth.init_teeth_mesh(vertices, vertices_size, triangles, triangles_size);

    if (m_lower_init)
    {
        compute_bounding_box();
        update_teeth_tree();
    }

    return m_lower_init;
}

void Function::compute_bounding_box()
{
    PointList points;

    if (m_upper_init)
    {
        m_upper_teeth.read_polygon_vertices(points);
    }

    if (m_lower_init)
    {
        m_lower_teeth.read_polygon_vertices(points);
    }

    if (points.size() > 0)
    {
        Bbox bounding_box = CGAL::bbox_3(points.begin(), points.end());
        m_center = Point(0.5 * (bounding_box.xmin() + bounding_box.xmax()),
            0.5 * (bounding_box.ymin() + bounding_box.ymax()),
            0.5 * (bounding_box.zmin() + bounding_box.zmax()));
        Point corner(bounding_box.xmin(), bounding_box.ymin(), bounding_box.zmin());
        m_radius = std::sqrt(CGAL::squared_distance(m_center, corner));
    }
}

void Function::update_teeth_tree()
{
    m_teeth_tree.clear();

    if (m_upper_init)
    {
        Polyhedron& teeth_poly = m_upper_teeth.get_teeth_poly();
        m_teeth_tree.insert(faces(teeth_poly).first, faces(teeth_poly).second, teeth_poly);
    }

    if (m_lower_init)
    {
        Polyhedron& teeth_poly = m_lower_teeth.get_teeth_poly();
        m_teeth_tree.insert(faces(teeth_poly).first, faces(teeth_poly).second, teeth_poly);
    }

    if (!m_teeth_tree.empty())
        m_teeth_tree.accelerate_distance_queries();
}

const double Function::get_radius() { return m_radius; }
const Point& Function::get_center() { return m_center; }
int Function::number_of_upper_teeth_faces() { return m_upper_teeth.size_of_facets(); }
int Function::number_of_upper_teeth_vertices() { return m_upper_teeth.size_of_vertices(); }
int Function::number_of_lower_teeth_faces() { return m_lower_teeth.size_of_facets(); }
int Function::number_of_lower_teeth_vertices() { return m_lower_teeth.size_of_vertices(); }

//----------------------------- Selection -----------------------------//

bool Function::update_selection(Point camera_pos, Vector camera_dir, PointList& selects, DataList& selected_points)
{
    if (m_teeth_tree.empty())
        return false;

    Ray ray(camera_pos, camera_dir);
    auto intersection = m_teeth_tree.first_intersection(ray);
    if (intersection)
    {
        if (boost::get<Point>(&(intersection->first))) {
            Point select = *(boost::get<Point>(&(intersection->first)));
            selects.push_back(select);
            selected_points.push_back((float)select.x());
            selected_points.push_back((float)select.y());
            selected_points.push_back((float)select.z());
            std::cout << select.x() << "," << select.y() << "," << select.z() << std::endl;
            return true;
        }
    }

    return false;
}
bool Function::add_selected_point(const double x, const double y, const double z, PointList& selects, DataList& selected_points)
{
	Point select(x, y, z);
	selects.push_back(select);
	selected_points.push_back((float)select.x());
	selected_points.push_back((float)select.y());
	selected_points.push_back((float)select.z());
	std::cout << select.x() << "," << select.y() << "," << select.z() << std::endl;
	return true;
}
//----------------------------- Algorithm -----------------------------//

bool Function::compute_selected_segmentation(PointList& selected_points, double isovalue)
{
    if (selected_points.size() == 0)
    {
        // std::cout << "No available selection!" << std::endl;
        spdlog::error("No available selection!");
        return false;
    }

    double upper_dist = m_upper_teeth.compute_distance_to_inside_teeth(selected_points[0]);
    double lower_dist = m_lower_teeth.compute_distance_to_inside_teeth(selected_points[0]);

    if (upper_dist < 0 || lower_dist < 0)
    {
        //std::cout << "The inside surfaces are not available, please compute the inside surface!" << std::endl;
        spdlog::error("The inside surfaces are not available, please compute the inside surface!");
        return false;
    }

    clear_segmentation();
    spdlog::trace("Begin extracting segmentation!");
    //std::cout << "Begin extracting segmentation!" << std::endl;
    bool flag_lower = (upper_dist < lower_dist) ? false : true;
    if (!flag_lower)
        m_upper_teeth.compute_selected_segmentation(selected_points, isovalue, m_segment);
    else
        m_lower_teeth.compute_selected_segmentation(selected_points, isovalue, m_segment);

    if (m_segment.size_of_vertices() == 0)
    {
        //std::cout << "Segmentation failed!" << std::endl;
        spdlog::error("Segmentation failed!");
        return false;
    }

    m_segment_status = (flag_lower) ? -1 : 1;

    return true;
}

int Function::validate_selected_segmentation()
{
    if (m_segment_status == 0)
    {
        //std::cout << "No available segmentation!" << std::endl;
        spdlog::error("No available selection!");
        return 0;
    }

    int curr_stat = m_segment_status;
    if (curr_stat == -1)
        m_lower_teeth.validate_selected_segmentation(m_segment, true);
    else
        m_upper_teeth.validate_selected_segmentation(m_segment, false);

    clear_segmentation();

    return curr_stat;
}

void Function::compute_geodesic_pad_outlines(int index, bool flag_lower, int discretize_step, double isovalue)
{
    if (flag_lower)
        m_lower_teeth.compute_geodesic_pad_outlines(index, discretize_step, isovalue, true);
    else
        m_upper_teeth.compute_geodesic_pad_outlines(index, discretize_step, isovalue, false);
}

void Function::compute_circular_pad_outlines(int index, bool flag_lower, int discretize_step, double isovalue, double min_area)
{
    if (flag_lower)
        m_lower_teeth.compute_circular_pad_outlines(index, discretize_step, isovalue, min_area);
    else
        m_upper_teeth.compute_circular_pad_outlines(index, discretize_step, isovalue, min_area);
}

void Function::recompute_geodesic_pad_outlines(int index, bool flag_lower, double isovalue)
{
    if (flag_lower)
        m_lower_teeth.recompute_geodesic_pad_outlines(index, isovalue, true);
    else
        m_upper_teeth.recompute_geodesic_pad_outlines(index, isovalue, false);
}

void Function::compute_tooth_cut_plane(int index, bool flag_lower)
{
    if (flag_lower)
        m_lower_teeth.compute_tooth_cut_plane(index);
    else
        m_upper_teeth.compute_tooth_cut_plane(index);
}

void Function::compute_paired_wires(int index, bool flag_lower, double distance, double tol, int max_bend, double radius, double offset)
{
    if (flag_lower)
        m_lower_teeth.compute_paired_wires(index, distance, tol, max_bend, radius, offset);
    else
        m_upper_teeth.compute_paired_wires(index, distance, tol, max_bend, radius, offset);
}

void Function::compute_coplanar_wires(bool flag_lower, double distance)
{
    if (flag_lower)
        m_lower_teeth.compute_coplanar_wires(distance);
    else
        m_upper_teeth.compute_coplanar_wires(distance);
}

void Function::compute_wire_shapes(int index, bool flag_lower, int geodesic_step, double width, double enlarge_ratio)
{
    if (flag_lower)
    {
         m_lower_teeth.compute_wire_shapes(index, geodesic_step, width);
         m_lower_teeth.compute_enlarge_wire_shapes(index, geodesic_step, width, enlarge_ratio);
    }
    else
    {
        m_upper_teeth.compute_wire_shapes(index, geodesic_step, width);
        m_upper_teeth.compute_enlarge_wire_shapes(index, geodesic_step, width, enlarge_ratio);
    }
        
}

// void Function::compute_pad_shapes(int index, bool flag_lower, double pad_height, double wire_height)
// {
//     if (flag_lower)
//         m_lower_teeth.compute_pad_shapes(index, pad_height, wire_height);
//     else
//         m_upper_teeth.compute_pad_shapes(index, pad_height, wire_height);
// }

void Function::compute_minimum_surface_pad(int index, bool flag_lower, double pad_surface)
{
    if (flag_lower)
        m_lower_teeth.compute_minimum_surface_pad(index, pad_surface);
    else
        m_upper_teeth.compute_minimum_surface_pad(index, pad_surface);
}

void Function::compute_teeth_distances()
{
    if (!m_upper_init || !m_upper_init)
    {
        spdlog::error("Teeth not loaded!");
        //std::cout << "Teeth not loaded!" << std::endl;
        return;
    }

    Polyhedron& lower_teeth_poly = m_lower_teeth.get_teeth_poly();
    Polyhedron& upper_teeth_poly = m_upper_teeth.get_teeth_poly();
    m_lower_teeth.compute_teeth_distance(upper_teeth_poly);
    m_upper_teeth.compute_teeth_distance(lower_teeth_poly);
}

//----------------------------- Visualization -----------------------------//

void Function::save_pad_statistics()
{
    StringList names;
    DoubleList left_areas, right_areas;

    m_lower_teeth.get_pad_statistics(names, left_areas, right_areas);
    m_upper_teeth.get_pad_statistics(names, left_areas, right_areas);

    std::ofstream myfile;
    myfile.open("pad_stats.csv");
    myfile << "Tooth,Left,Right\n";
    for(int i = 0; i < names.size(); i++)
    {
        myfile << names[i] << "," << left_areas[i] << "," << right_areas[i] << ",\n";
        std::cout << names[i] << "," << left_areas[i] << "," << right_areas[i] << ",\n";
    }
        

    myfile.close();
}

//----------------------------- Visualization -----------------------------//

void Function::generate_mesh_visu(Polyhedron &teeth_poly, DataList &teeth_facets,
                                  DataList &teeth_normals)
{
    FaceVecMap fnormals;
    CGAL::Polygon_mesh_processing::compute_face_normals(teeth_poly, boost::make_assoc_property_map(fnormals));

    teeth_facets.reserve(teeth_poly.size_of_facets() * 9);
    teeth_normals.reserve(teeth_poly.size_of_facets() * 9);

    for (Facet_iterator face = teeth_poly.facets_begin(); face != teeth_poly.facets_end(); ++face) {
        Halfedge_facet_circulator he = face->facet_begin();
        Vector face_normal = fnormals[face];
        do {
            Point p = he->vertex()->point();
            teeth_facets.push_back((float)p.x());
            teeth_facets.push_back((float)p.y());
            teeth_facets.push_back((float)p.z());
            teeth_normals.push_back((float)face_normal.x());
            teeth_normals.push_back((float)face_normal.y());
            teeth_normals.push_back((float)face_normal.z());
        } while (++he != face->facet_begin());
    }
}

void Function::update_teeth_mesh(DataList& teeth_facets, DataList& teeth_normals, bool flag_lower)
{
    Polyhedron& teeth_poly = (flag_lower) ? m_lower_teeth.get_teeth_poly() : m_upper_teeth.get_teeth_poly();
    generate_mesh_visu(teeth_poly, teeth_facets, teeth_normals);
}

void Function::update_inside_function(DataList& inside_ray_facets,
    DataList& inside_ray_normals,
    DataList& inside_ray_colors,
    DataList& inside_func_colors,
    DataList& inside_facets,
    DataList& inside_normals,
    bool flag_lower)
{
    Polyhedron& teeth_poly = (flag_lower) ? m_lower_teeth.get_teeth_poly() : m_upper_teeth.get_teeth_poly();
    VertDoubleMap& smoothed_ray_map = (flag_lower) ? m_lower_teeth.get_smoothed_ray_map() : m_upper_teeth.get_smoothed_ray_map();
    VertDoubleMap& inside_func_map = (flag_lower) ? m_lower_teeth.get_vertex_inside_function_map() : m_upper_teeth.get_vertex_inside_function_map();

    // Find min and max value
    double fmin = 0., fmax = 1.;
    for (Vertex_handle vd : vertices(teeth_poly))
    {
        fmin = ortho_min(inside_func_map[vd], fmin);
        fmax = ortho_max(inside_func_map[vd], fmax);
    }

    spdlog::trace("Function::update_inside_function {0} teeth fmin:{1}, teeth fmax:{2}", flag_lower?"Lower":"Upper", fmin, fmax);

    //if (flag_lower)
    //    std::cout << "  Lower teeth fmin: " << fmin << ", lower teeth fmax: " << fmax << std::endl;
    //else
    //    std::cout << "  Upper teeth fmin: " << fmin << ", Upper teeth fmax: " << fmax << std::endl;

    // Fill Datalists
    FaceVecMap fnormals;
    CGAL::Polygon_mesh_processing::compute_face_normals(teeth_poly, boost::make_assoc_property_map(fnormals));
    inside_ray_facets.reserve(teeth_poly.size_of_facets() * 9);
    inside_ray_normals.reserve(teeth_poly.size_of_facets() * 9);
    inside_ray_colors.reserve(teeth_poly.size_of_facets() * 9);
    inside_func_colors.reserve(teeth_poly.size_of_facets() * 9);
    Ramp ramp;
    ramp.set_range(fmin, fmax);

    for (Facet_iterator face = teeth_poly.facets_begin(); face != teeth_poly.facets_end(); ++face)
    {
        Halfedge_facet_circulator he = face->facet_begin();
        Vector face_normal = fnormals[face];
        do {
            // smoothed ray field
            Point p = he->vertex()->point();
            inside_ray_facets.push_back((float)p.x());
            inside_ray_facets.push_back((float)p.y());
            inside_ray_facets.push_back((float)p.z());
            inside_ray_normals.push_back((float)face_normal.x());
            inside_ray_normals.push_back((float)face_normal.y());
            inside_ray_normals.push_back((float)face_normal.z());
            double r, g, b;
            ramp.get_color(smoothed_ray_map[he->vertex()], r, g, b);
            inside_ray_colors.push_back((float)r);
            inside_ray_colors.push_back((float)g);
            inside_ray_colors.push_back((float)b);
            ramp.get_color(inside_func_map[he->vertex()], r, g, b);
            inside_func_colors.push_back((float)r);
            inside_func_colors.push_back((float)g);
            inside_func_colors.push_back((float)b);
        } while (++he != face->facet_begin());
    }

    Polyhedron& inside_teeth_poly = (flag_lower) ? m_lower_teeth.get_inside_teeth_poly() : m_upper_teeth.get_inside_teeth_poly();

    generate_mesh_visu(inside_teeth_poly, inside_facets, inside_normals);
}

void Function::update_convex_function(DataList& local_convex_facets,
    DataList& local_convex_normals,
    DataList& local_convex_colors,
    DataList& convex_func_colors,
    bool flag_lower)
{
    Polyhedron& inside_teeth_poly = (flag_lower) ? m_lower_teeth.get_inside_teeth_poly() : m_upper_teeth.get_inside_teeth_poly();
    VertDoubleMap& smoothed_convex_map = (flag_lower) ? m_lower_teeth.get_smoothed_convex_map() : m_upper_teeth.get_smoothed_convex_map();
    VertDoubleMap& convex_func_map = (flag_lower) ? m_lower_teeth.get_vertex_convex_function_map() : m_upper_teeth.get_vertex_convex_function_map();

    // Find min and max value
    double fmin = 1e10, fmax = -1e10;
    for (Vertex_handle vd : vertices(inside_teeth_poly))
    {
        fmin = ortho_min(smoothed_convex_map[vd], fmin);
        fmin = ortho_min(convex_func_map[vd], fmin);
        fmax = ortho_max(smoothed_convex_map[vd], fmax);
        fmax = ortho_max(convex_func_map[vd], fmax);
    }

    spdlog::trace("Function::update_convex_function {0} teeth fmin:{1}, teeth fmax:{2}", flag_lower ? "Lower" : "Upper", fmin, fmax);

    //if (flag_lower)
    //    std::cout << "  Lower teeth fmin: " << fmin << ", lower teeth fmax: " << fmax << std::endl;
    //else
    //    std::cout << "  Upper teeth fmin: " << fmin << ", Upper teeth fmax: " << fmax << std::endl;

    // Fill Datalists
    FaceVecMap fnormals;
    CGAL::Polygon_mesh_processing::compute_face_normals(inside_teeth_poly, boost::make_assoc_property_map(fnormals));
    local_convex_facets.reserve(inside_teeth_poly.size_of_facets() * 9);
    local_convex_normals.reserve(inside_teeth_poly.size_of_facets() * 9);
    local_convex_colors.reserve(inside_teeth_poly.size_of_facets() * 9);
    convex_func_colors.reserve(inside_teeth_poly.size_of_facets() * 9);
    Ramp ramp;
    ramp.set_range(fmin, fmax);

    for (Facet_iterator face = inside_teeth_poly.facets_begin(); face != inside_teeth_poly.facets_end(); ++face)
    {
        Halfedge_facet_circulator he = face->facet_begin();
        Vector face_normal = fnormals[face];
        do {
            // smoothed ray field
            Point p = he->vertex()->point();
            local_convex_facets.push_back((float)p.x());
            local_convex_facets.push_back((float)p.y());
            local_convex_facets.push_back((float)p.z());
            local_convex_normals.push_back((float)face_normal.x());
            local_convex_normals.push_back((float)face_normal.y());
            local_convex_normals.push_back((float)face_normal.z());
            double r, g, b;
            ramp.get_color(smoothed_convex_map[he->vertex()], r, g, b);
            local_convex_colors.push_back((float)r);
            local_convex_colors.push_back((float)g);
            local_convex_colors.push_back((float)b);
            ramp.get_color(convex_func_map[he->vertex()], r, g, b);
            convex_func_colors.push_back((float)r);
            convex_func_colors.push_back((float)g);
            convex_func_colors.push_back((float)b);
        } while (++he != face->facet_begin());
    }
}

void Function::update_convexity_isovalue(DataList& convex_isoedges, double value, bool flag_lower)
{
    spdlog::trace("Function::update_convexity_isovalue value {0} flag_lower:{1}", value, flag_lower);

    Polyhedron& inside_teeth_poly = (flag_lower) ? m_lower_teeth.get_inside_teeth_poly() : m_upper_teeth.get_inside_teeth_poly();
    VertDoubleMap& convex_func_map = (flag_lower) ? m_lower_teeth.get_vertex_convex_function_map() : m_upper_teeth.get_vertex_convex_function_map();

    for (Facet_iterator face = inside_teeth_poly.facets_begin(); face != inside_teeth_poly.facets_end(); ++face)
    {
        Halfedge_facet_circulator he = face->facet_begin();
        PointList iso_points;
        do {
            double a = convex_func_map[he->vertex()] - value;
            double b = convex_func_map[he->opposite()->vertex()] - value;
            if (a * b < 0)
            {
                Point pa = he->vertex()->point();
                Point pb = he->opposite()->vertex()->point();
                Point iso;
                bool flag = find_level_set_point(pa, a, pb, b, 0., iso);
                if (flag)
                    iso_points.push_back(iso);
            }
        } while (++he != face->facet_begin());

        if (iso_points.size() == 2)
        {
            convex_isoedges.push_back((float)iso_points[0].x());
            convex_isoedges.push_back((float)iso_points[0].y());
            convex_isoedges.push_back((float)iso_points[0].z());
            convex_isoedges.push_back((float)iso_points[1].x());
            convex_isoedges.push_back((float)iso_points[1].y());
            convex_isoedges.push_back((float)iso_points[1].z());
        }
    }
}

void Function::update_segmentation(DataList& seg_facets, DataList& seg_normals)
{
    if (m_segment.size_of_vertices() == 0)
        return;

    generate_mesh_visu(m_segment, seg_facets, seg_normals);
}

void Function::update_validated_segmentation(int index, bool flag_lower, DataList& seg_facets, DataList& seg_normals)
{
    Polyhedron& tooth_poly = (flag_lower) ? m_lower_teeth.get_tooth_poly(index) : m_upper_teeth.get_tooth_poly(index);
    if (tooth_poly.empty())
        return;

    generate_mesh_visu(tooth_poly, seg_facets, seg_normals);
}

void Function::update_tooth_geodesic_function(int index, bool flag_lower, DataList& poly_colors)
{
    Polyhedron& tooth_poly = (flag_lower) ? m_lower_teeth.get_tooth_poly(index) : m_upper_teeth.get_tooth_poly(index);
    if (tooth_poly.empty())
        return;

    VertDoubleMap& tooth_func = (flag_lower) ? m_lower_teeth.get_geodesic_distance_map(index) : m_upper_teeth.get_geodesic_distance_map(index);
    if (tooth_func.empty())
        return;

    // Find min and max value
    double fmin = 1e10, fmax = -1e10;
    for (Vertex_handle vd : vertices(tooth_poly))
    {
        fmin = ortho_min(tooth_func[vd], fmin);
        fmax = ortho_max(tooth_func[vd], fmax);
    }
    spdlog::trace("Function::update_tooth_geodesic_function: eodesic distance fmin:{0}, fmax:{1}", 
        fmin, fmax);

    //std::cout << "  Geodesic distance fmin: " << fmin << ", geodesic distance fmax: " << fmax << std::endl;


    poly_colors.reserve(tooth_poly.size_of_facets() * 9);
    Ramp ramp;
    ramp.set_range(fmin, fmax);

    for (Facet_iterator face = tooth_poly.facets_begin(); face != tooth_poly.facets_end(); ++face)
    {
        Halfedge_facet_circulator he = face->facet_begin();
        do {
            double r, g, b;
            ramp.get_color(tooth_func[he->vertex()], r, g, b);
            poly_colors.push_back((float)r);
            poly_colors.push_back((float)g);
            poly_colors.push_back((float)b);
        } while (++he != face->facet_begin());
    }
}

void Function::update_pad_outline(int index, bool flag_lower, DataList& poly_facets, DataList& poly_normals)
{
    Polyhedron &pad_outline_poly = (flag_lower) ? m_lower_teeth.get_pad_outline_poly(index) : m_upper_teeth.get_pad_outline_poly(index);
    if (pad_outline_poly.empty())
        return;
    generate_mesh_visu(pad_outline_poly, poly_facets, poly_normals);
}

void Function::update_tooth_cut_plane(int index, bool flag_lower, DataList& upper_bound, DataList& upper_points, DataList& cut_plane_points, DataList& cut_plane_normals)
{
    HEdgeSet& upper_bound_set = (flag_lower) ? m_lower_teeth.get_upper_bound_set(index) : m_upper_teeth.get_upper_bound_set(index);
    PointList& cut_points = (flag_lower) ? m_lower_teeth.get_cut_points(index) : m_upper_teeth.get_cut_points(index);
    Plane& cut_plane = (flag_lower) ? m_lower_teeth.get_cut_plane(index) : m_upper_teeth.get_cut_plane(index);

    if (upper_bound_set.size() == 0 || cut_points.size() == 0)
        return;

    for (int i = 0; i < cut_points.size(); i++)
    {
        Point p = cut_points[i];
        upper_points.push_back((float)p.x());
        upper_points.push_back((float)p.y());
        upper_points.push_back((float)p.z());
    }

    for (auto& elem : upper_bound_set)
    {
        Point source = elem->vertex()->point();
        upper_bound.push_back((float)source.x());
        upper_bound.push_back((float)source.y());
        upper_bound.push_back((float)source.z());
        Point target = elem->opposite()->vertex()->point();
        upper_bound.push_back((float)target.x());
        upper_bound.push_back((float)target.y());
        upper_bound.push_back((float)target.z());
    }

    Vector normal;
    PointList corners;
    compute_plane_corners(cut_plane, cut_points[0], cut_points[1], corners, normal);

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

void Function::compute_plane_corners(Plane& plane, Point& center, Point& right, PointList& corners, Vector& normal)
{
    double scale = 3.;

    Point plane_center = plane.projection(center);
    Point plane_right = plane.projection(right);
    Vector horizontal = plane_right - plane_center;
    horizontal = horizontal / std::sqrt(horizontal.squared_length());
    normal = plane.orthogonal_vector();
    Vector vertical = CGAL::cross_product(horizontal, normal);
    vertical = vertical / std::sqrt(vertical.squared_length());

    Point left_down = plane_center - scale * horizontal - scale * vertical;
    Point left_up = plane_center - scale * horizontal + scale * vertical;
    Point right_up = plane_center + scale * horizontal + scale * vertical;
    Point right_down = plane_center + scale * horizontal - scale * vertical;

    corners.push_back(left_down);
    corners.push_back(left_up);
    corners.push_back(right_up);
    corners.push_back(left_down);
    corners.push_back(right_up);
    corners.push_back(right_down);
}

void Function::update_pair_wire_shape(int index, bool flag_lower,
    DataList& wire_shape_points,
    DataList& wire_shape_normals,
    DataList& wire_extrude_shape_points,
    DataList& wire_extrude_shape_normals)
{
    Wire& wire = (flag_lower) ? m_lower_teeth.get_wire(index) : m_upper_teeth.get_wire(index);
    Polyhedron& wire_shape = wire.get_wire_shape();

    generate_mesh_visu(wire_shape, wire_shape_points, wire_shape_normals);

    Polyhedron& wire_extrude_shape = wire.get_wire_extrude_shape();

    if (!wire_extrude_shape.empty())
    {
        generate_mesh_visu(wire_extrude_shape, wire_extrude_shape_points, wire_extrude_shape_normals);
    }
}

void Function::update_pad_shape(int index, bool flag_lower,
    DataList& pad_shape_points,
    DataList& pad_shape_normals)
{
    Tooth& tooth = (flag_lower) ? m_lower_teeth.get_tooth(index) : m_upper_teeth.get_tooth(index);
    bool flag_split = tooth.get_split_flag();

    Polyhedron pad_shape;
    if (flag_split)
        merge_polyehdral_surface(tooth.get_left_pad_shape_poly(), tooth.get_right_pad_shape_poly(), pad_shape);
    else
        copy_polyehdral_surface(tooth.get_pad_shape_poly(), pad_shape);

    if (!pad_shape.empty())
    {
        generate_mesh_visu(pad_shape, pad_shape_points, pad_shape_normals);
    }
}

void Function::update_teeth_distance_colors(bool flag_lower, double threshold, DataList& poly_colors)
{
    poly_colors.clear();

    VertDoubleMap& teeth_dist_map = (flag_lower) ? m_lower_teeth.get_teeth_distance_map() : m_upper_teeth.get_teeth_distance_map();
    Polyhedron& teeth_poly = (flag_lower) ? m_lower_teeth.get_teeth_poly() : m_upper_teeth.get_teeth_poly();

    double fmax = -1e10;
    for (Vertex_handle vd : vertices(teeth_poly))
    {
        fmax = ortho_max(teeth_dist_map[vd], fmax);
    }

    poly_colors.reserve(teeth_poly.size_of_facets() * 9);
    Ramp ramp;
    ramp.set_range(0., ortho_min(fmax, threshold));

    for (Facet_iterator face = teeth_poly.facets_begin(); face != teeth_poly.facets_end(); ++face)
    {
        Halfedge_facet_circulator he = face->facet_begin();
        do {
            double r, g, b;
            double value = ortho_min(teeth_dist_map[he->vertex()], threshold);
            ramp.get_color(value, r, g, b);
            poly_colors.push_back((float)r);
            poly_colors.push_back((float)g);
            poly_colors.push_back((float)b);
        } while (++he != face->facet_begin());
    }
}


void Function::close_mesh(const std::string &pad, bool shell_solid)
{
    Polyhedron _m_pad_shape; 
    mesh_import(pad, _m_pad_shape);
    hole_filling(_m_pad_shape, shell_solid);
    mesh_export(_m_pad_shape, pad+"_closed.stl");
}

void Function::Expand_mesh(const std::string &pad, double thickness)
{
    Polyhedron _m_pad_shape; 
    mesh_import(pad, _m_pad_shape);
    Polyhedron _m_pad_shape_expanded;
    compute_enlarge_msh(_m_pad_shape, _m_pad_shape_expanded, thickness);
    mesh_export(_m_pad_shape_expanded, pad+"_expand.stl");
}

} // namespace OrthoLab
