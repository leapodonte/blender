#include "tooth.h"

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <boost/property_map/property_map.hpp>
#include <CGAL/bounding_box.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/convex_hull_2.h>
#include <CGAL/Subdivision_method_3.h>
#include <functional>

#include "teeth.h"
#include "mesh_io.h"
#include "meshing.h"
#include "save.h"
#include "extrusion.h"
#include "project.h"
#include "smooth_mesh.h"
#include "base_helpers.h"
#include <CGAL/alpha_wrap_3.h>
#include "file_services.h"

namespace OrthoLab {

Tooth::Tooth(const std::string& id, const Teeth& teeth):
SceneObject(id, teeth),
m_pad(id, *this),
m_teeth(teeth)
{
}

Tooth::~Tooth()
{
    reset();
}

void Tooth::reset()
{
    m_name = "";
    m_tooth_poly.clear();
    m_vert_ind_map.clear();
    m_vert_geodesic_map.clear();
    m_vert_smooth_geodesic_map.clear();
    m_center = Point(0., 0., 0.);
    reset_cut_plane();
    reset_pad();
}

void Tooth::reset_cut_plane()
{
    m_upper_bound.clear();
    m_upper_points.clear();
    m_cut_points.clear();
    m_cut_plane = Plane(0., 0., 0., 1.);
}

void Tooth::reset_pad()
{
    m_pad_outlines.clear();
    m_pad_shape.clear();
    m_pad_left.clear();
    m_pad_right.clear();
    m_pad_split_flag = false;
    m_pad.reset();
}

std::string Tooth::get_base_path() const
{
    return m_teeth.get_base_path()+
           PATH_SEPERATION + "t_" + get_id();           
}


void Tooth::save_to_db()
{
	std::string tooth_key = get_base_path();
	std::string tooth_path = ensure_directory_exists(tooth_key);
	tooth_path += PATH_SEPERATION;
    get_project().get_store().store_vector3(
		tooth_key + "pad_extrusion_direction", m_pad_extrusion_direction);
    get_project().get_store().store_transformation(
		tooth_key + "m_transform_left", m_transform_left);
    get_project().get_store().store_transformation(
		tooth_key + "m_transform_right", m_transform_right);
	std::vector<double> transformation_status;
	transformation_status.push_back(m_transform_status);
    get_project().get_store().store_vector(
		tooth_key + "m_transform_status", transformation_status);

	if (!get_tooth_poly().empty())
		mesh_export(m_tooth_poly, tooth_path + "m_tooth_poly" + get_id() + ".stl");
	if (!get_pad_shape_poly().empty())
		mesh_export(m_pad_shape, tooth_path + "m_pad_shape" + get_id() + ".stl");
	if (!get_pad_outline_poly().empty())
		mesh_export(m_pad_outlines, tooth_path + "m_pad_outlines" + get_id() + ".stl");
    if (!get_left_pad_shape_poly().empty())
        mesh_export(m_pad_left, tooth_path + "m_pad_left" + get_id() + ".stl");
    if (!get_right_pad_shape_poly().empty())
        mesh_export(m_pad_right, tooth_path + "m_pad_right" + get_id() + ".stl");	
    if (!get_left_pad_outline_poly().empty())
        mesh_export(m_pad_outlines_left, tooth_path + "m_pad_outlines_left" + get_id() + ".stl");
    if (!get_right_pad_outline_poly().empty())
        mesh_export(m_pad_outlines_right, tooth_path + "m_pad_outlines_right" + get_id() + ".stl");

	if (!get_tooth_poly(true).empty())
		mesh_export(m_tooth_poly_t, tooth_path + "m_tooth_poly_t" + get_id() + ".stl");
	if (!get_pad_shape_poly(true).empty())
		mesh_export(m_pad_shape_t, tooth_path + "m_pad_shape_t" + get_id() + ".stl");
	if (!get_pad_outline_poly(true).empty())
		mesh_export(m_pad_outlines_t, tooth_path + "m_pad_outlines_t" + get_id() + ".stl");
	if (!get_left_pad_shape_poly(true).empty())
		mesh_export(m_pad_left_t, tooth_path + "m_pad_left_t" + get_id() + ".stl");
	if (!get_right_pad_shape_poly(true).empty())
		mesh_export(m_pad_right_t, tooth_path + "m_pad_right_t" + get_id() + ".stl");
    if (!get_left_pad_outline_poly(true).empty())
        mesh_export(m_pad_outlines_left_t, tooth_path + "m_pad_outlines_left_t" + get_id() + ".stl");
    if (!get_right_pad_outline_poly(true).empty())
        mesh_export(m_pad_outlines_right_t, tooth_path + "m_pad_outlines_right_t" + get_id() + ".stl");

 }

void Tooth::load_from_db()
{
    std::string tooth_key = get_base_path();
    std::string tooth_path = ensure_directory_exists(tooth_key);
    tooth_path += PATH_SEPERATION;

    get_project().get_store().retrieve_vector3(tooth_key + "pad_extrusion_direction", m_pad_extrusion_direction);
    std::cout << "load_from_db:: "  << tooth_key + "pad_extrusion_direction" << m_pad_extrusion_direction << std::endl;
    std::vector<double> transformation_status;
    int rt = get_project().get_store().retrieve_vector(tooth_key + "m_transform_status",transformation_status);
    if (transformation_status.size())
        m_transform_status = transformation_status[0];
    else
        m_transform_status = 1; // hack

    std::string shape_file_path = tooth_path + "m_tooth_poly" + get_id() + ".stl";
    if(file_exists(shape_file_path))
        mesh_import(shape_file_path, m_tooth_poly);

    shape_file_path = tooth_path + "m_pad_shape" + get_id() + ".stl";
    if(file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_shape);

    shape_file_path = tooth_path + "m_pad_outlines" + get_id() + ".stl";
    if(file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_outlines);

    shape_file_path = tooth_path + "m_pad_left" + get_id() + ".stl";
    if(file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_left);
        
    shape_file_path = tooth_path + "m_pad_right" + get_id() + ".stl";
    if(file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_right);

    shape_file_path = tooth_path + "m_pad_outlines_left" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_outlines_left);

    shape_file_path = tooth_path + "m_pad_outlines_right" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_outlines_right);

    // transformed
    shape_file_path = tooth_path + "m_tooth_poly_t" + get_id() + ".stl";
    if(file_exists(shape_file_path))
        mesh_import(shape_file_path, m_tooth_poly_t);

    shape_file_path = tooth_path + "m_pad_shape_t" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_shape_t);

    shape_file_path = tooth_path + "m_pad_outlines_t" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_outlines_t);

    shape_file_path = tooth_path + "m_pad_left_t" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_left_t);

    shape_file_path = tooth_path + "m_pad_right_t" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_right_t);

    shape_file_path = tooth_path + "m_pad_outlines_left_t" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_outlines_left_t);

    shape_file_path = tooth_path + "m_pad_outlines_right_t" + get_id() + ".stl";
    if (file_exists(shape_file_path))
        mesh_import(shape_file_path, m_pad_outlines_right_t);
}

/*void Tooth::set_tooth_mesh(Polyhedron& seg_poly, Point centroid)
{
    reset();

    m_center = centroid;
    copy_polyehdral_surface(seg_poly, m_tooth_poly);
    std::cout << "    Load tooth mesh with " << m_tooth_poly.size_of_vertices() << " vertices and " << m_tooth_poly.size_of_facets() << " faces." << std::endl;

    // hole filling
    HEdgeList border_cycles;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m_tooth_poly, std::back_inserter(border_cycles));

    if (border_cycles.size() > 1) // contain holes
    {
        DoubleList border_lengths;
        for (int i = 0; i < border_cycles.size(); i++)
        {
            double length = compute_length_for_hole(border_cycles[i], m_tooth_poly);
            border_lengths.push_back(length);
        }

        unsigned int nb_holes = 0;
        std::vector<int> indices(border_lengths.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool { return border_lengths[A] > border_lengths[B]; });

        for (int i = 1; i < indices.size(); i++)
        {
            FacetList  patch_facets;
            VertexList patch_vertices;
            bool success = std::get<0>(CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(m_tooth_poly,
                border_cycles[indices[i]],
                CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
                .vertex_output_iterator(std::back_inserter(patch_vertices))
                .vertex_point_map(get(CGAL::vertex_point, m_tooth_poly))
                .geom_traits(Kernel())));
            ++nb_holes;
        }

        std::cout << "    Filled " << nb_holes << "holes!" << std::endl;
    }

    initialize_vertex_indices();
}*/

void Tooth::set_tooth_mesh(Polyhedron& seg_poly, Polyhedron& teeth_poly)
{
    reset();

    m_center = compute_centroid(seg_poly);
    std::cout << "Recompute center: " << m_center << std::endl;

    // fit plane
    TriangleList candidate_faces;
    copy_faces_to_triangles(seg_poly, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    
    // project to 2d
    Point plane_center = fitting_plane.projection(m_center);
    Vector plane_base_1 = fitting_plane.base1();
    Vector plane_base_2 = fitting_plane.base2();
    plane_base_1 = plane_base_1 / std::sqrt(plane_base_1.squared_length());
    plane_base_2 = plane_base_2 / std::sqrt(plane_base_2.squared_length());

    // compute boundary
    HEdgeList border_cycles;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(seg_poly, std::back_inserter(border_cycles));

    double max_length = 0.;
    int max_ind = -1;

    for(int i = 0; i < border_cycles.size(); i++)
    {
        double length = compute_length_for_hole(border_cycles[i], seg_poly);
        if(length > max_length)
        {
            max_length = length;
            max_ind = i;
        }
    }

    Point2dList tooth_pts_2d;
    for (Halfedge_handle hc : CGAL::halfedges_around_face(border_cycles[max_ind], seg_poly))
    {
        Point point = hc->vertex()->point();
        tooth_pts_2d.push_back(to_2d(fitting_plane, plane_center, plane_base_1, plane_base_2, point));
    }

    // Convex hull
    Point2dList convex_tooth_pts_2d;
    CGAL::convex_hull_2(tooth_pts_2d.begin(), tooth_pts_2d.end(), std::back_inserter(convex_tooth_pts_2d));

    // Build segment tree
    int num_points = convex_tooth_pts_2d.size();
    SegmentList convex_tooth_semgents;
    for(int j = 0; j < num_points - 1; j++)
        convex_tooth_semgents.push_back(Segment(Point(convex_tooth_pts_2d[j].x()  , convex_tooth_pts_2d[j].y()  , 0.),
                                              Point(convex_tooth_pts_2d[j+1].x(), convex_tooth_pts_2d[j+1].y(), 0.)));

    convex_tooth_semgents.push_back(Segment(Point(convex_tooth_pts_2d[num_points-1].x(), convex_tooth_pts_2d[num_points-1].y(), 0.),
                                          Point(convex_tooth_pts_2d[0].x()           , convex_tooth_pts_2d[0].y()           , 0.)));

    SegmentTree tree(convex_tooth_semgents.begin(),convex_tooth_semgents.end());
    tree.accelerate_distance_queries();

    // Construct implicit function
    VertDoubleMap vert_dist_map;
    VertIntMap vert_ind_map;
    int index = 0;
    for(Vertex_handle vd: vertices(teeth_poly))
    {
        Point_2 p2d = to_2d(fitting_plane, plane_center, plane_base_1, plane_base_2, vd->point());
        Point p3d(p2d.x(), p2d.y(), 0.);
        double dist = std::sqrt(tree.squared_distance(p3d));
        Point proj = tree.closest_point(p3d);
        if((p3d-CGAL::ORIGIN).squared_length() > (proj-CGAL::ORIGIN).squared_length())
            dist = -dist;
        vert_dist_map.insert({vd, dist - 0.05});
        vert_ind_map.insert({vd, index});
        index++;
    }

    // extract isofaces
    PointList poly_points;
    extract_inside_vertices_from_function(teeth_poly, poly_points);

    HEdgeIntMap he_point_map;
    int count = extract_isovertices_from_function(teeth_poly, vert_dist_map, 0., poly_points, he_point_map);

    FacetList candidate_facets;
    for (Facet_handle fd : faces(teeth_poly))
    {
        Halfedge_facet_circulator he = fd->facet_begin();
        Vertex_handle va = he->vertex();
        Vertex_handle vb = he->next()->vertex();
        Vertex_handle vc = he->next()->next()->vertex();

        if(CGAL::squared_distance(va->point(), m_center) < 30. && CGAL::squared_distance(vb->point(), m_center) < 30. && CGAL::squared_distance(vc->point(), m_center) < 30.)
            candidate_facets.push_back(fd);
    }

    PolyList  poly_faces;
    extract_positive_faces_from_function(vert_dist_map, 0., vert_ind_map, he_point_map, candidate_facets, poly_faces);

    m_tooth_poly.clear();
    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, m_tooth_poly);
    CGAL::Polygon_mesh_processing::keep_largest_connected_components(m_tooth_poly, 1);

    m_center = compute_centroid(m_tooth_poly);

    std::cout << "    Load tooth mesh with " << m_tooth_poly.size_of_vertices() << " vertices and " << m_tooth_poly.size_of_facets() << " faces." << std::endl;

    initialize_vertex_indices();
}

void Tooth::set_name(std::string name) { m_name = name; }

Point Tooth::compute_centroid(Polyhedron& poly)
{
    double xsum = 0., ysum = 0., zsum = 0.;
    int sum = 0;

    for(Vertex_handle vd: vertices(poly))
    {
        Point p = vd->point();
        xsum += p.x();
        ysum += p.y();
        zsum += p.z();
        sum += 1;
    }

    Point centroid(xsum / sum, ysum / sum, zsum / sum);

    return centroid;
}

double Tooth::compute_length_for_hole(Halfedge_handle he, Polyhedron& mesh)
{
    double length = 0.;

    for (Halfedge_handle hc : CGAL::halfedges_around_face(he, mesh))
        length += std::sqrt(CGAL::squared_distance(hc->vertex()->point(), hc->opposite()->vertex()->point()));

    return length;
}

void Tooth::initialize_vertex_indices()
{
    m_vert_ind_map.clear();
    int index = 0;

    for (Vertex_handle vd : vertices(m_tooth_poly))
    {
        m_vert_ind_map.insert({ vd, index });
        index++;
    }
}

std::string Tooth::get_name() { return m_name; }
bool Tooth::get_init_flag() { return (m_tooth_poly.size_of_vertices() > 0); }
const Point& Tooth::get_center() { return m_center; }
VertDoubleMap& Tooth::get_geodesic_distance_map() { return m_vert_geodesic_map; }

void Tooth::set_left_transformation(const Transformation t)
{
    m_transform_left = t;
    m_transform_status = 1;
    m_tooth_poly_t.clear();
    m_pad_outlines_t.clear();
    m_pad_shape_t.clear();
    m_pad_left_t.clear();
    m_pad_right_t.clear();
}
Polyhedron &Tooth::get_tooth_poly(bool transformed)
{
    if (!transformed)
        return m_tooth_poly;
    if (m_tooth_poly_t.empty() && m_transform_status == 1)
    {
        copy_polyehdral_surface(m_tooth_poly, m_tooth_poly_t);
        transform_mesh(m_tooth_poly_t, m_transform_left); 
    }
    return m_tooth_poly_t;
}

Polyhedron &Tooth::get_pad_outline_poly(bool transformed)
{
    if (!transformed)
        return m_pad_outlines;
    if (m_pad_outlines_t.empty())
    {
        copy_polyehdral_surface(m_pad_outlines, m_pad_outlines_t);
        transform_mesh(m_pad_outlines_t, m_transform_left);
    }
    return m_pad_outlines_t;
}

Polyhedron &Tooth::get_pad_shape_poly(bool transformed)
{
    if (!transformed)
        return m_pad_shape;
    if (m_pad_shape_t.empty())
    {
        copy_polyehdral_surface(m_pad_shape, m_pad_shape_t);
        transform_mesh(m_pad_shape_t, m_transform_left);
    }
    return m_pad_shape_t;
}

Polyhedron &Tooth::get_left_pad_shape_poly(bool transformed)
{
    if (!transformed)
        return m_pad_left;
    if (m_pad_left_t.empty())
    {
        copy_polyehdral_surface(m_pad_left,m_pad_left_t);
        transform_mesh(m_pad_left_t, m_transform_left);
    }
    return m_pad_left_t;    
}

Polyhedron &Tooth::get_right_pad_shape_poly(bool transformed)
{
    if (!transformed)
        return m_pad_right;
    if (m_pad_right_t.empty())
    {
        copy_polyehdral_surface(m_pad_right, m_pad_right_t);
        transform_mesh(m_pad_right_t, m_transform_right);
    }
    return m_pad_right_t;
}

Polyhedron& Tooth::get_left_pad_outline_poly(bool transformed)
{
    if (!transformed)
        return m_pad_outlines_left;
    if (m_pad_outlines_left_t.empty())
    {
        copy_polyehdral_surface(m_pad_outlines_left, m_pad_outlines_left_t);
        transform_mesh(m_pad_outlines_left_t, m_transform_left);
    }
    return m_pad_outlines_left_t;
}

Polyhedron& Tooth::get_right_pad_outline_poly(bool transformed)
{
	if (!transformed)
		return m_pad_outlines_right;
	if (m_pad_outlines_right_t.empty())
	{
		copy_polyehdral_surface(m_pad_outlines_right, m_pad_outlines_right_t);
		transform_mesh(m_pad_outlines_right_t, m_transform_right);
	}
	return m_pad_outlines_right_t;
}

HEdgeSet& Tooth::get_upper_bound_set() { return m_upper_bound; }
PointList& Tooth::get_cut_points() { return m_cut_points; }
Plane& Tooth::get_cut_plane() { return m_cut_plane; }
bool Tooth::get_split_flag() { return m_pad_split_flag; }

Polyhedron& Tooth::get_bonding_surface() { return m_bonding_surface; }

bool Tooth::compute_geodedic_pad_outlines(int discretize_step, double isovalue, Polyhedron& teeth_poly, bool flag_smooth)
{
    if (m_tooth_poly.size_of_vertices() == 0)
    {
        std::cout << "    No available tooth polyhedral surface!" << std::endl;
        return false;
    }

    compute_geodesic_distance(discretize_step);
    update_geodesic_pad_outlines(isovalue, teeth_poly, flag_smooth);

    return true;
}

bool Tooth::compute_circular_pad_outlines(int discretize_step, double isovalue, double min_area)
{
    if (m_tooth_poly.size_of_vertices() == 0)
    {
        std::cout << "    No available tooth polyhedral surface!" << std::endl;
        return false;
    }

    compute_geodesic_distance(discretize_step);

    double radius = std::sqrt(min_area / 3.141592653);
    update_circular_pad_outlines(isovalue, radius);

    return true;
}

void Tooth::compute_geodesic_distance(int discretize_step)
{
    // clear map
    m_vert_geodesic_map.clear();

    Geodesic_tree geodesic_tree(m_tooth_poly,
        get(boost::vertex_external_index, m_tooth_poly),
        get(CGAL::halfedge_external_index, m_tooth_poly),
        get(CGAL::face_external_index, m_tooth_poly),
        get(CGAL::vertex_point, m_tooth_poly));

    // discretize tooth outline
    double offset = 1. / discretize_step;
    std::vector<Face_location> faceLocations;

    for (Halfedge_iterator hedge = m_tooth_poly.halfedges_begin(); hedge != m_tooth_poly.halfedges_end(); ++hedge)
    {
        if (hedge->is_border()) // infinite face
        {
            Facet_handle face = hedge->opposite()->facet(); // opposite face
            Vertex_handle va = hedge->vertex();
            Vertex_handle vb = hedge->opposite()->vertex();
            Halfedge_handle face_hedge = face->halfedge()->prev();
            int index_a = -1, index_b = -1;

            for (int i = 0; i < 3; i++)
            {
                if (face_hedge->vertex() == va)
                    index_a = i;
                if (face_hedge->vertex() == vb)
                    index_b = i;
                face_hedge = face_hedge->next();
            }

            for (int i = 0; i < discretize_step; i++)  // discretize edge
            {
                DoubleList weights(3, 0);
                weights[index_a] = offset * i;
                weights[index_b] = 1 - weights[index_a];

                Geodesic_coordinates face_location = { {weights[0], weights[1], weights[2]} };
                faceLocations.push_back(Face_location(face, face_location));
            }
        }
    }

    geodesic_tree.add_source_points(faceLocations.begin(), faceLocations.end());

    for (Vertex_handle vd : vertices(m_tooth_poly))
    {
        double dist = geodesic_tree.shortest_distance_to_source_points(vd).first;
        m_vert_geodesic_map.insert({ vd, dist });
    }
}

bool Tooth::update_geodesic_pad_outlines(double isovalue, Polyhedron& teeth_poly, bool flag_smooth)
{
    if (m_vert_geodesic_map.empty())
    {
        std::cout << "    No available geodesic distance map!" << std::endl;
        return false;
    }

    Teeth::compute_smoothed_field(pad_geodesic_smooth_range, m_vert_geodesic_map, m_vert_smooth_geodesic_map);

    reset_pad();

    PointList poly_points;
    extract_inside_vertices_from_function(m_tooth_poly, poly_points);

    HEdgeIntMap he_point_map;
    int count = extract_isovertices_from_function(m_tooth_poly, m_vert_smooth_geodesic_map, isovalue, poly_points, he_point_map);

    FacetList candidate_facets;
    for (Facet_handle fd : faces(m_tooth_poly))
        candidate_facets.push_back(fd);

    PolyList  poly_faces;
    extract_positive_faces_from_function(m_vert_smooth_geodesic_map, isovalue, m_vert_ind_map, he_point_map, candidate_facets, poly_faces);

    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, m_pad_outlines);
    CGAL::Polygon_mesh_processing::keep_largest_connected_components(m_pad_outlines, 1);

    // smooth pad boundary
    if(flag_smooth)
        smooth_geodesic_pad_outlines(teeth_poly);

    double area = CGAL::Polygon_mesh_processing::area(m_pad_outlines);
    std::cout << "    Compute pad mesh with " << m_pad_outlines.size_of_vertices() << " vertices and " << m_pad_outlines.size_of_facets() << " faces with area " << area << "." << std::endl;
    return true;
}

void Tooth::smooth_geodesic_pad_outlines(Polyhedron& teeth_poly)
{
    if(m_pad_outlines.size_of_vertices() == 0)
    {
        std::cout << "    No available pad outline to smooth!" << std::endl;
        return;
    }

    // fit plane
    TriangleList candidate_faces;
    copy_faces_to_triangles(m_tooth_poly, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());

    // project to 2d
    Point plane_center = fitting_plane.projection(m_center);
    Vector plane_base_1 = fitting_plane.base1();
    Vector plane_base_2 = fitting_plane.base2();
    plane_base_1 = plane_base_1 / std::sqrt(plane_base_1.squared_length());
    plane_base_2 = plane_base_2 / std::sqrt(plane_base_2.squared_length());

    Point2dList pad_outlines_2d;
    for(Vertex_handle vd: vertices(m_pad_outlines))
    {
        Halfedge_vertex_circulator he = vd->vertex_begin();
        do {
            if(he->is_border_edge())
            {
                Point point = vd->point();
                pad_outlines_2d.push_back(to_2d(fitting_plane, plane_center, plane_base_1, plane_base_2, point));
                break;
            }
        } while (++he != vd->vertex_begin());
    }
    std::cout << "      Found " << pad_outlines_2d.size() << " outline points in 2d." << std::endl;

    // Convex hull
    Point2dList convex_pad_outlines_2d;
    CGAL::convex_hull_2(pad_outlines_2d.begin(), pad_outlines_2d.end(), std::back_inserter(convex_pad_outlines_2d));
    std::cout << "      Found " << pad_outlines_2d.size() << " convex hull points in 2d." << std::endl;

    // Smooth
    int num_points = convex_pad_outlines_2d.size();
    /*for(int i = 0; i < 1; i++)
    {
        Point2dList new_pad_outlines_2d;
        new_pad_outlines_2d.reserve(num_points);
        new_pad_outlines_2d.push_back(CGAL::midpoint(convex_pad_outlines_2d[1], convex_pad_outlines_2d[num_points-1]));

        for(int j = 1; j < num_points - 1; j++)
        {
            Point_2 mid = CGAL::midpoint(convex_pad_outlines_2d[j-1], convex_pad_outlines_2d[j+1]);
            new_pad_outlines_2d.push_back(mid);
        }
        new_pad_outlines_2d.push_back(CGAL::midpoint(convex_pad_outlines_2d[0], convex_pad_outlines_2d[num_points-2]));

        convex_pad_outlines_2d.swap(new_pad_outlines_2d);
        new_pad_outlines_2d.clear();
    }*/

    // build segment tree
    SegmentList convex_pad_semgents;
    for(int j = 0; j < num_points - 1; j++)
        convex_pad_semgents.push_back(Segment(Point(convex_pad_outlines_2d[j].x()  , convex_pad_outlines_2d[j].y()  , 0.),
                                              Point(convex_pad_outlines_2d[j+1].x(), convex_pad_outlines_2d[j+1].y(), 0.)));

    convex_pad_semgents.push_back(Segment(Point(convex_pad_outlines_2d[num_points-1].x(), convex_pad_outlines_2d[num_points-1].y(), 0.),
                                          Point(convex_pad_outlines_2d[0].x()           , convex_pad_outlines_2d[0].y()           , 0.)));

    SegmentTree tree(convex_pad_semgents.begin(),convex_pad_semgents.end());
    tree.accelerate_distance_queries();

    // construct implicit function
    VertDoubleMap vert_dist_map;
    VertIntMap vert_ind_map;
    int index = 0;
    for(Vertex_handle vd: vertices(teeth_poly))
    {
        Point_2 p2d = to_2d(fitting_plane, plane_center, plane_base_1, plane_base_2, vd->point());
        Point p3d(p2d.x(), p2d.y(), 0.);
        double dist = std::sqrt(tree.squared_distance(p3d));
        Point proj = tree.closest_point(p3d);
        if((p3d-CGAL::ORIGIN).squared_length() > (proj-CGAL::ORIGIN).squared_length())
            dist = -dist;
        vert_dist_map.insert({vd, dist});
        vert_ind_map.insert({vd, index});
        index++;
    }

    // extract isofaces
    PointList poly_points;
    extract_inside_vertices_from_function(teeth_poly, poly_points);

    HEdgeIntMap he_point_map;
    int count = extract_isovertices_from_function(teeth_poly, vert_dist_map, 0., poly_points, he_point_map);

    FacetList candidate_facets;
    for (Facet_handle fd : faces(teeth_poly))
    {
        Halfedge_facet_circulator he = fd->facet_begin();
        Vertex_handle va = he->vertex();
        Vertex_handle vb = he->next()->vertex();
        Vertex_handle vc = he->next()->next()->vertex();

        if(CGAL::squared_distance(va->point(), m_center) < 50. && CGAL::squared_distance(vb->point(), m_center) < 50. && CGAL::squared_distance(vc->point(), m_center) < 50.)
            candidate_facets.push_back(fd);
    }
        

    PolyList  poly_faces;
    extract_positive_faces_from_function(vert_dist_map, 0., vert_ind_map, he_point_map, candidate_facets, poly_faces);

    m_pad_outlines.clear();
    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, m_pad_outlines);
    CGAL::Polygon_mesh_processing::keep_largest_connected_components(m_pad_outlines, 1);

    // hole filling
    HEdgeList border_cycles;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m_pad_outlines, std::back_inserter(border_cycles));

    if (border_cycles.size() > 1) // contain holes
    {
        DoubleList border_lengths;
        for (int i = 0; i < border_cycles.size(); i++)
        {
            double length = compute_length_for_hole(border_cycles[i], m_pad_outlines);
            border_lengths.push_back(length);
        }

        unsigned int nb_holes = 0;
        std::vector<int> indices(border_lengths.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool { return border_lengths[A] > border_lengths[B]; });

        for (int i = 1; i < indices.size(); i++)
        {
            FacetList  patch_facets;
            VertexList patch_vertices;
            bool success = std::get<0>(CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(m_pad_outlines,
                border_cycles[indices[i]],
                CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
                .vertex_output_iterator(std::back_inserter(patch_vertices))
                .vertex_point_map(get(CGAL::vertex_point, m_pad_outlines))
                .geom_traits(Kernel())));
            ++nb_holes;
        }

        std::cout << "    Filled " << nb_holes << "holes!" << std::endl;
    }


    std::cout << "Finish!" << std::endl;
}

bool Tooth::update_circular_pad_outlines(double isovalue, double radius)
{
    if (m_vert_geodesic_map.empty())
    {
        std::cout << "    No available geodesic distance map!" << std::endl;
        return false;
    }

    Teeth::compute_smoothed_field(pad_geodesic_smooth_range, m_vert_geodesic_map, m_vert_smooth_geodesic_map);

    Point geodesic_center = m_center;
    //compute_geodesic_center(geodesic_center);

    Geodesic_tree geodesic_tree(m_tooth_poly,
        get(boost::vertex_external_index, m_tooth_poly),
        get(CGAL::halfedge_external_index, m_tooth_poly),
        get(CGAL::face_external_index, m_tooth_poly),
        get(CGAL::vertex_point, m_tooth_poly));

    Face_location center_location = geodesic_tree.locate<AABBTraits>(geodesic_center);
    geodesic_tree.add_source_point(center_location);

    VertDoubleMap vert_value_map;
    for (Vertex_handle vd : vertices(m_tooth_poly))
    {
        double center_dist = geodesic_tree.shortest_distance_to_source_points(vd).first;
        double bound_dist = m_vert_smooth_geodesic_map[vd];
        double value = std::min(bound_dist - isovalue, radius - center_dist);
        vert_value_map.insert({vd, value});
    }

    reset_pad();

    PointList poly_points;
    extract_inside_vertices_from_function(m_tooth_poly, poly_points);

    HEdgeIntMap he_point_map;
    int count = extract_isovertices_from_function(m_tooth_poly, vert_value_map, 0., poly_points, he_point_map);

    FacetList candidate_facets;
    for (Facet_handle fd : faces(m_tooth_poly))
        candidate_facets.push_back(fd);

    PolyList  poly_faces;
    extract_positive_faces_from_function(vert_value_map, 0., m_vert_ind_map, he_point_map, candidate_facets, poly_faces);

    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, m_pad_outlines);
    CGAL::Polygon_mesh_processing::keep_largest_connected_components(m_pad_outlines, 1);
    double area = CGAL::Polygon_mesh_processing::area(m_pad_outlines);

    std::cout << "    Compute pad mesh with " << m_pad_outlines.size_of_vertices() << " vertices and " << m_pad_outlines.size_of_facets() << " faces with area " << area << "." << std::endl;
    return true;
}

bool Tooth::compute_tooth_cut_plane(KDTree& teeth_boundary_tree, AABBTree& teeth_tree)
{
    std::cout << "compute_tooth_cut_plane " << std::endl;
        
    if (m_vert_smooth_geodesic_map.empty() || m_pad_outlines.empty())
    {
        std::cout << "    No available geodesic distance map or pad outlines!" << std::endl;
        return false;
    }

    reset_cut_plane();

    // Find three points to fix the cut plane
    //Point geodesic_center;
    //compute_geodesic_center(geodesic_center);
    //m_cut_points.push_back(geodesic_center);
    PointList pad_verts;
    extract_inside_vertices_from_function(m_pad_outlines, pad_verts);
    Point pad_center = CGAL::centroid(pad_verts.begin(), pad_verts.end());
    m_cut_points.push_back(pad_center);


    SegmentList upper_bound_segs;
    compute_boundary_halfedges(teeth_boundary_tree, upper_bound_segs);
    Line fitting_line;
    CGAL::linear_least_squares_fitting_3(upper_bound_segs.begin(), upper_bound_segs.end(), fitting_line, CGAL::Dimension_tag<1>());
    Point upper_bound_center;
    compute_boundary_midpoint(fitting_line, upper_bound_center);
    //m_cut_points.push_back(upper_bound_center);

    //Point mid = CGAL::midpoint(m_cut_points[0], m_cut_points[1]);
    //Point mid_proj = teeth_tree.closest_point(mid);

    Point pzmin, pzmax;
    double zmin = 1e10, zmax = -1e10;
    for(Vertex_handle vd: vertices(m_pad_outlines))
    {
        Point p = vd->point();
        if(p.z() > zmax)
        {
            zmax = p.z();
            pzmax = p;
        }
        if(p.z() < zmin)
        {
            zmin = p.z();
            pzmin = p;
        }
    }
    if(CGAL::squared_distance(upper_bound_center, pzmin) > CGAL::squared_distance(upper_bound_center, pzmax))
        m_cut_points.push_back(pzmin);
    else
        m_cut_points.push_back(pzmax);

    Point mid = CGAL::midpoint(m_cut_points[0], m_cut_points[1]);
    Point mid_proj = teeth_tree.closest_point(mid);
    m_cut_points.push_back(mid_proj);

    // Fit the plane
    m_cut_plane = Plane(m_cut_points[0], m_cut_points[1], m_cut_points[2]);
    
    std::cout << "compute_tooth_cut_plane finished" << std::endl;
     
    return true;
}

void Tooth::compute_geodesic_center(Point& geodesic_center)
{
    double value = -1e10;

    for (Vertex_handle vd : vertices(m_tooth_poly))
    {
        double dist = m_vert_smooth_geodesic_map[vd];
        if (dist > value)
        {
            value = dist;
            geodesic_center = vd->point();
        }
    }
}

// Function to calculate the Euclidean distance between two points
double calculate_distance(const Point& p1, const Point& p2) 
{
    return std::sqrt(CGAL::squared_distance(p1, p2));
}

// Function to calculate the distance from a point to the center of a face
double distance_from_point_to_face(const Point& p, const Polyhedron::Facet& f) 
{
    auto h = f.facet_begin();
    int count = 0;

    double total_dist = 0;
    // Iterate around the face
    do {
        total_dist += calculate_distance(h->vertex()->point(), p);
        count++;
    } while (++h != f.facet_begin());

    // Compute the average
    Kernel::FT scale = 1.0 / count; 
    return total_dist * scale;
}

double Compute_area(Polyhedron::Facet_handle f)
{
    return Kernel::Compute_area_3()(
        (*f).halfedge()->vertex()->point(),
        (*f).halfedge()->next()->vertex()->point(),
        (*f).halfedge()->opposite()->vertex()->point());
}

void Tooth::compute_minimum_surface_pad(const double surface)
{
    std::cout << "compute_minimum_surface_pad" << std::endl;
    Point geodesic_center;
    compute_geodesic_center(geodesic_center);
    //
    if (m_pad_outlines.empty())
        return;

    AABBTree pad_tree;
    // initialize aabbtree
    pad_tree.insert(faces(m_pad_outlines).first, faces(m_pad_outlines).second, m_pad_outlines);
    pad_tree.accelerate_distance_queries();

    // Vector to hold faces and distances
    std::vector<std::pair<Polyhedron::Facet_handle, double>> face_distances;

    // Calculate distances and fill the vector
    for(auto fi = m_pad_outlines.facets_begin(); fi != m_pad_outlines.facets_end(); ++fi) {
        double dist = distance_from_point_to_face(geodesic_center, *fi);
        face_distances.push_back(std::make_pair(fi, dist));
    }

    // Sort based on distances
    std::sort(face_distances.begin(), face_distances.end(),
              [](const std::pair<Polyhedron::Facet_handle, double> &a,
                 const std::pair<Polyhedron::Facet_handle, double> &b) -> bool
              {
                  return a.second < b.second;
              });
    // Now face_distances is sorted based on the distance from the given point
    double area = 0;
    std::vector<Polyhedron::Facet_handle> extracted_facets;
    VertIntMap vert_value_map;
    PointList extracted_points;
    PolyList extracted_faces;
    int vertext_index = 0;
    for (auto fi2 = face_distances.begin(); fi2 != face_distances.end(); ++fi2)
    {
        auto one_facet = fi2->first;
        area += Compute_area(one_facet);
        IntList face;
        auto he = one_facet->facet_begin();
        do
        {
            int face_vertext_index = 0;
            auto map_index = vert_value_map.find(he->vertex());
            if (map_index == vert_value_map.end())
            {
                extracted_points.push_back(he->vertex()->point());
                face_vertext_index = vert_value_map[he->vertex()] = vertext_index++;
            }
            else
            {
                face_vertext_index = map_index->second;
            }
            face.push_back(face_vertext_index);
        } while (++he != one_facet->facet_begin());
        extracted_faces.push_back(face);

        if (area >= surface)
            break;
    }
    std::cout << "found vertex:" << vertext_index << std::endl;
    Polyhedron extracted_poly;
    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(extracted_points, extracted_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(extracted_points, extracted_faces, extracted_poly);
}

void Tooth::compute_boundary_halfedges(KDTree& boundary_tree, SegmentList& upper_bound_segs)
{
    for (Halfedge_iterator hedge = m_tooth_poly.halfedges_begin(); hedge != m_tooth_poly.halfedges_end(); ++hedge)
    {
        if (hedge->is_border()) // infinite face
        {
            PointList neighbor_points;
            Fuzzy_circle target_query(hedge->vertex()->point(), 1e-1);
            boundary_tree.search(std::back_inserter(neighbor_points), target_query);
            bool flag_target = (neighbor_points.size() > 0);

            neighbor_points.clear();
            Fuzzy_circle source_query(hedge->opposite()->vertex()->point(), 1e-1);
            boundary_tree.search(std::back_inserter(neighbor_points), source_query);
            bool flag_source = (neighbor_points.size() > 0);

            if (flag_target && flag_source)
            {
                m_upper_bound.insert(hedge);
                upper_bound_segs.push_back(Segment(hedge->vertex()->point(), hedge->opposite()->vertex()->point()));
                m_upper_points.insert(hedge->vertex());
                m_upper_points.insert(hedge->opposite()->vertex());
            }
        }
    }
    std::cout << "compute_boundary_halfedges m_upper_points.size:" << m_upper_points.size() << std::endl;
}

void Tooth::compute_boundary_midpoint(Line& fitting_line, Point& upper_bound_center)
{
    Point center = fitting_line.point(0);
    Vector direction = fitting_line.to_vector();
    direction = direction / std::sqrt(direction.squared_length());
    double min_scale = 1e10, max_scale = -1e10;

    for (auto& elem : m_upper_points)
    {
        Point proj = fitting_line.projection(elem->point());
        double value = std::sqrt(CGAL::squared_distance(proj, center));
        if ((proj - center) * direction < 0)
            value = -value;
        min_scale = ortho_min(min_scale, value);
        max_scale = ortho_max(max_scale, value);
    }

    // Cut with plane
    double avg_scale = (min_scale + max_scale) / 2.;
    upper_bound_center = center + avg_scale * direction;
}


// Vector compute_average_normal(const Polyhedron& polyhedron) {
//     Vector sum(0, 0, 0); // Initialize sum of normals
//     for (Polyhedron::Facet_const_iterator fit = polyhedron.facets_begin(); fit != polyhedron.facets_end(); ++fit) {
//         Polyhedron::Halfedge_around_facet_const_circulator hfc = fit->facet_begin();
//         // Assuming triangular facets for simplicity
//         Point a = hfc->vertex()->point();
//         Point b = (++hfc)->vertex()->point();
//         Point c = (++hfc)->vertex()->point();
        
//         // Compute the normal of the triangle
//         Vector normal = CGAL::cross_product(b - a, c - a);
//         normal = normal / std::sqrt(normal * normal); // Normalize
//         sum = sum + normal; // Add to the sum
//     }
    
//     if (sum.squared_length() > 0) {
//         sum = sum / std::sqrt(sum * sum); // Normalize the sum
//     }
//     return sum;
// }

void construct_cone(Polyhedron &cones_collection, Point base, Vector base_direction, Vector base_normal, Vector cone_direction, 
                    double radius, double height, int num_segments, Polyhedron& glue_zone_pad)
{ 
    Polyhedron cone;
    construct_cone(cone, base, base_direction, base_normal, cone_direction, radius, height, num_segments);
    // mesh_export(cones_collection, "cones_collection.stl");
    // mesh_export(cone, "cones_collection_coned_to_add.stl");
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(cone, cones_collection, cones_collection);
    
}

void _compute_bonding_surface(Polyhedron& m_pad_shape, Polyhedron& m_pad_outlines, Polyhedron& output_bonding_surface)
{
    // CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(m_pad_shape, m_pad_outlines, m_bonding_surface);
    // import the pad    
    // CGAL::Polygon_mesh_processing::IO::read_polygon_mesh("dataoutput_pad_outline_smooth/lowerpad3.stl", m_pad_shape);
    // CGAL::Polygon_mesh_processing::IO::read_polygon_mesh("dataoutput_pad_outline_smooth/lowerpad_outline3.stl", m_pad_outlines);

    std::cout << "m_pad_outlines facets: " << m_pad_outlines.size_of_facets() << std::endl;
    std::cout << "m_pad_shape facets: " << m_pad_shape.size_of_facets() << std::endl;
    // import the pad outline
    //


    // 1. find the glue fill zone
    Polyhedron glue_zone;
    double isovalue = get_variable_float(BONDING_ZONE_MARGINE, 0.25);
    int discretize_step = get_variable_float(BONDING_DISCRETIZE_STEP, 10);
    Polyhedron  m_tooth_poly = m_pad_outlines; // tricky
    Geodesic_tree geodesic_tree(m_tooth_poly,
                                get(boost::vertex_external_index, m_tooth_poly),
                                get(CGAL::halfedge_external_index, m_tooth_poly),
                                get(CGAL::face_external_index, m_tooth_poly),
                                get(CGAL::vertex_point, m_tooth_poly));

    // ----------------------compute_geodesic_distance---------------------------
    VertDoubleMap m_vert_geodesic_map;
    // discretize tooth outline
    double offset = 1. / discretize_step;
    std::vector<Face_location> faceLocations;
    for (Halfedge_iterator hedge = m_tooth_poly.halfedges_begin(); hedge != m_tooth_poly.halfedges_end(); ++hedge)
    {
        if (hedge->is_border()) // infinite face
        {
            Facet_handle face = hedge->opposite()->facet(); // opposite face
            Vertex_handle va = hedge->vertex();
            Vertex_handle vb = hedge->opposite()->vertex();
            Halfedge_handle face_hedge = face->halfedge()->prev();
            int index_a = -1, index_b = -1;

            for (int i = 0; i < 3; i++)
            {
                if (face_hedge->vertex() == va)
                    index_a = i;
                if (face_hedge->vertex() == vb)
                    index_b = i;
                face_hedge = face_hedge->next();
            }

            for (int i = 0; i < discretize_step; i++)  // discretize edge
            {
                DoubleList weights(3, 0);
                weights[index_a] = offset * i;
                weights[index_b] = 1 - weights[index_a];

                Geodesic_coordinates face_location = { {weights[0], weights[1], weights[2]} };
                faceLocations.push_back(Face_location(face, face_location));
            }
        }
    }

    geodesic_tree.add_source_points(faceLocations.begin(), faceLocations.end());

    for (Vertex_handle vd : vertices(m_tooth_poly))
    {
        double dist = geodesic_tree.shortest_distance_to_source_points(vd).first;
        m_vert_geodesic_map.insert({ vd, dist });
    }

    // ---------------update_geodesic_pad_outlines----------------------------------------
    int pad_geodesic_smooth_range = 0;
    VertDoubleMap m_vert_smooth_geodesic_map;
    Teeth::compute_smoothed_field(pad_geodesic_smooth_range, m_vert_geodesic_map,
                                  m_vert_smooth_geodesic_map);

    PointList poly_points;
    extract_inside_vertices_from_function(m_tooth_poly, poly_points);
    
    HEdgeIntMap he_point_map;
    int count = extract_isovertices_from_function(m_tooth_poly, m_vert_smooth_geodesic_map,
    isovalue, poly_points, he_point_map);

    FacetList candidate_facets;
    for (Facet_handle fd : faces(m_tooth_poly))
        candidate_facets.push_back(fd);

    int index = 0;
    VertIntMap m_vert_ind_map;
    for (Vertex_handle vd : vertices(m_tooth_poly))
    {
        m_vert_ind_map.insert({vd, index});
        index++;
    }

    PolyList poly_faces;
    extract_positive_faces_from_function(m_vert_smooth_geodesic_map, isovalue,
    m_vert_ind_map, he_point_map, candidate_facets, poly_faces);

    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, glue_zone);
    CGAL::Polygon_mesh_processing::keep_largest_connected_components(glue_zone, 1);

    m_pad_outlines = glue_zone;
    // 2. compute the fitting plan
    TriangleList candidate_faces;
    copy_faces_to_triangles(m_pad_outlines, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    Vector normal = fitting_plane.orthogonal_vector();
    normal = normal / std::sqrt(normal.squared_length());

    //decide the direction of the normal
    Vector avg = compute_average_normal(m_pad_outlines);
    if(avg*normal<0)
    {
        normal = -normal;
    }
    // find the tooth center and the projection to the fitting plan
    Point tooth_center = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());
    Point plane_projected_center = fitting_plane.projection(tooth_center);

    // 3. find the farthest point of the outline w.r.t the center
    AABBTree  distance_tree;
    distance_tree.insert(faces(m_pad_outlines).first, faces(m_pad_outlines).second, m_pad_outlines);
    distance_tree.accelerate_distance_queries();
    // auto farthest = distance_tree.farest(tooth_center);

    double max_distance = -1;
    Point farthest_point;
    // Iterate over the vertices
    for (auto v = m_pad_outlines.vertices_begin(); v != m_pad_outlines.vertices_end(); ++v)
    {
        double distance = CGAL::squared_distance(tooth_center, v->point());
        if (distance > max_distance)
        {
            max_distance = distance;
            farthest_point = v->point();
        }
    }

    // 4. project the farthest point to the fitting plan
    Point plane_projected_farthest = fitting_plane.projection(farthest_point);
    Vector principal_direction = plane_projected_farthest - plane_projected_center;
    Vector second_direction = CGAL::cross_product(principal_direction, normal);
    principal_direction = principal_direction / std::sqrt(principal_direction.squared_length());
    second_direction = second_direction / std::sqrt(second_direction.squared_length());

    // 5. --extrusion of glue zone--
    double glue_zone_depth = get_variable_float(BONDING_ZONE_DEPTH, 0.15);
    Polyhedron glue_zone_pad;
    compute_extrusion(glue_zone, glue_zone_pad, normal, glue_zone_depth);
    //mesh_export(glue_zone_pad, "glue_zone_pad.stl");

    // 6. construct a point grid
    double grid_size = get_variable_float(BONDING_GRID_SIZE, 0.5);
    double cone_base_radius = get_variable_float(BONDING_CONE_BASE_RADIUS, 0.15);
    double cone_height = get_variable_float(BONDING_CONE_HEIGHT, 1.5);
    int cone_seg = get_variable_float(BONDING_CONE_SEGMENTS, 20);
    double cone_surface_tolerence = get_variable_float(BONDING_CONE_SURFACE_TELERENCE, 0.05);
    Polyhedron cones_collection;
    max_distance = std::sqrt(max_distance);
    for(double i = -max_distance; i <= max_distance; i += grid_size)
    {
        for (double j = -max_distance; j <= max_distance; j += grid_size)
        {
            Point point_on_plane = plane_projected_center + i * principal_direction + j * second_direction;
            // Line line(point_on_plane, point_on_plane + normal);
            Ray ray(point_on_plane, normal);
            auto intersection = distance_tree.first_intersection(ray);
            if (!intersection)
            {
                Ray ray2(point_on_plane, -normal);
                intersection = distance_tree.first_intersection(ray2);
            }
            if (intersection)
            {
                if (boost::get<Point>(&(intersection->first)))
                {
                    Point intersection_point = *(boost::get<Point>(&(intersection->first)));

                    // Retrieve the face handle of the intersected triangle
                    auto closest = distance_tree.closest_point_and_primitive(intersection_point);

                    Facet_handle face = closest.second; // Primitive containing the closest point

                    // Assuming the mesh is triangulated
                    Halfedge_facet_circulator he = face->facet_begin();
                    Vertex_handle va = he->vertex();
                    Vertex_handle vb = he->next()->vertex();
                    Vertex_handle vc = he->next()->next()->vertex();

                    // Calculate the normal
                    Vector direction_face = CGAL::cross_product(vb->point() - va->point(), 
                    vc->point() - va->point());

                    construct_cone(cones_collection, intersection_point - cone_surface_tolerence * normal,
                                   principal_direction, direction_face, normal, cone_base_radius, cone_height, cone_seg, glue_zone_pad);
                    std::cout << "cones_collection facets: " << cones_collection.size_of_facets() << std::endl;
                }
            }
        }
    }
    
    // save the cone matrix
    //mesh_export(cones_collection, "cones_collection.stl");
    CGAL::alpha_wrap_3(cones_collection, 0.05, 0.05, cones_collection);
    hole_filling(m_pad_shape, 0);
    if (CGAL::Polygon_mesh_processing::does_self_intersect(m_pad_shape))
        CGAL::alpha_wrap_3(m_pad_shape, 0.05, 0.01, m_pad_shape);

    Polyhedron glue_zone_pad_intersection;
    bool result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(glue_zone_pad, cones_collection, glue_zone_pad_intersection);
    hole_filling(glue_zone_pad_intersection, 0);
    //mesh_export(glue_zone_pad_intersection, "glue_zone_pad_intersection.stl");
    std::cout << "glue_zone_pad_intersection facets: " << glue_zone_pad_intersection.size_of_facets() << " is manifold: " << result_code << std::endl;
    result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(m_pad_shape, glue_zone_pad_intersection, glue_zone_pad_intersection);
    hole_filling(glue_zone_pad_intersection, 0);
    //mesh_export(glue_zone_pad_intersection, "glue_zone_pad_intersection2.stl");
    std::cout << "glue_zone_pad_intersection facets: " << glue_zone_pad_intersection.size_of_facets() << " is manifold: " << result_code << std::endl;
    //mesh_export(glue_zone_pad_intersection, "glue_zone_pad_intersection.stl"); 
    // Polyhedron poly_out_wrap;
    //CGAL::alpha_wrap_3(glue_zone_pad_intersection, 0.01, 0.01, glue_zone_pad_intersection);
    //std::cout << "new_pad glue_zone_pad_intersection: " << glue_zone_pad_intersection.size_of_facets() << std::endl;

    // 7. substract pad by the glue zone = new pad
    Polyhedron new_pad;
    // move glue zone with small delta to the normal direction
    for (Polyhedron::Vertex_iterator it = glue_zone_pad.vertices_begin(); it != glue_zone_pad.vertices_end(); ++it) {
        it->point() = it->point() - normal*0.02;
    }
    result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(m_pad_shape, glue_zone_pad, new_pad);
     hole_filling(new_pad, 0);
    //Polyhedron poly_out_wrap;
    std::cout << "new_pad facets: " << new_pad.size_of_facets() << " is manifold: " << result_code << std::endl;
    //mesh_export(new_pad, "new_pad.stl");

    // hole_filling(new_pad, 0);
    //if(result_code != true)
   
    //std::cout << "new_pad facets: " << new_pad.size_of_facets() << std::endl;
    // save the cone matrix
    // mesh_export(new_pad, "new_pad.stl");

    // // 8 intersection of cone collection and the glue area = trimmed cone structure
    // Polyhedron trimmed_glue_zone_cone_structure;
    // CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(cones_collection, glue_zone_pad,
    //                          trimmed_glue_zone_cone_structure);    

    //std::cout << "trimmed_glue_zone_cone_structure facets: " << trimmed_glue_zone_cone_structure.size_of_facets() << std::endl;
    // 9. union of new pad and trimmed cone structure = final pad
    //mesh_export(glue_zone_pad_intersection, "glue_zone_pad_intersection.stl");
    //mesh_export(new_pad, "new_pad.stl");
    CGAL::alpha_wrap_3(new_pad, 0.05, 0.01, new_pad);
    //CGAL::alpha_wrap_3(glue_zone_pad_intersection, 0.05, 0.01, glue_zone_pad_intersection);
    result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_union(glue_zone_pad_intersection,
                                                                            new_pad, output_bonding_surface);

    hole_filling(output_bonding_surface, 0);
    std::cout << "output_bonding_surface facets: " << output_bonding_surface.size_of_facets() << " is manifold: " << result_code << std::endl;
    
    if(output_bonding_surface.size_of_facets() == 0)
    {
        std::cout << "output_bonding_surface is empty" << std::endl;
        CGAL::alpha_wrap_3(new_pad, 0.05, 0.01, new_pad);
        CGAL::Polygon_mesh_processing::corefine_and_compute_union(glue_zone_pad_intersection,
                                                                  new_pad, output_bonding_surface);
    }

    // if (solve_selfintersection(output_bonding_surface, poly_out_wrap))
    // {
    //     copy_polyehdral_surface(poly_out_wrap, output_bonding_surface);
    // }

    std::cout << "output_bonding_surface facets: " << output_bonding_surface.size_of_facets() << std::endl;
    // save the cone matrix
    //mesh_export(output_bonding_surface, "output_bonding_surface.stl");
}


/*void Tooth::compute_pad_shape(Polyhedron& wire_outlines, Polyhedron& wire_shape, double pad_height, double wire_height, Point& teeth_center)
{
    if (m_pad_outlines.empty())
        return;

    if (wire_height <= pad_height)
    {
        std::cout << "Bad parameter!" << std::endl;
        return;
    }

    m_pad_shape.clear();
    m_pad_left.clear();
    m_pad_right.clear();
    m_pad_split_flag = false;

    TriangleList candidate_faces;
    copy_faces_to_triangles(m_pad_outlines, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    Vector normal = fitting_plane.orthogonal_vector();
    normal = normal / std::sqrt(normal.squared_length());
    Point tooth_center = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());
    if((teeth_center - tooth_center) * normal < 0)
        normal = -normal;
    
    compute_extrusion(m_pad_outlines, m_pad_shape, normal, pad_height);
    
    //_compute_bonding_surface(m_pad_shape, m_pad_outlines);

    Polyhedron wire_extrusion;
    //compute_extrusion(m_pad_outlines, m_pad_shape, pad_height);
    compute_extrusion(wire_outlines, wire_extrusion, wire_height);
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(m_pad_shape, wire_extrusion, m_pad_shape);

    Polyhedron poly_out_wrap;
    if (solve_selfintersection(m_pad_shape, poly_out_wrap))
    {
        copy_polyehdral_surface(poly_out_wrap, m_pad_shape);
    }

    //debug
    // char buffer[10];
    // _itoa(rand()*1000, buffer, 10);
    // std::string filename_i = "DEBUG_tooth";
    // filename_i = filename_i + buffer + ".stl";
    // mesh_export(wire_extrusion, filename_i);

    if (solve_selfintersection(wire_shape, poly_out_wrap))
    {
        copy_polyehdral_surface(poly_out_wrap, wire_shape);
    }
    //smooth the pad
    //CGAL::Polygon_mesh_processing::angle_and_area_smoothing(m_pad_shape, CGAL::parameters::number_of_iterations(50)
    //    .use_safety_constraints(false));
    //std::cout << "Polygon_mesh_processing::angle_and_area_smoothing" << std::endl;

    
    if (solve_selfintersection(m_pad_shape, poly_out_wrap))
    {
        copy_polyehdral_surface(poly_out_wrap, m_pad_shape);
    }
    // mesh_export(m_pad_shape, "pad_output.stl");
    // smooth the mesh
    smooth_mesh(m_pad_shape, get_variable_float(PAD_SMOOTH_ITERATIONS,400), m_pad_outlines);

    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(m_pad_shape, wire_shape, m_pad_shape);

    std::cout << "Computed pad shape with " << m_pad_shape.size_of_vertices() << " vertices and " << m_pad_shape.size_of_facets() << " faces!" << std::endl;
}*/

void Tooth::compute_pad_shape(Polyhedron& wire_shape, Polyhedron& enlarge_wire_shape, double pad_height, double wire_height, Point& teeth_center, double wire_radius)
{
    if (m_pad_outlines.empty())
        return;

    if (wire_height <= pad_height)
    {
        std::cout << "Bad parameter!" << std::endl;
        return;
    }

    m_pad_shape.clear();
    m_pad_left.clear();
    m_pad_right.clear();
    m_pad_split_flag = false;

    TriangleList candidate_faces;
    copy_faces_to_triangles(m_pad_outlines, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    Vector normal = fitting_plane.orthogonal_vector();
    normal = normal / std::sqrt(normal.squared_length());
    Point tooth_center = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());
    if((teeth_center - tooth_center) * normal < 0)
        normal = -normal;

    AABBTree enlarge_wire_tree;
    enlarge_wire_tree.insert(faces(enlarge_wire_shape).first, faces(enlarge_wire_shape).second, enlarge_wire_shape);
    enlarge_wire_tree.accelerate_distance_queries();

    VertDoubleMap pad_heights;
    for(Vertex_handle vd: vertices(m_pad_outlines))
    {
        Point source = vd->point() + normal * 3.;
        Ray ray(source, vd->point());
        int number = (int)enlarge_wire_tree.number_of_intersected_primitives(ray);

        if (number <= 0)
            pad_heights.insert({vd, pad_height});
        else
        {
            auto intersection = enlarge_wire_tree.first_intersection(ray);
            if(intersection && boost::get<Point>(&(intersection->first)))
            {
                const Point* p =  boost::get<Point>(&(intersection->first));
                double height = std::max(pad_height, std::sqrt(CGAL::squared_distance(*p, vd->point())));
                pad_heights.insert({vd, height});
            }
            else
                pad_heights.insert({vd, pad_height});
        }
    }
    m_pad_extrusion_direction = normal;
    compute_extrusion(m_pad_outlines, m_pad_shape, normal, pad_heights);

    // smooth the mesh
    smooth_mesh(m_pad_shape, get_variable_float(PAD_SMOOTH_ITERATIONS, 10), m_pad_outlines);

    // swept volume
    Polyhedron curr_wire, offset_wire, swept_wire;
    copy_polyehdral_surface(wire_shape, curr_wire);
    copy_polyehdral_surface(wire_shape, offset_wire);
    int step = 3;
    double delta = -2. * wire_radius / step;
    Vector offset = normal * delta;
    for(int i = 1; i < step; i++)
    {
        for(Vertex_handle vd: vertices(offset_wire))
        {
            Point p = vd->point() + offset;
            vd->point() = p ;
        }
            
        
        merge_polyehdral_surface(curr_wire, offset_wire, swept_wire);
        curr_wire.clear();
        copy_polyehdral_surface(swept_wire, curr_wire);
    }

    Polyhedron swept_wire_wrap;
    solve_selfintersection(swept_wire, swept_wire_wrap);

    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(m_pad_shape, swept_wire_wrap, m_pad_shape);

    Polyhedron pad_shape_wrap;
    if (solve_selfintersection(m_pad_shape, pad_shape_wrap))
    {
        copy_polyehdral_surface(pad_shape_wrap, m_pad_shape);
    }

    // compute pad height map
    /*SegmentTree wire_tree(wire_outlines.begin(), wire_outlines.end());
    VertDoubleMap pad_heights;
    double max_dist = wire_radius * 1.732;
    double min_dist = wire_radius;
    for(Vertex_handle vd: vertices(m_pad_outlines))
    {
        double dist = std::min(max_dist, std::max(min_dist, std::sqrt(wire_tree.squared_distance(vd->point()))));
        double height = (max_dist - dist) / (max_dist - min_dist) * wire_height + (dist - min_dist) / (max_dist - min_dist) * pad_height; 
        pad_heights.insert({vd, height});
    }

    compute_extrusion(m_pad_outlines, m_pad_shape, normal, pad_heights);

    // smooth the mesh
    smooth_mesh(m_pad_shape, get_variable_float(PAD_SMOOTH_ITERATIONS, 400), m_pad_outlines);

    Polyhedron swept_wire;
    compute_swept_extrusion(wire_shape, swept_wire, normal, -10. * wire_radius);

    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(m_pad_shape, swept_wire, m_pad_shape);

    Polyhedron pad_shape_wrap;
    if (solve_selfintersection(m_pad_shape, pad_shape_wrap))
    {
        copy_polyehdral_surface(pad_shape_wrap, m_pad_shape);
    }*/

    std::cout << "Computed pad shape with " << m_pad_shape.size_of_vertices() << " vertices and " << m_pad_shape.size_of_facets() << " faces!" << std::endl;
}


void Tooth::compute_bonding_surface(const std::string& pad, const std::string& padoutline , const std::string& padout)
{
    Polyhedron _m_pad_shape;
    Polyhedron _m_pad_outlines;
    Polyhedron _m_pad_out;
    if(!mesh_import(pad, _m_pad_shape)) 
        return; 
    if(!mesh_import(padoutline, _m_pad_outlines)) 
        return; 

    std::cout << "m_pad_outlines facets: " << _m_pad_outlines.size_of_facets() << std::endl;
    std::cout << "m_pad_shape facets: " << _m_pad_shape.size_of_facets() << std::endl;

    _compute_bonding_surface(_m_pad_shape, _m_pad_outlines, _m_pad_out);
    if (padout.length() > 0)
        mesh_export(_m_pad_out, padout);
}

/*void Tooth::compute_splited_pad_shape(Polyhedron& left_wire_outlines, Polyhedron& left_wire_shape,
    Polyhedron& right_wire_outlines, Polyhedron& right_wire_shape,
    double pad_height, double wire_height, Point& teeth_center)
{
    if (m_pad_outlines.empty())
        return;

    if (wire_height <= pad_height)
    {
        std::cout << "Bad parameter!" << std::endl;
        return;
    }

    m_pad_shape.clear();
    m_pad_left.clear();
    m_pad_right.clear();
    m_pad_split_flag = true;

    TriangleList candidate_faces;
    copy_faces_to_triangles(m_pad_outlines, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    Vector normal = fitting_plane.orthogonal_vector();
    normal = normal / std::sqrt(normal.squared_length());
    Point tooth_center = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());
    if((teeth_center - tooth_center) * normal < 0)
        normal = -normal;

    compute_extrusion(m_pad_outlines, m_pad_shape, normal, pad_height);
    //compute_extrusion(m_pad_outlines, m_pad_shape, pad_height);
    //_compute_bonding_surface(m_pad_shape, m_pad_outlines);// bonding structures

    Polyhedron left_wire_extrusion, right_wire_extrusion;
    compute_extrusion(left_wire_outlines, left_wire_extrusion, wire_height);
    compute_extrusion(right_wire_outlines, right_wire_extrusion, wire_height);

    CGAL::Polygon_mesh_processing::corefine_and_compute_union(m_pad_shape, left_wire_extrusion, m_pad_shape);
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(m_pad_shape, right_wire_extrusion, m_pad_shape);

    smooth_mesh(m_pad_shape, get_variable_float(PAD_SMOOTH_ITERATIONS,400), m_pad_outlines);

    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(m_pad_shape, left_wire_shape, m_pad_shape);
    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(m_pad_shape, right_wire_shape, m_pad_shape);

    // after smooth, need to solve intersection
    Polyhedron poly_out_wrap;
    if (solve_selfintersection(m_pad_shape, poly_out_wrap))
    {
        copy_polyehdral_surface(poly_out_wrap, m_pad_shape);
    }

    Point left_center = CGAL::ORIGIN;
    double ratio = 1. / left_wire_outlines.size_of_vertices();
    for (Vertex_handle vd : vertices(left_wire_outlines))
        left_center = left_center + ratio * (vd->point() - CGAL::ORIGIN);

    copy_polyehdral_surface(m_pad_shape, m_pad_left);
    copy_polyehdral_surface(m_pad_shape, m_pad_right);

    if ((left_center - m_cut_plane.projection(left_center)) * m_cut_plane.orthogonal_vector() < 0.)
    {
        CGAL::Polygon_mesh_processing::clip(m_pad_left, m_cut_plane, CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
        CGAL::Polygon_mesh_processing::clip(m_pad_right, m_cut_plane.opposite(),CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
    }
    else
    {
        CGAL::Polygon_mesh_processing::clip(m_pad_right, m_cut_plane, CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
        CGAL::Polygon_mesh_processing::clip(m_pad_left, m_cut_plane.opposite(), CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
    }

    std::cout << "Computed left pad shape with " << m_pad_left.size_of_vertices() << " vertices and " << m_pad_left.size_of_facets() << " faces!" << std::endl;
    std::cout << "Computed right pad shape with " << m_pad_right.size_of_vertices() << " vertices and " << m_pad_right.size_of_facets() << " faces!" << std::endl;
}*/

void Tooth::compute_splited_pad_shape(  Polyhedron& left_wire_shape, Polyhedron& left_enlarge_wire_shape, 
                                        Polyhedron& right_wire_shape, Polyhedron& right_enlarge_wire_shape,
                                        double pad_height, double wire_height, Point& teeth_center, double wire_radius)
{
    if (m_pad_outlines.empty())
        return;

    m_pad_shape.clear();
    m_pad_left.clear();
    m_pad_right.clear();
    m_pad_split_flag = false;

    TriangleList candidate_faces;
    copy_faces_to_triangles(m_pad_outlines, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    Vector normal = fitting_plane.orthogonal_vector();
    normal = normal / std::sqrt(normal.squared_length());
    Point tooth_center = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());
    if((teeth_center - tooth_center) * normal < 0)
        normal = -normal;

    // compute pad height map
    Polyhedron enlarge_wire_shape;
    merge_polyehdral_surface(left_enlarge_wire_shape, right_enlarge_wire_shape, enlarge_wire_shape);

    AABBTree enlarge_wire_tree;
    enlarge_wire_tree.insert(faces(enlarge_wire_shape).first, faces(enlarge_wire_shape).second, enlarge_wire_shape);
    enlarge_wire_tree.accelerate_distance_queries();

    VertDoubleMap pad_heights;
    for(Vertex_handle vd: vertices(m_pad_outlines))
    {
        Point source = vd->point() + normal * 3.;
        Ray ray(source, vd->point());
        int number = (int)enlarge_wire_tree.number_of_intersected_primitives(ray);

        if (number <= 0)
            pad_heights.insert({vd, pad_height});
        else
        {
            auto intersection = enlarge_wire_tree.first_intersection(ray);
            if(intersection && boost::get<Point>(&(intersection->first)))
            {
                const Point* p =  boost::get<Point>(&(intersection->first));
                double height = std::max(pad_height, std::sqrt(CGAL::squared_distance(*p, vd->point())));
                pad_heights.insert({vd, height});
            }
            else
                pad_heights.insert({vd, pad_height});
        }
    }
    m_pad_extrusion_direction = normal;
    compute_extrusion(m_pad_outlines, m_pad_shape, normal, pad_heights);

    // smooth the mesh
    smooth_mesh(m_pad_shape, get_variable_float(PAD_SMOOTH_ITERATIONS, 10), m_pad_outlines);
    

    std::cout << "Computed pad shape with " << m_pad_shape.size_of_vertices() << " vertices and " << m_pad_shape.size_of_facets() << " faces!" << std::endl;

    // swept volume
    Polyhedron wire_shape;
    merge_polyehdral_surface(left_wire_shape, right_wire_shape, wire_shape);
    Polyhedron curr_wire, offset_wire, swept_wire;
    copy_polyehdral_surface(wire_shape, curr_wire);
    copy_polyehdral_surface(wire_shape, offset_wire);
    int step = 3;
    double delta = -2. * wire_radius / step;
    Vector offset = normal * delta;
    for(int i = 1; i < step; i++)
    {
        for(Vertex_handle vd: vertices(offset_wire))
        {
            Point p = vd->point() + offset;
            vd->point() = p ;
        }
            
        
        merge_polyehdral_surface(curr_wire, offset_wire, swept_wire);
        curr_wire.clear();
        copy_polyehdral_surface(swept_wire, curr_wire);
    }

    Polyhedron swept_wire_wrap;
    solve_selfintersection(swept_wire, swept_wire_wrap);

    CGAL::Polygon_mesh_processing::corefine_and_compute_difference(m_pad_shape, swept_wire_wrap, m_pad_shape);

    Polyhedron pad_shape_wrap;
    if (solve_selfintersection(m_pad_shape, pad_shape_wrap))
    {
        copy_polyehdral_surface(pad_shape_wrap, m_pad_shape);
    }

    copy_polyehdral_surface(m_pad_shape, m_pad_left);
    copy_polyehdral_surface(m_pad_shape, m_pad_right);
    copy_polyehdral_surface(m_pad_outlines, m_pad_outlines_left);
    copy_polyehdral_surface(m_pad_outlines, m_pad_outlines_right);

    PointList verts;
    extract_inside_vertices_from_function(left_wire_shape, verts);
    Point left_center = CGAL::centroid(verts.begin(), verts.end());

    if ((left_center - m_cut_plane.projection(left_center)) * m_cut_plane.orthogonal_vector() < 0.)
    {
        CGAL::Polygon_mesh_processing::clip(m_pad_left, m_cut_plane, CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
        CGAL::Polygon_mesh_processing::clip(m_pad_right, m_cut_plane.opposite(), CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
        CGAL::Polygon_mesh_processing::clip(m_pad_outlines_left, m_cut_plane, CGAL::Polygon_mesh_processing::parameters::clip_volume(false));
        CGAL::Polygon_mesh_processing::clip(m_pad_outlines_right, m_cut_plane.opposite(), CGAL::Polygon_mesh_processing::parameters::clip_volume(false));
    }
    else
    {
        CGAL::Polygon_mesh_processing::clip(m_pad_right, m_cut_plane, CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
        CGAL::Polygon_mesh_processing::clip(m_pad_left, m_cut_plane.opposite(), CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
        CGAL::Polygon_mesh_processing::clip(m_pad_outlines_left, m_cut_plane.opposite(), CGAL::Polygon_mesh_processing::parameters::clip_volume(false));
        CGAL::Polygon_mesh_processing::clip(m_pad_outlines_right, m_cut_plane, CGAL::Polygon_mesh_processing::parameters::clip_volume(false));
    }

    std::cout << "Computed left pad shape with " << m_pad_left.size_of_vertices() << " vertices and " << m_pad_left.size_of_facets() << " faces!" << std::endl;
    std::cout << "Computed right pad shape with " << m_pad_right.size_of_vertices() << " vertices and " << m_pad_right.size_of_facets() << " faces!" << std::endl;
}



void Tooth::save_pad(const std::string &filename)
{
    mesh_export(m_pad_shape, filename);
}

void Tooth::get_pad_statistics(StringList& names, DoubleList& left, DoubleList& right)
{
    names.push_back(m_name);

    if(m_pad_outlines.size_of_vertices() == 0)
    {
        left.push_back(0.);
        right.push_back(0.);
    }
    else if(m_pad_split_flag)
    {
        Polyhedron pad_1, pad_2;
        copy_polyehdral_surface(m_pad_outlines, pad_1);
        copy_polyehdral_surface(m_pad_outlines, pad_2);

        CGAL::Polygon_mesh_processing::clip(pad_1, m_cut_plane);
        CGAL::Polygon_mesh_processing::clip(pad_2, m_cut_plane.opposite());
        double area_1 = CGAL::Polygon_mesh_processing::area(pad_1);
        double area_2 = CGAL::Polygon_mesh_processing::area(pad_2);

        PointList verts;
        extract_inside_vertices_from_function(pad_1, verts);
        Point center_1 = CGAL::centroid(verts.begin(), verts.end());
        verts.clear();
        extract_inside_vertices_from_function(pad_2, verts);
        Point center_2 = CGAL::centroid(verts.begin(), verts.end());

        if(center_1.x() < center_2.x())
        {
            left.push_back(area_1);
            right.push_back(area_2);
        }
        else
        {
            left.push_back(area_2);
            right.push_back(area_1);
        }
    }
    else
    {
        double area = CGAL::Polygon_mesh_processing::area(m_pad_outlines);
        left.push_back(area);
        right.push_back(0.);
    }
}

} // namespace OrthoLab
