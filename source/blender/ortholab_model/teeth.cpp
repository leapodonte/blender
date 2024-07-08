#include "teeth.h"

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/centroid.h>
#include <CGAL/intersections.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>

#include "types.h"
#include "ramp.h"
#include "solver.h"
#include "meshing.h"
#include "tooth.h"
#include "wire.h"
#include "extrusion.h"
#include "save.h"
#include "mesh_io.h"
#include "project.h"
#include "base_helpers.h"
#include "file_services.h"

#include "tooth.h"
#include "wire.h"

#include "manu_asm_guide.h"
#include "spdlog/spdlog.h"

namespace OrthoLab {

//----------------------------- Initialization -----------------------------//

Teeth::Teeth(const std::string& id, const Project& prj):
    SceneObject(id, prj),
    m_manu_asm_guide(*this)
{
    // std::cout << "Teeth address:" << this << std::endl;
    for (int i = 0; i < 6; i++){
        char buffer[10];
        _itoa(i, buffer, 10);
        m_tooth.push_back(std::make_unique<Tooth>(buffer,*this));
    }
    for (int i = 0; i < 5; i++){
        char buffer[10];
        _itoa(i, buffer, 10);
        m_wires.push_back(std::make_unique<Wire>(buffer,*this));
    }
}

Teeth::~Teeth()
{
    reset();
}

void Teeth::save(FileWriter& file)
{
    //    Polyhedron      m_teeth_poly;
    //    VertIntMap      m_vert_ind_map;
    //    VertDoubleMap   m_vert_ray_map;
    //    VertDoubleMap   m_vert_smooth_ray_map;
    //    VertDoubleMap   m_vert_inside_func_map;
    //
    //    Polyhedron      m_teeth_inside;
    //    VertIntMap      m_inside_vert_ind_map;
    //    VertDoubleMap   m_inside_vert_convex_map;
    //    VertDoubleMap   m_inside_vert_smooth_convex_map;
    //    VertDoubleMap   m_inside_vert_convex_func_map;
    //    AABBTree        m_inside_tree;
    //    KDTree          m_inside_boundary_tree;
    //
    //    ToothList       m_tooth;
    //    WireList        m_wires;
}

void Teeth::load(FileReader& file)
{

}

std::string Teeth::get_base_path() const // for db save and load
{
    return get_project().get_case_name() +
           PATH_SEPERATION + get_id();
}

void Teeth::save_to_db()
{
    std::string teeth_path = get_base_path();
    teeth_path = ensure_directory_exists(teeth_path);
    // export meshes
    teeth_path += PATH_SEPERATION;
    
    mesh_export(m_teeth_poly_expanded, teeth_path + "teeth_poly_expanded.stl");
    
    savePointsToXYZ(m_flattened_wire, teeth_path + "flattened_wire.xyz");

#pragma omp parallel for
    for (int i = 0; i < m_tooth.size(); i++)
        m_tooth[i]->save_to_db();
#pragma omp parallel for
    for (int i = 0; i < m_wires.size(); i++)
        m_wires[i]->save_to_db();


}
void Teeth::load_from_db()
{
    std::string teeth_path = get_base_path();
    teeth_path = ensure_directory_exists(teeth_path);
    // export meshes
    teeth_path += PATH_SEPERATION;
    
    std::string mesh_file = teeth_path + "teeth_poly_expanded.stl";
    // mesh_export(m_teeth_poly_expanded, mesh_file);
    mesh_import(mesh_file, m_teeth_poly);

    mesh_file = teeth_path + "flattened_wire.xyz";
	readPointsFromXYZ(mesh_file, m_flattened_wire);

#pragma omp parallel for
    for (int i = 0; i < m_tooth.size(); i++)
        m_tooth[i]->load_from_db();
#pragma omp parallel for
    for (int i = 0; i < m_wires.size(); i++)
        m_wires[i]->load_from_db();
}

void Teeth::save_teeth_model(const std::string &filename)
{
    mesh_export(m_teeth_poly, filename);
}

void Teeth::reset()
{
    m_teeth_poly.clear();
    m_teeth_distance_map.clear();
    m_vert_ray_map.clear();
    m_vert_ind_map.clear();
    m_vert_smooth_ray_map.clear();
    m_vert_inside_func_map.clear();

    reset_inside();

    for (int i = 0; i < m_tooth.size(); i++)
        m_tooth[i]->reset();

    for (int i = 0; i < m_wires.size(); i++)
        m_wires[i]->reset();
}

void Teeth::reset_inside()
{
    m_teeth_inside.clear();
    m_inside_vert_ind_map.clear();
    m_inside_vert_convex_map.clear();
    m_inside_vert_smooth_convex_map.clear();
    m_inside_vert_convex_func_map.clear();
    m_inside_tree.clear();
    m_inside_boundary_tree.clear();
}

void Teeth::reset_convexity()
{
    m_inside_vert_convex_map.clear();
    m_inside_vert_smooth_convex_map.clear();
    m_inside_vert_convex_func_map.clear();
}

bool Teeth::load_teeth_mesh(std::string filename)
{
    reset();

    if(!mesh_import(filename, m_teeth_poly)) return false;

    spdlog::trace("Load teeth mesh with {0}  vertices and {1} faces.", m_teeth_poly.size_of_vertices(),
        m_teeth_poly.size_of_facets()); 

    initialize_vertex_indices();

    return true;
}

bool Teeth::init_teeth_mesh(const float* vertices, size_t vertices_size, const int* triangles, size_t triangles_size)
{
    reset();

    PointList poly_points;
    PolyList poly_faces;

    for (size_t vi = 0; vi < vertices_size; ++vi)
    {
        const float* v = vertices + (vi * 3);
        Point p(v[0], v[2], v[3]);
        poly_points.push_back(p);
    }

    IntList face_inds;
    face_inds.resize(3);
    for (size_t ti = 0; ti < triangles_size; ++ti)
    {
        const int* tri = triangles + ti * 3;
        face_inds[0] = tri[0];
        face_inds[1] = tri[1];
        face_inds[2] = tri[2];
        poly_faces.push_back(face_inds);
    }

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, m_teeth_poly);

    spdlog::trace("Init teeth mesh with {0}  vertices and {1} faces.", m_teeth_poly.size_of_vertices(),
        m_teeth_poly.size_of_facets());
    initialize_vertex_indices();

    return true;
}

double compute_length_for_hole(Halfedge_handle he, Polyhedron& mesh)
{
    double length = 0.;

    for (Halfedge_handle hc : CGAL::halfedges_around_face(he, mesh))
        length += std::sqrt(CGAL::squared_distance(hc->vertex()->point(), hc->opposite()->vertex()->point()));

    return length;
}

void Teeth::expandTeethMesh(double dist)
{   
    // copy_polyehdral_surface(m_teeth_poly, m_teeth_poly_expanded);
     
    // expand_mesh(m_teeth_poly_expanded, dist);
    compute_enlarge_msh(m_teeth_poly, m_teeth_poly_expanded,dist);
}



void Teeth::flattenWires(double dist)
{
    spdlog::trace("Teeth::flattenWires:{}", dist);
    // std::cout << "flattenWires first wire " << std::endl;
    std::vector<Transformation> trs;
    Transformation prev_trans = Transformation(CGAL::TRANSLATION, Vector(0, 0, 0));
    m_flattened_wire.clear();

    // push back the frist wire
    auto first_poly_line = m_wires[0]->get_polylines();
    if (first_poly_line.size() == 0)
        return; // in case the wires are absent
    for (int i = first_poly_line.size()-1; i>=0; --i)
    {
        auto point = first_poly_line[i].target();
        m_flattened_wire.push_back(point);
        spdlog::trace("first_poly_line target:{0} {1} {2}", point.x(), point.y(), point.z());
    }
    m_flattened_wire.push_back(first_poly_line[0].source());
    // std::cout << m_flattened_wire.back() << std::endl;
     m_wires[0]->setTransformation(prev_trans);
     int valid_wire_num = 1;
    for (int i = 0; i < m_wires.size() - 1; ++i)
    {
        spdlog::trace("flattenWires wire:{}", i + 1);
        // std::cout << "flattenWires wire " << i + 1 << ":"<<std::endl; 
        if (m_wires[i + 1]->get_polylines().size() == 0)
            break; // the last wire is empty
        // prolonge the last segment of the first polyline
        auto first_segment = m_wires[i]->get_polylines()[0];
        auto first_seg_dir = first_segment.source() - first_segment.target();
        first_seg_dir = first_seg_dir / std::sqrt(first_seg_dir.squared_length());
        Segment prolonged = Segment(first_segment.target(),
                                    first_segment.source() + first_seg_dir * dist);
        auto second_segment = m_wires[i + 1]->get_polylines().back();
        Segment connection_next = Segment(second_segment.target(), second_segment.source());
        // local transformation
        auto tr = compute_alignment_transform(m_wires[i]->get_cut_plane().m_plane,
                                              prolonged,
                                              m_wires[i + 1]->get_cut_plane().m_plane,
                                              connection_next);
        // global transformation
        prev_trans = prev_trans * tr;
        auto polyline = m_wires[i + 1]->get_polylines();
        for (int j = polyline.size() - 1; j >= 0; --j)
        {
            m_flattened_wire.push_back(prev_trans.transform(polyline[j].target()));
            // std::cout << m_flattened_wire.back() << std::endl;
        }
        m_flattened_wire.push_back(prev_trans.transform(polyline[0].source()));
        m_wires[i + 1]->setTransformation(prev_trans);
        valid_wire_num++;
    }
    
    spdlog::trace("valid_wire_num:{0}", valid_wire_num);
    // flatten pads  
	for (int i = 0; i < valid_wire_num + 1; i++)
    {
        if (m_tooth[i]->get_init_flag())
        { 
            if (i != 0) // i=1,2,3,4,5
            {
                m_tooth[i]->set_left_transformation(m_wires[i - 1]->getTransformation());                
            }
            if (i == 0)
            {
                m_tooth[i]->set_left_transformation(Transformation());
            }
			/*else if (i == m_tooth.size() - 1)
            {
                m_tooth[i]->set_left_transformation( m_wires[i - 1]->getTransformation());
            }*/
            else if(i != valid_wire_num)// i=1,2,3,4
            {
                m_tooth[i]->set_left_transformation( m_wires[i - 1]->getTransformation());
                m_tooth[i]->set_right_transformation( m_wires[i]->getTransformation());
            }
        }
    }
}

void Teeth::close_teeth_mesh()
{ 
    // if not empty then return
    if(!m_teeth_poly_closed.is_empty())
        return;
    // check if m_teeth_poly is empty;
    copy_polyehdral_surface(m_teeth_poly, m_teeth_poly_closed);
    // hole filling
    HEdgeList border_cycles;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m_teeth_poly_closed, std::back_inserter(border_cycles));

    if (border_cycles.size() > 0) // contain holes
    {
        DoubleList border_lengths;
        for (int i = 0; i < border_cycles.size(); i++)
        {
            double length = compute_length_for_hole(border_cycles[i], m_teeth_poly_closed);
            border_lengths.push_back(length);
        }

        unsigned int nb_holes = 0;
        std::vector<int> indices(border_lengths.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool { return border_lengths[A] > border_lengths[B]; });

        for (int i = 0; i < indices.size(); i++)
        {
            FacetList  patch_facets;
            VertexList patch_vertices;
            bool success = std::get<0>(CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(m_teeth_poly,
                border_cycles[indices[i]],
                CGAL::parameters::face_output_iterator(std::back_inserter(patch_facets))
                .vertex_output_iterator(std::back_inserter(patch_vertices))
                .vertex_point_map(get(CGAL::vertex_point, m_teeth_poly_closed))
                .geom_traits(Kernel())));
            ++nb_holes;
        }
        spdlog::trace("Filled:{0} holes!", nb_holes);
        // std::cout << "    Filled " << nb_holes << "holes!" << std::endl;
    }

    initialize_vertex_indices();
}

Polyhedron& Teeth::get_teeth_poly() { return m_teeth_poly; }
Polyhedron& Teeth::get_closed_teeth_poly()
{
    if (m_teeth_poly_closed.is_empty())
    {
        close_teeth_mesh();
    }
    return m_teeth_poly_closed;
}
VertDoubleMap& Teeth::get_teeth_distance_map() { return m_teeth_distance_map; }
VertDoubleMap& Teeth::get_smoothed_ray_map() { return m_vert_smooth_ray_map; }
VertDoubleMap& Teeth::get_vertex_inside_function_map() { return m_vert_inside_func_map; }
Polyhedron& Teeth::get_inside_teeth_poly() { return m_teeth_inside; }
VertDoubleMap& Teeth::get_smoothed_convex_map() { return m_inside_vert_smooth_convex_map; }
VertDoubleMap& Teeth::get_vertex_convex_function_map() { return m_inside_vert_convex_func_map; }

int Teeth::size_of_facets() { return static_cast<int>(m_teeth_poly.size_of_facets()); }
int Teeth::size_of_vertices() { return static_cast<int>(m_teeth_poly.size_of_vertices()); }
int Teeth::size_of_inside_facets() { return static_cast<int>(m_teeth_inside.size_of_facets()); }
int Teeth::size_of_inside_vertices() { return static_cast<int>(m_teeth_inside.size_of_vertices()); }

void Teeth::read_polygon_vertices(PointList& points)
{
    for (Vertex_handle vd : vertices(m_teeth_poly))
        points.push_back(vd->point());
}

void Teeth::initialize_vertex_indices()
{
    m_vert_ind_map.clear();
    int index = 0;

    for (Vertex_handle vd : vertices(m_teeth_poly))
    {
        m_vert_ind_map.insert({ vd, index });
        index++;
    }
}

double Teeth::get_zmax()
{
    double zmax = -1e10;

    for (Vertex_handle vd : vertices(m_teeth_poly))
        zmax = ortho_max(zmax, (vd->point()).z());

    return zmax;
}

void Teeth::compute_teeth_distance(Polyhedron& other_teeth)
{
    m_teeth_distance_map.clear();

    AABBTree other_tree;
    other_tree.insert(faces(other_teeth).first, faces(other_teeth).second, other_teeth);
    other_tree.accelerate_distance_queries();

    for(Vertex_handle vd : vertices(m_teeth_poly))
    {
        double dist = std::sqrt(other_tree.squared_distance(vd->point()));
        m_teeth_distance_map.insert({vd, dist});
    }
}

//----------------------------- Inside Function -----------------------------//

bool Teeth::compute_laplacian_based_inside(AABBTree& teeth_tree, Point& center, int smooth_range, double lambda)
{
    spdlog::trace("Teeth::compute_laplacian_based_inside center:{0} {1} {2} smooth_range: {3} lambda: {4}", 
        center.x(), center.y(), center.z(), smooth_range, lambda);
    // Initialize intersection field
    compute_intersection_ray(teeth_tree, center);

    // Smoothed intersection field
    compute_smoothed_field(smooth_range, m_vert_ray_map, m_vert_smooth_ray_map);

    // Initialize linear system
    int nb_variables = size_of_vertices();
    ESMatrix L(nb_variables * 2, nb_variables);
    EVector B(nb_variables * 2), X(nb_variables);

    ESTripleList LT;
    LT.reserve(12 * nb_variables);

    // Assemble system
    IntDoubleMap smoothed_ray_map;
    for (Vertex_handle vd : vertices(m_teeth_poly))
        smoothed_ray_map.insert({ m_vert_ind_map[vd], m_vert_smooth_ray_map[vd] });

    assemble_laplacian_matrix(m_teeth_poly, LT, B, m_vert_ind_map);
    assemble_constraint_matrix(lambda, LT, B, smoothed_ray_map, nb_variables);

    L.setFromTriplets(LT.begin(), LT.end());
    LT.clear();
    bool flag_solver = solve_laplacian(L, B, X);

    if (!flag_solver)
        return false;

    // Assign value
    m_vert_inside_func_map.clear();
    for (Vertex_handle vd : vertices(m_teeth_poly))
        m_vert_inside_func_map.insert({ vd, X[m_vert_ind_map[vd]] });

    return true;
}

void Teeth::compute_intersection_ray(AABBTree& teeth_tree, Point& center)
{
    if (!m_vert_ray_map.empty())
        return;

    for (Vertex_handle vd : vertices(m_teeth_poly))
    {
        Vector direction = vd->point() - center;
        Segment segment(center, center + direction * 1.00001);
        int number = (int)teeth_tree.number_of_intersected_primitives(segment);

        if (number == 1)
            m_vert_ray_map[vd] = 1.;
        else
            m_vert_ray_map[vd] = 0.;
    }
}

double Teeth::compute_distance_to_inside_teeth(Point& point)
{
    if (m_inside_tree.empty())
        return -1.;

    double dist = std::sqrt(m_inside_tree.squared_distance(point));
    return dist;
}

void Teeth::extract_inside_from_function(double value)
{
    if (m_vert_inside_func_map.empty())
        return;

    reset_inside();

    // Fill holes
    VertDoubleMap vert_inside_func_fill_map;
    fill_holes_in_vertmap(m_teeth_poly, m_vert_inside_func_map, value, vert_inside_func_fill_map);

    PointList poly_points;
    extract_inside_vertices_from_function(m_teeth_poly, poly_points);

    HEdgeIntMap he_point_map;
    int count = extract_isovertices_from_function(m_teeth_poly, vert_inside_func_fill_map, value, poly_points, he_point_map);

    FacetList candidate_facets;
    for (Facet_handle fd : faces(m_teeth_poly))
        candidate_facets.push_back(fd);

    PolyList  poly_faces;
    extract_positive_faces_from_function(vert_inside_func_fill_map, value, m_vert_ind_map, he_point_map, candidate_facets, poly_faces);

    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, m_teeth_inside);
    CGAL::Polygon_mesh_processing::keep_largest_connected_components(m_teeth_inside, 1);
    // std::cout << "  Compute inside mesh with " << m_teeth_inside.size_of_vertices() << " vertices and " << m_teeth_inside.size_of_facets() << " faces." << std::endl;
    spdlog::trace("Compute inside mesh with {0}  vertices and {1} faces.", m_teeth_inside.size_of_vertices(),
        m_teeth_inside.size_of_facets());
    // initialize vertex index map for inside teeth
    initialize_inside_vertex_index();
    // initialize aabbtree
    m_inside_tree.insert(faces(m_teeth_inside).first, faces(m_teeth_inside).second, m_teeth_inside);
    m_inside_tree.accelerate_distance_queries();
    // initialize boundary kdtree
    initialize_inside_boundary_kdtree();
}

void Teeth::initialize_inside_vertex_index()
{
    if (m_teeth_inside.size_of_vertices() > 0)
    {
        int index = 0;
        for (Vertex_handle vd : vertices(m_teeth_inside))
        {
            m_inside_vert_ind_map.insert({ vd, index });
            index++;
        }
    }
}

void Teeth::initialize_inside_boundary_kdtree()
{
    PointList boundary_verts;

    for (Halfedge_iterator hedge = m_teeth_inside.halfedges_begin(); hedge != m_teeth_inside.halfedges_end(); ++hedge)
    {
        if (hedge->is_border()) // infinite face
        {
            Point source = hedge->opposite()->vertex()->point();
            Vector offset = hedge->vertex()->point() - source;

            for (int i = 0; i < 20; i++)
                boundary_verts.push_back(source + i * 0.05 * offset);
        }
    }

    m_inside_boundary_tree.insert(boundary_verts.begin(), boundary_verts.end());
}

void Teeth::fill_holes_in_vertmap(Polyhedron& poly, VertDoubleMap& vert_map, double isovalue, VertDoubleMap& new_vert_map)
{
    VertexSet finished_set;
    std::vector<VertexList> connected_components;

    for (Vertex_handle vd : vertices(poly))
    {
        new_vert_map.insert({ vd, vert_map[vd] });

        if (finished_set.find(vd) != finished_set.end())
            continue;

        if (vert_map[vd] >= isovalue)
        {
            finished_set.insert(vd);
            continue;
        }

        VertexList component;
        VertexQueue to_finish_queue;
        to_finish_queue.push(vd);

        while (!to_finish_queue.empty())
        {
            Vertex_handle curr_vert = to_finish_queue.front();
            to_finish_queue.pop();

            if (finished_set.find(curr_vert) == finished_set.end())
            {
                finished_set.insert(curr_vert);
                component.push_back(curr_vert);

                Halfedge_vertex_circulator vd_begin = curr_vert->vertex_begin();
                do {
                    Vertex_handle nb_vert = vd_begin->opposite()->vertex();
                    if ((finished_set.find(nb_vert) == finished_set.end()) && vert_map[nb_vert] < isovalue)
                        to_finish_queue.push(nb_vert);
                } while (++vd_begin != curr_vert->vertex_begin());
            }
        }

        if (component.size() > 0)
            connected_components.push_back(component);
    }

    if (connected_components.size() <= 1)
        return;

    std::sort(connected_components.begin(), connected_components.end(), [&](VertexList A, VertexList B) -> bool { return A.size() > B.size(); });

    for (int i = 1; i < connected_components.size(); i++)
    {
        for (int j = 0; j < connected_components[i].size(); j++)
            new_vert_map[connected_components[i][j]] = isovalue + 1e-3;
    }
}

//----------------------------- Convexity Function -----------------------------//

bool Teeth::compute_laplacian_based_convexity(double radius_ratio, double bounding_radius, int smooth_range, double lambda)
{
    reset_convexity();

    // initialize convexity map
    compute_convexity_from_sampling(radius_ratio, bounding_radius);

    // smooth convexity map
    compute_smoothed_field(smooth_range, m_inside_vert_convex_map, m_inside_vert_smooth_convex_map);

    // Initialize linear system
    int nb_variables = size_of_inside_vertices();
    ESMatrix L(nb_variables * 2, nb_variables);
    EVector B(nb_variables * 2), X(nb_variables);

    ESTripleList LT;
    LT.reserve(12 * nb_variables);

    // Assemble system
    IntDoubleMap smoothed_convex_map;
    for (Vertex_handle vd : vertices(m_teeth_inside))
        smoothed_convex_map.insert({ m_inside_vert_ind_map[vd], m_inside_vert_smooth_convex_map[vd] });

    assemble_laplacian_matrix(m_teeth_inside, LT, B, m_inside_vert_ind_map);
    assemble_constraint_matrix(lambda, LT, B, smoothed_convex_map, nb_variables);

    L.setFromTriplets(LT.begin(), LT.end());
    LT.clear();
    bool flag_solver = solve_laplacian(L, B, X);

    if (!flag_solver)
        return false;

    // Assign value
    for (Vertex_handle vd : vertices(m_teeth_inside))
        m_inside_vert_convex_func_map.insert({ vd, X[m_inside_vert_ind_map[vd]] });

    return true;
}

void Teeth::compute_convexity_from_sampling(double radius_ratio, double bounding_radius)
{
    spdlog::trace("Teeth::compute_convexity_from_sampling  radius_ratio:{0} bounding_radius:{1}", 
        radius_ratio, bounding_radius);
    
    PointList samples;
    CGAL::Polygon_mesh_processing::sample_triangle_mesh(m_teeth_inside,
        std::back_inserter(samples),
        CGAL::Polygon_mesh_processing::parameters::use_random_uniform_sampling(false)
        .number_of_points_on_faces(10));
    // std::cout << "  Sampled " << samples.size() << " points on the surface." << std::endl;
    spdlog::trace("Sampled {0} points on the surface.", samples.size());
    KDTree tree(samples.begin(), samples.end());
    double nb_radius = bounding_radius * radius_ratio;

    for (Vertex_handle vd : vertices(m_teeth_inside))
    {
        PointList neighbor_points;
        Fuzzy_circle query(vd->point(), nb_radius);
        tree.search(std::back_inserter(neighbor_points), query);

        if (neighbor_points.size() == 0)
            continue;

        Vector mean_coord = CGAL::NULL_VECTOR;

        for (int i = 0; i < neighbor_points.size(); i++)
        {
            Point neighbor = neighbor_points[i];
            mean_coord += (vd->point() - neighbor);
        }

        mean_coord = mean_coord / (double)neighbor_points.size();
        Vector normal = CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, m_teeth_inside);
        double convexity = mean_coord * normal;
        m_inside_vert_convex_map.insert({ vd, convexity });
    }
}

//----------------------------- Segmentation -----------------------------//

void Teeth::compute_selected_segmentation(PointList& selected_points, double isovalue, Polyhedron& segment_poly)
{
    spdlog::trace("Teeth::compute_selected_segmentation selected_points {0} isovalue {1} segment_poly_size {2}", 
        selected_points.size(), isovalue, segment_poly.size_of_vertices());

    if (m_inside_vert_convex_func_map.empty())
    {
        spdlog::error("No available convexity function map for extracting segmentation.");
        //std::cout << "  No available convexity function map for extracting segmentation." << std::endl;
        return;
    }

    FacetSet finished_set;
    FacetList candidate_facets;
    FacetQueue to_finish_queue;

    for (int i = 0; i < selected_points.size(); i++)
    {
        Point_and_primitive_id pp = m_inside_tree.closest_point_and_primitive(selected_points[i]);
        Point point = pp.first;
        Facet_handle face = pp.second;
        double value = locate_and_evaluate_function(point, face, m_inside_vert_convex_func_map);
        if (value > isovalue)
            to_finish_queue.push(std::make_pair(face, true));
        else
            to_finish_queue.push(std::make_pair(face, false));
    }

    VertDoubleMap vert_value_map;

    while (!to_finish_queue.empty())
    {
        FacetPair elem = to_finish_queue.front();
        to_finish_queue.pop();
        Facet_handle curr_face = elem.first;
        bool flag_pos = elem.second;

        // skip if already processed
        if (finished_set.find(curr_face) != finished_set.end())
            continue;
        finished_set.insert(curr_face);

        // find pos and neg counts
        int pos_count = 0;
        HEdgeList hedges;
        Halfedge_facet_circulator he = curr_face->facet_begin();
        do {
            hedges.push_back(he);
            if (m_inside_vert_convex_func_map[he->vertex()] >= isovalue)
                pos_count++;
        } while (++he != curr_face->facet_begin());

        if (flag_pos && pos_count == 0)
            continue;
        if (!flag_pos && pos_count == 3)
            continue;

        // insert face to the candidate facet list
        candidate_facets.push_back(curr_face);

        // modify the isovalues for negative zones and keep the original isovalues for positive zones
        for (int i = 0; i < 3; i++)
        {
            Vertex_handle vd = hedges[i]->vertex();
            if (!flag_pos)
                vert_value_map[vd] = isovalue + 1.;
            else
            {
                if (vert_value_map.find(vd) == vert_value_map.end())
                    vert_value_map[vd] = m_inside_vert_convex_func_map[vd];
            }
        }

        // insert all neighbor faces into priority queue
        if (pos_count == 3 || pos_count == 0)
        {
            for (int i = 0; i < 3; i++)
            {
                Vertex_handle curr_vert = hedges[i]->vertex();
                Halfedge_vertex_circulator vd_begin = curr_vert->vertex_begin();
                do {
                    if (!vd_begin->is_border())
                    {
                        Facet_handle oppo_face = vd_begin->facet();
                        to_finish_queue.push(std::make_pair(oppo_face, flag_pos));
                    }
                } while (++vd_begin != curr_vert->vertex_begin());
            }
        }
    }

    PointList poly_points;
    extract_inside_vertices_from_function(m_teeth_inside, poly_points);

    HEdgeIntMap he_point_map;
    int count = extract_isovertices_from_function(m_teeth_inside, vert_value_map, isovalue, poly_points, he_point_map);

    PolyList  poly_faces;
    extract_positive_faces_from_function(vert_value_map, isovalue, m_inside_vert_ind_map, he_point_map, candidate_facets, poly_faces);

    segment_poly.clear();
    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, segment_poly);
    //CGAL::Polygon_mesh_processing::keep_largest_connected_components(segment_poly, 1);
    std::cout << "  Compute segment mesh with " << segment_poly.size_of_vertices() << " vertices and " << segment_poly.size_of_facets() << " faces." << std::endl;
}

bool Teeth::validate_selected_segmentation(Polyhedron& segment_poly, bool flag_lower)
{
    for (int i = 0; i < 6; i++)
        m_tooth[i]->reset();

    Polyhedron orig_poly;
    copy_polyehdral_surface(segment_poly, orig_poly);

    std::vector<Polyhedron> split_meshes;
    CGAL::Polygon_mesh_processing::split_connected_components(segment_poly, split_meshes);
    std::cout << "  Found " << split_meshes.size() << " components!" << std::endl;

    if (split_meshes.size() > 6)
    {
        std::cout << "  Found too many components!" << std::endl;
        return false;
    }

    PointList centroids;
    for (int i = 0; i < split_meshes.size(); i++)
    {
        Bbox bounding_box = CGAL::Polygon_mesh_processing::bbox(split_meshes[i]);
        Point center = Point(0.5 * (bounding_box.xmin() + bounding_box.xmax()),
            0.5 * (bounding_box.ymin() + bounding_box.ymax()),
            0.5 * (bounding_box.zmin() + bounding_box.zmax()));
        centroids.push_back(center);
    }

    std::vector<int> indices(centroids.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool { return centroids[A].x() > centroids[B].x(); });

    std::string upper_names[6] = { "23", "22", "21", "11", "12", "13" };
    std::string lower_names[6] = { "33", "32", "31", "41", "42", "43" };

    for (int i = 0; i < indices.size(); i++)
    {
        m_tooth[i]->set_tooth_mesh(split_meshes[indices[i]], m_teeth_inside);

        if(flag_lower)
            m_tooth[i]->set_name(lower_names[i]);
        else
            m_tooth[i]->set_name(upper_names[i]);
    }

    return true;
}

/******************* test **********************/

void Teeth::compute_manu_guide(bool flag_lower, const int start_index, const int end_index, const bool flattened)
{
    m_manu_asm_guide.update(flag_lower,start_index, end_index, flattened);
}

void Teeth::compute_bonding_surface(const std::string& pad, const std::string& padoutline, const std::string& padout)
{
    for (int i = 0; i < 1; i++)
    {
        //if(m_tooth[i]->get_init_flag())
        //{
            m_tooth[i]->compute_bonding_surface(pad, padoutline,padout);
        //}
    }
}

void Teeth::save_validated_segmentation(const std::string& filename)
{
    std::string filename_t = filename + "teeth.stl";
    std::string expanded_filename_t = filename + "teeth_expanded.stl";
    mesh_export(m_teeth_poly, filename_t);

    mesh_export(m_teeth_poly_expanded, expanded_filename_t);

    mesh_export(m_teeth_inside, filename_t + "inside.stl");

    savePointsToXYZ(m_flattened_wire, filename_t + "flattened_wire.xyz");

    for (int i = 0; i < 6; i++)
    {
        if (m_tooth[i]->get_init_flag())
        {
            char buffer[10];
            _itoa(i, buffer, 10);
            std::string filename_i = filename + "tooth" + buffer + ".stl";
            mesh_export(m_tooth[i]->get_tooth_poly(), filename_i);
            std::string filename_i_outline = filename + "pad_outline" + buffer + ".stl";
            mesh_export(m_tooth[i]->get_pad_outline_poly(), filename_i_outline);
            std::string filename_i_pad = filename + "pad" + buffer + ".stl";
            mesh_export(m_tooth[i]->get_pad_shape_poly(), filename_i_pad);
            std::string filename_i_pad_l = filename + "pad_l" + buffer + ".stl";
            mesh_export(m_tooth[i]->get_left_pad_shape_poly(), filename_i_pad_l);
            std::string filename_i_pad_r = filename + "pad_r" + buffer + ".stl";
            mesh_export(m_tooth[i]->get_right_pad_shape_poly(), filename_i_pad_r);
            std::string filename_i_bonding = filename + "pad_bonding" + buffer + ".stl";
            mesh_export(m_tooth[i]->get_bonding_surface(), filename_i_bonding);

            if (i != 0) // i=1,2,3,4,5
            {
                Polyhedron transformed_mesh;

                copy_polyehdral_surface(m_tooth[i]->get_tooth_poly(), transformed_mesh);
                transform_mesh(transformed_mesh, m_wires[i - 1]->getTransformation());
                filename_i = filename + "tooth" + buffer + "_t.stl";
                mesh_export(transformed_mesh, filename_i);

                copy_polyehdral_surface(m_tooth[i]->get_pad_outline_poly(), transformed_mesh);
                transform_mesh(transformed_mesh, m_wires[i - 1]->getTransformation());
                filename_i_outline = filename + "pad_outline" + buffer + "_t.stl";
                mesh_export(transformed_mesh, filename_i_outline);
            }
            if (i == 0)
            {

                filename_i = filename + "tooth" + buffer + "_t.stl";
                mesh_export(m_tooth[i]->get_tooth_poly(), filename_i);
                
                filename_i_outline = filename + "pad_outline" + buffer + "_t.stl";
                mesh_export(m_tooth[i]->get_pad_outline_poly(), filename_i_outline);

                filename_i_pad = filename + "pad" + buffer + "_t.stl";
                mesh_export(m_tooth[i]->get_pad_shape_poly(), filename_i_pad);
            }
            else if (i == 5)
            {
                Polyhedron transformed_mesh;
                copy_polyehdral_surface(m_tooth[i]->get_pad_shape_poly(), transformed_mesh);
                transform_mesh(transformed_mesh, m_wires[i - 1]->getTransformation());

                std::string filename_i_pad = filename + "pad" + buffer + "_t.stl";
                mesh_export(transformed_mesh, filename_i_pad);
            }
            else // i=1,2,3,4
            {
                Polyhedron transformed_mesh;
                copy_polyehdral_surface(m_tooth[i]->get_left_pad_shape_poly(), transformed_mesh);
                transform_mesh(transformed_mesh, m_wires[i - 1]->getTransformation());

                std::string filename_i_pad = filename + "pad_l" + buffer + "_t.stl";
                mesh_export(transformed_mesh, filename_i_pad);

                copy_polyehdral_surface(m_tooth[i]->get_right_pad_shape_poly(), transformed_mesh);
                transform_mesh(transformed_mesh, m_wires[i]->getTransformation());

                filename_i_pad = filename + "pad_r" + buffer + "_t.stl";
                mesh_export(transformed_mesh, filename_i_pad);
            }
        }
    }
    for (int i = 0; i < m_wires.size(); i++)
    {
        char buffer[10];
        _itoa(i, buffer, 10);
        std::string filename_i = filename + "wire" + buffer + ".stl";
        std::string filename_i_t = filename + "wire" + buffer + "_t.stl";
        // std::string filename_il0 = filename + "_wire_shape_left_" + buffer+"_0.stl";
        // std::string filename_il1 = filename + "_wire_shape_left_" + buffer+"_1.stl";
        // std::string filename_il2 = filename + "_wire_shape_left_" + buffer+"_2.stl";
        // std::string filename_ir0 = filename + "_wire_shape_right_" + buffer+"_0.stl";
        // std::string filename_ir1 = filename + "_wire_shape_right_" + buffer+"_1.stl";
        // std::string filename_ir2 = filename + "_wire_shape_right_" + buffer+"_2.stl";
        if (m_wires[i])
        {
            mesh_export(m_wires[i]->get_wire_extrude_shape(), filename_i);
            mesh_export(m_wires[i]->get_wire_extrude_shape(true), filename_i_t);
            // mesh_export(m_wires[i]->get_wire_left_shape(0),filename_il0);
            // mesh_export(m_wires[i]->get_wire_left_shape(1),filename_il1);
            // mesh_export(m_wires[i]->get_wire_left_shape(2),filename_il2);
            // mesh_export(m_wires[i]->get_wire_right_shape(0),filename_ir0);
            // mesh_export(m_wires[i]->get_wire_right_shape(1),filename_ir1);
            // mesh_export(m_wires[i]->get_wire_right_shape(2),filename_ir2);
        }
    }
}

// void Teeth::exportTeethMesh(QJsonObject &meshObject)
// {
//     exportPolyMesh(m_teeth_poly, meshObject);
// }

void Teeth::load_validated_segmentation(const std::string& filename)
{
    
}	

void Teeth::load_guide(const std::string& filename)
{
    
}

void Teeth::get_pad_statistics(StringList& names, DoubleList& left, DoubleList& right)
{
    for (int i = 0; i < 6; i++)
    {
        if(m_tooth[i]->get_init_flag())
            m_tooth[i]->get_pad_statistics(names, left, right);
    }
}

//----------------------------- Tooth Function -----------------------------//

Polyhedron& Teeth::get_tooth_poly(int index) { return m_tooth[index]->get_tooth_poly(); }
VertDoubleMap& Teeth::get_geodesic_distance_map(int index) { return m_tooth[index]->get_geodesic_distance_map(); }
Polyhedron& Teeth::get_pad_outline_poly(int index) { return m_tooth[index]->get_pad_outline_poly(); }

HEdgeSet& Teeth::get_upper_bound_set(int index) { return m_tooth[index]->get_upper_bound_set(); }
PointList& Teeth::get_cut_points(int index) { return m_tooth[index]->get_cut_points(); }
Plane& Teeth::get_cut_plane(int index) { return m_tooth[index]->get_cut_plane(); }

Wire& Teeth::get_wire(int index) { return *m_wires[index]; }

Tooth& Teeth::get_tooth(int index) { return *m_tooth[index]; }

void Teeth::compute_geodesic_pad_outlines(int index, int discretize_step, double isovalue, bool flag_smooth)
{
    if (m_tooth[index]->get_init_flag())
        m_tooth[index]->compute_geodedic_pad_outlines(discretize_step, isovalue, m_teeth_inside, flag_smooth);
}

void Teeth::compute_circular_pad_outlines(int index, int discretize_step, double isovalue, double min_area)
{
    if (m_tooth[index]->get_init_flag())
        m_tooth[index]->compute_circular_pad_outlines(discretize_step, isovalue, min_area);
}

void Teeth::recompute_geodesic_pad_outlines(int index, double isovalue, bool flag_smooth)
{
    if (m_tooth[index]->get_init_flag())
        m_tooth[index]->update_geodesic_pad_outlines(isovalue, m_teeth_inside, flag_smooth);
}

void Teeth::compute_tooth_cut_plane(int index)
{
    if (m_tooth[index]->get_init_flag())
        m_tooth[index]->compute_tooth_cut_plane(m_inside_boundary_tree, m_inside_tree);
}

void Teeth::compute_pad_shapes(int index, bool flag_left, bool flag_right, double pad_height, double wire_height, double wire_radius)
{
    Point teeth_center = CGAL::Polygon_mesh_processing::centroid(m_teeth_inside);	

    if (m_tooth[index]->get_init_flag())
    {
        if (flag_left && !flag_right)
        {
            //SegmentList& wire_outlines = m_wires[index-1]->get_polylines();
            Polyhedron& wire_shape = m_wires[index-1]->get_wire_extrude_shape();
            Polyhedron& enlarge_wire_shape = m_wires[index-1]->get_wire_shape(1);
            m_tooth[index]->compute_pad_shape(wire_shape, enlarge_wire_shape, pad_height, wire_height, teeth_center, wire_radius);
        }
        else if (!flag_left && flag_right)
        {
            //SegmentList& wire_outlines = m_wires[index]->get_polylines();
            Polyhedron& wire_shape = m_wires[index]->get_wire_extrude_shape();
            Polyhedron& enlarge_wire_shape = m_wires[index]->get_wire_shape(1);
            m_tooth[index]->compute_pad_shape(wire_shape, enlarge_wire_shape, pad_height, wire_height, teeth_center, wire_radius);
        }
        else
        {
            //SegmentList& left_wire_outlines = m_wires[index-1]->get_polylines();
            Polyhedron& left_wire_shape = m_wires[index-1]->get_wire_extrude_shape();
            Polyhedron& left_enlarge_wire_shape = m_wires[index-1]->get_wire_shape(1);
            //SegmentList& right_wire_outlines = m_wires[index]->get_polylines();
            Polyhedron& right_wire_shape = m_wires[index]->get_wire_extrude_shape();
            Polyhedron& right_enlarge_wire_shape = m_wires[index]->get_wire_shape(1);
            m_tooth[index]->compute_splited_pad_shape(  left_wire_shape, 
                                                        left_enlarge_wire_shape, 
                                                        right_wire_shape, 
                                                        right_enlarge_wire_shape, 
                                                        pad_height, 
                                                        wire_height, 
                                                        teeth_center,
                                                        wire_radius );
        }
    }
}

//----------------------------- Wire Function -----------------------------//

/*void Teeth::compute_paired_wires(int index, double distance)
{
    if (!m_tooth[index]->get_init_flag() || !m_tooth[index + 1]->get_init_flag())
        return;

    Polyhedron& right_tooth = get_tooth_poly(index);
    Polyhedron& left_tooth = get_tooth_poly(index + 1);

    // Pair fitting plane
    TriangleList candidate_faces;
    copy_faces_to_triangles(left_tooth, candidate_faces);
    Copy_faces_to_triangles(right_tooth, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    Point teeth_centroid = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());

    // Pair fitting line
    PointList left_boundaries, right_boundaries;
    SegmentList segment_boundaries;
    HEdgeSet& right_bound_hedge_set = get_upper_bound_set(index);
    HEdgeSet& left_bound_hedge_set = get_upper_bound_set(index + 1);
    copy_boundary_edges_to_segments(left_bound_hedge_set, left_boundaries, segment_boundaries);
    copy_boundary_edges_to_segments(right_bound_hedge_set, right_boundaries, segment_boundaries);
    Line fitting_line;
    CGAL::linear_least_squares_fitting_3(segment_boundaries.begin(), segment_boundaries.end(), fitting_line, CGAL::Dimension_tag<1>());

    // Pair fitting point
    Point left_left_point, left_right_point;
    project_points_on_line(left_boundaries, fitting_line, left_left_point, left_right_point);
    Point right_left_point, right_right_point;
    project_points_on_line(right_boundaries, fitting_line, right_left_point, right_right_point);
    left_left_point = fitting_line.projection(left_left_point);
    right_right_point = fitting_line.projection(right_right_point);
    Point fitting_point = CGAL::midpoint(left_right_point, right_left_point);

    // Pair cut plane
    Plane pair_plane;
    Vector line_direction = left_left_point-right_right_point;
    compute_pair_plane(distance, fitting_point, line_direction, teeth_centroid, fitting_plane, pair_plane);
    Polyhedron pair_plane_mesh;
    compute_plane_mesh(10., pair_plane, fitting_point, right_right_point, pair_plane_mesh);

    // Compute wire line
    Plane right_plane = get_cut_plane(index);
    Plane left_plane = get_cut_plane(index + 1);
    clip_mesh_by_plane(pair_plane_mesh, left_plane, fitting_point);
    clip_mesh_by_plane(pair_plane_mesh, right_plane, fitting_point);

    std::vector<PointList> polylines;
    CGAL::Polygon_mesh_processing::surface_intersection(m_teeth_inside, pair_plane_mesh, std::back_inserter(polylines));

    // Init Wire class
    m_wires[index]->reset();
    m_wires[index]->set_fitting_plane(fitting_plane, fitting_point, right_right_point);
    Segment fitting_segment(left_left_point, right_right_point);
    m_wires[index]->set_fitting_segment(fitting_segment);
    m_wires[index]->set_fitting_point(fitting_point);
    m_wires[index]->set_cut_plane(pair_plane, fitting_point, right_right_point);

    if (polylines.size() == 0)
    {
        std::cout << "Wrong intersection!" << std::endl;
        m_wires[index]->set_flag(false);
        return;
    }
    else if (polylines.size() > 1)
    {
        size_t max_size = -1;
        int max_id = -1;
        for (int i = 0; i < polylines.size(); i++)
        {
            if (polylines[i].size() > max_size)
            {
                max_size = polylines[i].size();
                max_id = i;
            }
        }
        m_wires[index]->set_polylines(polylines[max_id]);
        m_wires[index]->set_flag(true);
    }
    else
    {
        m_wires[index]->set_polylines(polylines[0]);
        m_wires[index]->set_flag(true);
    }

    return;
}*/

void Teeth::compute_paired_wires(int index, double distance, double tolerance, int max_bend, double radius, double offset)
{
    if (!m_tooth[index]->get_init_flag() || !m_tooth[index + 1]->get_init_flag())
        return;

    Polyhedron& right_tooth = get_tooth_poly(index);
    Polyhedron& left_tooth = get_tooth_poly(index + 1);

    // Pair fitting plane
    TriangleList candidate_faces;
    copy_faces_to_triangles(left_tooth, candidate_faces);
    copy_faces_to_triangles(right_tooth, candidate_faces);
    Plane fitting_plane;
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    Point teeth_centroid = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());

    // Pair fitting line
    PointList left_boundaries, right_boundaries;
    SegmentList segment_boundaries;
    HEdgeSet& right_bound_hedge_set = get_upper_bound_set(index);
    HEdgeSet& left_bound_hedge_set = get_upper_bound_set(index + 1);
    copy_boundary_edges_to_segments(left_bound_hedge_set, left_boundaries, segment_boundaries);
    copy_boundary_edges_to_segments(right_bound_hedge_set, right_boundaries, segment_boundaries);
    Line fitting_line;
    CGAL::linear_least_squares_fitting_3(segment_boundaries.begin(), segment_boundaries.end(), fitting_line, CGAL::Dimension_tag<1>());

    // Pair fitting point
    Point left_left_point, left_right_point;
    project_points_on_line(left_boundaries, fitting_line, left_left_point, left_right_point);
    Point right_left_point, right_right_point;
    project_points_on_line(right_boundaries, fitting_line, right_left_point, right_right_point);
    left_left_point = fitting_line.projection(left_left_point);
    right_right_point = fitting_line.projection(right_right_point);
    Point fitting_point = CGAL::midpoint(left_right_point, right_left_point);

    // Pair cut plane
    Plane pair_plane;
    Vector line_direction = left_left_point-right_right_point;
    compute_pair_plane(distance, fitting_point, line_direction, teeth_centroid, fitting_plane, pair_plane);
    Polyhedron pair_plane_mesh;
    compute_plane_mesh(10., pair_plane, fitting_point, right_right_point, pair_plane_mesh);

    // Compute wire line
    Polyhedron right_halftooth;
    Plane right_plane = get_cut_plane(index);
    copy_polyehdral_surface(right_tooth, right_halftooth);
    if(right_plane.orthogonal_vector() * (fitting_point - right_plane.point()) < 0.)
        CGAL::Polygon_mesh_processing::clip(right_halftooth, right_plane);
    else
        CGAL::Polygon_mesh_processing::clip(right_halftooth, right_plane.opposite());

    // compute polylines
    std::vector<PointList> right_all_polylines;
    PointList right_polylines;
    CGAL::Polygon_mesh_processing::surface_intersection(right_halftooth, pair_plane_mesh, std::back_inserter(right_all_polylines));
    if (right_all_polylines.size() == 0)
    {
        std::cout << "Wrong right intersection!" << std::endl;
        m_wires[index]->set_flag(false);
        return;
    }
    else if (right_all_polylines.size() > 1)
    {
        size_t right_max_size = 0;
        int right_max_id = -1;
        for (int i = 0; i < right_all_polylines.size(); i++)
        {
            if (right_all_polylines[i].size() > right_max_size)
            {
                right_max_size = right_all_polylines[i].size();
                right_max_id = i;
            }
        }
        right_polylines = right_all_polylines[right_max_id];
    }
    else
    {
        right_polylines = right_all_polylines[0];
    }

    Polyhedron left_halftooth;
    Plane left_plane = get_cut_plane(index + 1);
    copy_polyehdral_surface(left_tooth, left_halftooth);
    if(left_plane.orthogonal_vector() * (fitting_point - left_plane.point()) < 0.)
        CGAL::Polygon_mesh_processing::clip(left_halftooth, left_plane);
    else
        CGAL::Polygon_mesh_processing::clip(left_halftooth, left_plane.opposite());
    std::vector<PointList> left_all_polylines;
    PointList left_polylines;
    CGAL::Polygon_mesh_processing::surface_intersection(left_halftooth, pair_plane_mesh, std::back_inserter(left_all_polylines));
    
    if (left_all_polylines.size() == 0)
    {
        std::cout << "Wrong right intersection!" << std::endl;
        m_wires[index]->set_flag(false);
        return;
    }
    else if (left_all_polylines.size() > 1)
    {
        size_t left_max_size = 0;
        int left_max_id = -1;
        for (int i = 0; i < left_all_polylines.size(); i++)
        {
            if (left_all_polylines[i].size() > left_max_size)
            {
                left_max_size = left_all_polylines[i].size();
                left_max_id = i;
            }
        }
        left_polylines.insert(left_polylines.end(), left_all_polylines[left_max_id].begin(), left_all_polylines[left_max_id].end());
    }
    else
    {
        left_polylines = left_all_polylines[0];
    }

    // merge polylines
    PointList polylines;
    size_t size_right_poly = right_polylines.size();
    size_t size_left_poly = left_polylines.size();
    bool flag_left_inverse = (CGAL::squared_distance(left_polylines[0], right_polylines[0]) < CGAL::squared_distance(left_polylines[size_left_poly-1], right_polylines[0]));
    bool flag_right_inverse = (CGAL::squared_distance(left_polylines[0], right_polylines[0]) > CGAL::squared_distance(left_polylines[0], right_polylines[size_right_poly-1]));

    for(int i = 0; i < size_left_poly; i++)
    {
        if(flag_left_inverse)
            polylines.push_back(left_polylines[size_left_poly-i-1]);
        else
            polylines.push_back(left_polylines[i]);
    }

    for(int i = 0; i < size_right_poly; i++)
    {
        if(flag_right_inverse)
            polylines.push_back(right_polylines[size_right_poly-i-1]);
        else
            polylines.push_back(right_polylines[i]);
    }

    // Init Wire class
    m_wires[index]->reset();
    m_wires[index]->set_fitting_plane(fitting_plane, fitting_point, right_right_point);
    Segment fitting_segment(left_left_point, right_right_point);
    m_wires[index]->set_fitting_segment(fitting_segment);
    m_wires[index]->set_fitting_point(fitting_point);
    m_wires[index]->set_cut_plane(pair_plane, fitting_point, right_right_point);

    // Smooth and offset wire
    std::cout << "Smooth tolerance: " << tolerance << std::endl;
    compute_smooth_wire(tolerance, max_bend, polylines, pair_plane);
    Point teeth_center = CGAL::Polygon_mesh_processing::centroid(m_teeth_inside);
    Vector offset_corrector = teeth_center - candidate_faces[0][0];
    Vector polyline_direction = polylines[polylines.size()-1] - polylines[0];
    Vector offset_direction = CGAL::cross_product(polyline_direction, pair_plane.orthogonal_vector());
    radius = radius * offset;
    offset_direction = offset_direction / std::sqrt(offset_direction.squared_length());
    offset_corrector = offset_corrector / std::sqrt(offset_corrector.squared_length());
    if(offset_direction * offset_corrector < 0)
        offset_direction = -offset_direction;
    offset_direction = offset_direction * radius;
    for(int i = 0; i < polylines.size(); i++)
        polylines[i] = polylines[i] + offset_direction;
    
    // Set polyline
    m_wires[index]->set_polylines(polylines);
    m_wires[index]->set_flag(true);
    
    return;
}

void Teeth::compute_smooth_wire(double tolerance, int max_bend, PointList& polylines, Plane& wire_plane)
{
    Point center = polylines[0];
    Vector base_1 = wire_plane.base1();
    base_1 = base_1 / std::sqrt(base_1.squared_length());
    Vector base_2 = wire_plane.base2();
    base_2 = base_2 / std::sqrt(base_2.squared_length());

    Point2dList proj_polyline_pts;
    for(int i = 0; i < polylines.size(); i++)
    {
        Point_2 proj = to_2d(wire_plane, center, base_1, base_2, polylines[i]);
        proj_polyline_pts.push_back(proj);
    }

    std::vector<Point2dList> multi_proj_polyline_pts;
    multi_proj_polyline_pts.push_back(proj_polyline_pts);
    std::vector<Line_2> multi_fitting_lines;
    Line_2 first_fitting_line;
    double distortion = compute_fitting_line_2d(proj_polyline_pts, first_fitting_line);
    std::cout << "Distortion: " << distortion << std::endl;
    multi_fitting_lines.push_back(first_fitting_line);
    DoubleList multi_distortions;
    multi_distortions.push_back(distortion);

    bool flag = false;
    while(!flag)
    {
        flag = split_fitting_line_with_criterion(multi_proj_polyline_pts, multi_fitting_lines, multi_distortions, tolerance, max_bend);
    }

    Point2dList polylines_2d;
    compute_smooth_intersection(multi_proj_polyline_pts, multi_fitting_lines, polylines_2d);

    polylines.clear();
    for(int i = 0; i < polylines_2d.size(); i++)
    {
        Point proj = to_3d(center, base_1, base_2, polylines_2d[i]);
        polylines.push_back(proj);
    }
}

void Teeth::compute_smooth_intersection(std::vector<Point2dList>& multi_proj_polyline_pts, std::vector<Line_2>& multi_fitting_lines, Point2dList& polylines_2d)
{
    // Start point
    Point_2 start_proj = multi_fitting_lines[0].projection(multi_proj_polyline_pts[0][0]);
    polylines_2d.push_back(start_proj);

    // Compute intersection
    for(int i = 0; i < multi_fitting_lines.size(); i++)
    {
        Point_2 proj = multi_fitting_lines[i].projection(multi_proj_polyline_pts[i].back());
        polylines_2d.push_back(proj);
    }

    std::cout << "Fitting line size: " << multi_fitting_lines.size() << std::endl;
    std::cout << "Polyline size: " << polylines_2d.size() << std::endl;
}

bool Teeth::split_fitting_line_with_criterion(std::vector<Point2dList>& multi_proj_polyline_pts, std::vector<Line_2>& multi_fitting_lines, DoubleList& multi_distortions, double tol, int max_bend)
{
    auto min_distortion_ind = std::min_element(multi_distortions.begin(), multi_distortions.end());
    double min_distortion = *min_distortion_ind;
    if(min_distortion > tol || multi_distortions.size() >= max_bend)
    {
        std::cout << "Criterion satisfied!" << std::endl;
        return true;
    }

    std::cout << "Min quality: " << min_distortion << std::endl;

    int split_ind = min_distortion_ind - multi_distortions.begin();
    std::cout << "  split wire index: " << split_ind << std::endl;
    int min_split_ind = add_split_point_for_polyline(multi_proj_polyline_pts[split_ind], multi_fitting_lines[split_ind]);
    std::cout << "  split point index: " << min_split_ind << std::endl;

    std::vector<Point2dList> new_multi_proj_polyline_pts;
    std::vector<Line_2> new_multi_fitting_lines;
    DoubleList new_multi_distortions;

    for(int i = 0; i < multi_proj_polyline_pts.size(); i++)
    {
        if(i < split_ind || i > split_ind)
        {
            new_multi_proj_polyline_pts.push_back(multi_proj_polyline_pts[i]);
            new_multi_fitting_lines.push_back(multi_fitting_lines[i]);
            new_multi_distortions.push_back(multi_distortions[i]);
        }
        else
        {
            Point2dList leftpart, rightpart;
            for(int j = 0; j < multi_proj_polyline_pts[i].size(); j++)
            {
                if(j <= min_split_ind)
                    leftpart.push_back(multi_proj_polyline_pts[i][j]);
                if(j >= min_split_ind)
                    rightpart.push_back(multi_proj_polyline_pts[i][j]);
            }

            if(leftpart.size() > 0)
            {
                Line_2 left_line;
                double left_distortion = compute_fitting_line_2d(leftpart, left_line);
                new_multi_proj_polyline_pts.push_back(leftpart);
                new_multi_fitting_lines.push_back(left_line);
                new_multi_distortions.push_back(left_distortion);
            }

            if(rightpart.size() > 0)
            {
                Line_2 right_line;
                double right_distortion = compute_fitting_line_2d(rightpart, right_line);
                new_multi_proj_polyline_pts.push_back(rightpart);
                new_multi_fitting_lines.push_back(right_line);
                new_multi_distortions.push_back(right_distortion);
            }
        }
    }

    multi_proj_polyline_pts.swap(new_multi_proj_polyline_pts);
    new_multi_proj_polyline_pts.clear();
    multi_fitting_lines.swap(new_multi_fitting_lines);
    new_multi_fitting_lines.clear();
    multi_distortions.swap(new_multi_distortions);
    new_multi_distortions.clear();

    return false;
}

double Teeth::compute_fitting_line_2d(Point2dList& proj_seg_points, Line_2& proj_line)
{
    Segment2dList proj_polylines;
    for(int i = 0; i < proj_seg_points.size()-1; i++)
    {
        Segment_2 seg(proj_seg_points[i], proj_seg_points[i+1]);
        proj_polylines.push_back(seg);
    }

    double distortion = CGAL::linear_least_squares_fitting_2(proj_polylines.begin(), proj_polylines.end(), proj_line, CGAL::Dimension_tag<1>());
    return distortion;
}

int Teeth::add_split_point_for_polyline(Point2dList& proj_seg_points, Line_2& proj_line)
{
    double max_quality = 0.;
    int max_ind = 0;

    Point2dList leftpart, rightpart;
    Line_2 left_fitting, right_fitting;

    int min_ind = 5;
    

    for(int i = 0; i < proj_seg_points.size(); i++)
    {
        if(i <= min_ind)
            leftpart.push_back(proj_seg_points[i]);
        if(i >= min_ind)
            rightpart.push_back(proj_seg_points[i]);
    }

    for(int i = min_ind; i < proj_seg_points.size()-min_ind; i++)
    {
        double left_distortion = compute_fitting_line_2d(leftpart, left_fitting);
        double right_distortion = compute_fitting_line_2d(rightpart, right_fitting);
        double quality = std::min(left_distortion, right_distortion);

        if(quality > max_quality)
        {
            max_quality = quality;
            max_ind = i;
        }

        rightpart.erase(rightpart.begin());
        leftpart.push_back(rightpart[0]);
    }

    return max_ind;
}

void Teeth::compute_coplanar_wires(double distance)
{
    for(int index = 0; index < m_tooth.size(); index++)
    {
        if(!m_tooth[index]->get_init_flag())
        {
            std::cout << "Tooth not initialized!" << std::endl;
            return;
        }
    }
        
    // Fitting plane for all teeth
    TriangleList candidate_faces;
    for(int index = 0; index < m_tooth.size(); index++)
        copy_faces_to_triangles(get_tooth_poly(index), candidate_faces);
    //Plane fitting_plane;
    //CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), candidate_faces.end(), fitting_plane, CGAL::Dimension_tag<2>());
    Point plane_centroid = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());
    Plane fitting_plane(plane_centroid, Vector(0., 1., 0.));

    // Fitting line for all teeth
    SegmentList segment_boundaries;
    PointList point_boundaries;
    for(int index = 0; index < m_tooth.size(); index++)
    {
        HEdgeSet& bound_hedge_set = get_upper_bound_set(index);
        copy_boundary_edges_to_segments(bound_hedge_set, point_boundaries, segment_boundaries);
    }
    Line fitting_line;
    CGAL::linear_least_squares_fitting_3(segment_boundaries.begin(), segment_boundaries.end(), fitting_line, CGAL::Dimension_tag<1>());
    Point fitting_point = CGAL::centroid(point_boundaries.begin(), point_boundaries.end());

    // Coplanar cut plane
    Plane coplane;
    assert(0);
    // compute_pair_plane(distance, fitting_point, fitting_line, plane_centroid, fitting_plane, coplane);
    Polyhedron coplane_mesh;
    compute_plane_mesh(100., coplane, fitting_point, point_boundaries[0], coplane_mesh);

    for(int index = 0; index < m_tooth.size() - 1; index++)
    {
        std::cout << "    Begin computing wire " << index << std::endl;
        Polyhedron current_plane;
        copy_polyehdral_surface(coplane_mesh, current_plane);

        // Compute wire line
        Plane right_plane = get_cut_plane(index);
        Plane left_plane = get_cut_plane(index + 1);

        Point right_proj_point, left_proj_point;
        if(!compute_intersection_plane_line(right_plane, fitting_line, right_proj_point))
            right_proj_point = right_plane.projection(fitting_point);
        if(!compute_intersection_plane_line(left_plane, fitting_line, left_proj_point))
            left_proj_point = left_plane.projection(fitting_point);

        Point mid = CGAL::midpoint(right_proj_point, left_proj_point);
        clip_mesh_by_plane(current_plane, left_plane, mid);
        clip_mesh_by_plane(current_plane, right_plane, mid); 

        std::vector<PointList> polylines;
        CGAL::Polygon_mesh_processing::surface_intersection(m_teeth_inside, current_plane, std::back_inserter(polylines));      

        // Init Wire class
        m_wires[index]->reset();
        m_wires[index]->set_fitting_plane(fitting_plane, mid, right_proj_point);
        Segment fitting_segment(fitting_line.projection(left_proj_point), fitting_line.projection(right_proj_point));
        m_wires[index]->set_fitting_segment(fitting_segment);
        m_wires[index]->set_fitting_point(mid);
        m_wires[index]->set_cut_plane(coplane, mid, right_proj_point);

        if (polylines.size() == 0)
        {
            std::cout << "Wrong intersection!" << std::endl;
            m_wires[index]->set_flag(false);
        }
        else if (polylines.size() > 1)
        {
            size_t max_size = -1;
            int max_id = -1;
            for (int i = 0; i < polylines.size(); i++)
            {
                if (polylines[i].size() > max_size)
                {
                    max_size = polylines[i].size();
                    max_id = i;
                }
            }
            m_wires[index]->set_polylines(polylines[max_id]);
            m_wires[index]->set_flag(true);
        }
        else
        {
            m_wires[index]->set_polylines(polylines[0]);
            m_wires[index]->set_flag(true);
        } 
    }
}

void Teeth::copy_boundary_edges_to_segments(HEdgeSet& bound_hedge_set, PointList& points, SegmentList& segments)
{
    for (auto& hedge : bound_hedge_set)
    {
        points.push_back(hedge->vertex()->point());
        points.push_back(hedge->opposite()->vertex()->point());
        segments.push_back(Segment(hedge->vertex()->point(), hedge->opposite()->vertex()->point()));
    }
}

bool Teeth::compute_intersection_plane_line(Plane& plane, Line& line, Point& point)
{
    const auto result = CGAL::intersection(plane, line);

    if(result) 
    {
        if (const Point* s = boost::get<Point>(&*result)) 
        {
            point = Point(s->x(), s->y(), s->z());
            return true;
        } 
        else 
            return false;
    }

    return false;
}

void Teeth::project_points_on_line(PointList& candidate_points, Line& fitting_line, Point& left_point, Point& right_point)
{
    Point center = fitting_line.projection(candidate_points[0]);
    Vector direction = fitting_line.to_vector();
    Vector xaxis(1., 0., 0.);
    if (direction * xaxis < 0.)
        direction = -direction;

    double fmin = 0;
    double fmax = 0;
    int min_index = 0;
    int max_index = 0;

    for (int i = 0; i < candidate_points.size(); i++)
    {
        Point proj = fitting_line.projection(candidate_points[i]);
        double value = CGAL::squared_distance(proj, center);
        if ((proj - center) * direction < 0)
            value = -value;

        if (value < fmin)
        {
            fmin = value;
            min_index = i;
        }
        if (value > fmax)
        {
            fmax = value;
            max_index = i;
        }
    }

    left_point = candidate_points[min_index];
    right_point = candidate_points[max_index];
}

void Teeth::compute_pair_plane(double distance, Point& plane_center, Vector& line_direction, Point& teeth_center, Plane& fitting_plane, Plane& pair_plane)
{
    Vector vertical = fitting_plane.orthogonal_vector();

    Vector normal = CGAL::cross_product(line_direction, vertical);
    normal = normal / std::sqrt(normal.squared_length());
    if (normal * (teeth_center - plane_center) < 0.)
        normal = -normal;

    Point moved_center = plane_center + distance * normal;
    pair_plane = Plane(moved_center, normal);
}



/*
void Teeth::compute_wire_shapes(int index, int discretize_step, double width)
{
    if (m_teeth_inside.size_of_vertices() == 0 || !m_wires[index]->get_flag())
        return;

    Polyhedron& right_tooth = get_pad_outline_poly(index);
    Polyhedron& left_tooth = get_pad_outline_poly(index + 1);
    Plane right_plane = get_cut_plane(index);
    Plane left_plane = get_cut_plane(index + 1);
    Point fitting_point = m_wires[index]->get_fitting_point();

    Polyhedron inter_teeth;
    copy_polyehdral_surface(m_teeth_inside, inter_teeth);
    clip_mesh_by_plane(inter_teeth, left_plane, fitting_point);
    clip_mesh_by_plane(inter_teeth, right_plane, fitting_point);

    Geodesic_tree geodesic_tree(inter_teeth,
        get(boost::vertex_external_index, inter_teeth),
        get(CGAL::halfedge_external_index, inter_teeth),
        get(CGAL::face_external_index, inter_teeth),
        get(CGAL::vertex_point, inter_teeth));
    AABBTree geodesic_aabb_tree;
    geodesic_tree.build_aabb_tree(geodesic_aabb_tree);
    geodesic_aabb_tree.accelerate_distance_queries();

    // discretize segments
    double offset = 1. / discretize_step;
    std::vector<Face_location> faceLocations;

    SegmentList wire_outlines = m_wires[index]->get_polylines();
    for (int j = 0; j < wire_outlines.size(); j++)
    {
        Segment& seg = wire_outlines[j];
        Point source = seg.source();
        Vector direction = seg.target() - seg.source();
        for (int i = 0; i < discretize_step; i++)  // discretize edge
        {
            Point p = source + i * offset * direction;
            Face_location source_loc = geodesic_tree.locate<AABBTraits>(p, geodesic_aabb_tree);
            if (source_loc.first != boost::graph_traits<Polyhedron>::null_face())
                faceLocations.push_back(source_loc);
        }
    }
    geodesic_tree.add_source_points(faceLocations.begin(), faceLocations.end());

    Polyhedron wire_input, wire_shape;
    copy_polyehdral_surface(inter_teeth, wire_input);
    bool flag_wire = extract_wire_shape_for_halftooth(width, geodesic_tree, geodesic_aabb_tree, wire_input, wire_shape);

    if (!flag_wire)
        return;

    // restore the left and right half wire
    Polyhedron left_halftooth, right_halftooth;
    copy_polyehdral_surface(right_tooth, right_halftooth);
    clip_mesh_by_plane(right_halftooth, right_plane, fitting_point);
    copy_polyehdral_surface(left_tooth, left_halftooth);
    clip_mesh_by_plane(left_halftooth, left_plane, fitting_point);
    Polyhedron left_halfwire, right_halfwire;
    Point left_upper_point, right_upper_point, left_lower_point, right_lower_point;
    bool flag_right = extract_wire_shape_for_halftooth(width, right_plane, geodesic_tree, geodesic_aabb_tree,
                                                       right_halftooth, right_halfwire, right_upper_point, right_lower_point);
    bool flag_left = extract_wire_shape_for_halftooth(width, left_plane, geodesic_tree, geodesic_aabb_tree,
                                                      left_halftooth, left_halfwire, left_upper_point, left_lower_point);

    if (!flag_right || !flag_left)
        return;
    m_wires[index]->set_wire_left_shape(left_halfwire);
    m_wires[index]->set_wire_right_shape(right_halfwire);

    m_wires[index]->set_wire_shape(wire_shape);

    // Extrude
    Polyhedron extrude_mesh, wrap_mesh;
    compute_extrusion(wire_shape, extrude_mesh, width * 2.);
    alpha_wrap_mesh(extrude_mesh, wrap_mesh);

    m_wires[index]->set_wire_extrude_shape(wrap_mesh);

    std::cout << "compute_wire_shapes Finish! width:" << width << std::endl;
}*/

/*void Teeth::compute_wire_shapes(int index, int discretize_step, double width)
{
    if (m_teeth_inside.size_of_vertices() == 0 || !m_wires[index]->get_flag())
        return;

    Polyhedron& right_tooth = get_pad_outline_poly(index);
    Polyhedron& left_tooth = get_pad_outline_poly(index + 1);
    Plane right_plane = get_cut_plane(index);
    Plane left_plane = get_cut_plane(index + 1);
    Point fitting_point = m_wires[index]->get_fitting_point();

    Polyhedron inter_teeth;
    copy_polyehdral_surface(m_teeth_inside, inter_teeth);
    clip_mesh_by_plane(inter_teeth, left_plane, fitting_point);
    clip_mesh_by_plane(inter_teeth, right_plane, fitting_point);

    Geodesic_tree geodesic_tree(inter_teeth,
        get(boost::vertex_external_index, inter_teeth),
        get(CGAL::halfedge_external_index, inter_teeth),
        get(CGAL::face_external_index, inter_teeth),
        get(CGAL::vertex_point, inter_teeth));
    AABBTree geodesic_aabb_tree;
    geodesic_tree.build_aabb_tree(geodesic_aabb_tree);
    geodesic_aabb_tree.accelerate_distance_queries();

    // discretize segments
    double offset = 1. / discretize_step;
    std::vector<Face_location> faceLocations;

    SegmentList wire_outlines = m_wires[index]->get_polylines();
    for (int j = 0; j < wire_outlines.size(); j++)
    {
        Segment& seg = wire_outlines[j];
        Point source = seg.source();
        Vector direction = seg.target() - seg.source();
        for (int i = 0; i < discretize_step; i++)  // discretize edge
        {
            Point p = source + i * offset * direction;
            Face_location source_loc = geodesic_tree.locate<AABBTraits>(p, geodesic_aabb_tree);
            if (source_loc.first != boost::graph_traits<Polyhedron>::null_face())
                faceLocations.push_back(source_loc);
        }
    }
    geodesic_tree.add_source_points(faceLocations.begin(), faceLocations.end());

    Polyhedron left_halftooth, right_halftooth;
    copy_polyehdral_surface(right_tooth, right_halftooth);
    clip_mesh_by_plane(right_halftooth, right_plane, fitting_point);
    copy_polyehdral_surface(left_tooth, left_halftooth);
    clip_mesh_by_plane(left_halftooth, left_plane, fitting_point);

    Polyhedron left_halfwire, right_halfwire;
    Point left_upper_point, right_upper_point, left_lower_point, right_lower_point;
    
    bool flag_right = false; bool flag_left = false;


    flag_right = extract_wire_shape_for_halftooth(width, right_plane, geodesic_tree, geodesic_aabb_tree, 
    right_halftooth, right_halfwire, right_upper_point, right_lower_point);
    flag_left = extract_wire_shape_for_halftooth(width, left_plane, geodesic_tree, geodesic_aabb_tree, 
    left_halftooth, left_halfwire, left_upper_point, left_lower_point);

    if (!flag_right || !flag_left)
        return;

    Polyhedron wire_shape;

    Polyhedron middle_wire;
    PointList bound_points;
    bound_points.push_back(right_upper_point);
    bound_points.push_back(right_lower_point);
    bound_points.push_back(left_lower_point);
    bound_points.push_back(left_upper_point);
    PolyList bound_polys;
    bound_polys.push_back(IntList({ 0, 1, 2 }));
    bound_polys.push_back(IntList({ 0, 2, 3 }));
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(bound_points, bound_polys, middle_wire);
    Polyhedron middle_result;
    merge_polyehdral_surface(left_halfwire, right_halfwire, middle_result);
    merge_polyehdral_surface(middle_result, middle_wire, wire_shape);

    m_wires[index]->set_wire_left_shape(left_halfwire);
    m_wires[index]->set_wire_middle_shape(middle_wire);
    m_wires[index]->set_wire_right_shape(right_halfwire);
    m_wires[index]->set_wire_shape(wire_shape);

    // Extrude
    Polyhedron extrude_mesh, extrude_leftwire, extrude_rightwire, extrude_middle;

    compute_extrusion(left_halfwire, extrude_leftwire, width * 2.);
    compute_extrusion(right_halfwire, extrude_rightwire, width * 2.);
    compute_extrusion(middle_wire, extrude_middle, width * 2.);

    CGAL::Polygon_mesh_processing::corefine_and_compute_union(extrude_leftwire, extrude_middle, extrude_mesh);
    CGAL::Polygon_mesh_processing::corefine_and_compute_union(extrude_mesh, extrude_rightwire, extrude_mesh);


    m_wires[index]->set_wire_extrude_shape(extrude_mesh);

    // expanded wire shape expansion ration
    float expansion_ratio = get_variable_float(WIRE_EXPANSION_PAD_RATIO, 2.0);

    flag_right = extract_wire_shape_for_halftooth(width*expansion_ratio, right_plane, geodesic_tree, geodesic_aabb_tree, 
    right_halftooth, right_halfwire, right_upper_point, right_lower_point);
    flag_left = extract_wire_shape_for_halftooth(width*expansion_ratio, left_plane, geodesic_tree, geodesic_aabb_tree, 
    left_halftooth, left_halfwire, left_upper_point, left_lower_point);
    // if (!flag_right || !flag_left)
    //     return;
    m_wires[index]->set_wire_left_shape(left_halfwire, 1);
    m_wires[index]->set_wire_right_shape(right_halfwire, 1);

    flag_right = extract_wire_shape_for_halftooth(1, right_plane, geodesic_tree, geodesic_aabb_tree, 
    right_halftooth, right_halfwire, right_upper_point, right_lower_point);
    flag_left = extract_wire_shape_for_halftooth(1, left_plane, geodesic_tree, geodesic_aabb_tree, 
    left_halftooth, left_halfwire, left_upper_point, left_lower_point);
    // if (!flag_right || !flag_left)
    //     return;
    m_wires[index]->set_wire_left_shape(left_halfwire, 2);
    m_wires[index]->set_wire_right_shape(right_halfwire, 2);
    std::cout << "Finish!" << std::endl;
}*/

void Teeth::compute_wire_shapes(int index, int discretize_step, double radius)
{
    std::cout << "Compute wire shape: " << index << std::endl;

    if (m_teeth_inside.size_of_vertices() == 0 || !m_wires[index]->get_flag())
        return;

    SegmentList wire_outlines = m_wires[index]->get_polylines();
    SegmentTree wire_tree(wire_outlines.begin(), wire_outlines.end());
    auto cylinder = [&](Point p) { return (std::sqrt(wire_tree.squared_distance(p)) - radius); };

    Point bounding_center = CGAL::midpoint(wire_outlines[0].source(), wire_outlines[wire_outlines.size()-1].target());
    double bounding_radius = std::sqrt(CGAL::squared_distance(bounding_center, wire_outlines[0].source()));

    Tr tr;            // 3D-Delaunay triangulation
    C2t3 c2t3(tr);   // 2D-complex in 3D-Delaunay triangulation
    // defining the surface
    Surface_3 surface(cylinder,             // pointer to function
                      Sphere(bounding_center, bounding_radius + 20.)); // bounding sphere
    // Note that "2." above is the *squared* radius of the bounding sphere!
    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<Tr> criteria( 30.,  // angular bound
                                                        0.1,  // radius bound
                                                        0.1); // distance bound
    // meshing surface
    CGAL::make_surface_mesh(c2t3, surface, criteria, CGAL::Non_manifold_tag());
    Polyhedron wire_shape;
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, wire_shape);
    
    // Cut left and right side
    Plane right_plane = get_cut_plane(index);
    Plane left_plane = get_cut_plane(index + 1);
    Point fitting_point = m_wires[index]->get_fitting_point();
    clip_volume_by_plane(wire_shape, left_plane, fitting_point);
    clip_volume_by_plane(wire_shape, right_plane, fitting_point);

    m_wires[index]->set_wire_extrude_shape(wire_shape);
    std::cout << "Finish for " << index << std::endl;
}

void Teeth::compute_enlarge_wire_shapes(int index, int discretize_step, double radius, double enlarge_ratio)
{
    std::cout << "Compute wire shape: " << index << "enlarge ration:" << enlarge_ratio << "radius:" << radius << std::endl;

    
    if (!m_wires[index]->get_flag())
        return;

    Polyhedron& wire_shape = m_wires[index]->get_wire_extrude_shape();
    if (wire_shape.size_of_vertices() == 0)
        return;

    Polyhedron enlarge_wire;
    copy_polyehdral_surface(wire_shape, enlarge_wire);

    VertVecMap vert_normals;
    double offset = radius * (enlarge_ratio - 1.);

    for(Vertex_handle vd: vertices(enlarge_wire))
    {
        Vector normal = CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, enlarge_wire);
        vert_normals.insert({vd, normal});
    }
    
    for(Vertex_handle vd: vertices(enlarge_wire))
    {
        Point p = vd->point() + vert_normals[vd] * offset;
        vd->point() = p;
    }

    m_wires[index]->set_wire_shape(enlarge_wire, 1);
    std::cout << "Finish for " << index << std::endl;
}

void Teeth::clip_mesh_by_plane(Polyhedron& mesh, Plane& plane, Point& center)
{
    if ((center - plane.point()) * plane.orthogonal_vector() > 0.)
        CGAL::Polygon_mesh_processing::clip(mesh, plane.opposite());
    else
        CGAL::Polygon_mesh_processing::clip(mesh, plane);
}

void Teeth::clip_volume_by_plane(Polyhedron& mesh, Plane& plane, Point& center)
{
    if ((center - plane.point()) * plane.orthogonal_vector() > 0.)
        CGAL::Polygon_mesh_processing::clip(mesh, plane.opposite(), CGAL::parameters::clip_volume(true));
    else
        CGAL::Polygon_mesh_processing::clip(mesh, plane, CGAL::parameters::clip_volume(true));
}

bool Teeth::extract_wire_shape_for_halftooth(double width, Geodesic_tree& geodesic_tree, AABBTree& geodesic_aabb_tree, Polyhedron& halftooth,  // input
    Polyhedron& halfwire)   // output
{
    // isotropic remeshing
    std::cout << "  Before isotropic remeshing: " << halftooth.size_of_vertices() << " vertices and " << halftooth.size_of_facets() << std::endl;
    CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(halftooth), 0.8 * width, halftooth);
    std::cout << "  After isotropic remeshing: " << halftooth.size_of_vertices() << " vertices and " << halftooth.size_of_facets() << std::endl;
    std::cout << "  extract_wire_shape_for_halftooth, wire width: " << width << std::endl;

    VertDoubleMap vert_geodesic_map;
    VertIntMap vert_ind_map;
    int count = 0;
    for (Vertex_handle vd : vertices(halftooth))
    {
        Point proj = geodesic_aabb_tree.closest_point(vd->point());
        Face_location proj_loc = geodesic_tree.locate<AABBTraits>(proj, geodesic_aabb_tree);
        double dist = geodesic_tree.shortest_distance_to_source_points(proj_loc.first, proj_loc.second).first;
        vert_geodesic_map.insert({ vd, -dist });
        vert_ind_map.insert({ vd, count });
        count++;
    }

    // Extract surface
    PointList poly_points;
    extract_inside_vertices_from_function(halftooth, poly_points);

    HEdgeIntMap he_point_map;
    extract_isovertices_from_function(halftooth, vert_geodesic_map, -width, poly_points, he_point_map);

    FacetList candidate_facets;
    for (Facet_handle fd : faces(halftooth))
        candidate_facets.push_back(fd);

    PolyList  poly_faces;
    extract_positive_faces_from_function(vert_geodesic_map, -width, vert_ind_map, he_point_map, candidate_facets, poly_faces);
    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, halfwire);

    return true;
}

bool Teeth::extract_wire_shape_for_halftooth(double width, Plane& cut_plane, Geodesic_tree& geodesic_tree, AABBTree& geodesic_aabb_tree, Polyhedron& halftooth,  // input
    Polyhedron& halfwire, Point& upper_bound, Point& lower_bound)   // output
{
    // isotropic remeshing
    std::cout << "  Before isotropic remeshing: " << halftooth.size_of_vertices() << " vertices and " << halftooth.size_of_facets() << std::endl;
    CGAL::Polygon_mesh_processing::isotropic_remeshing(faces(halftooth), 0.8 * width, halftooth);
    std::cout << "  After isotropic remeshing: " << halftooth.size_of_vertices() << " vertices and " << halftooth.size_of_facets() << std::endl;

    VertDoubleMap vert_geodesic_map;
    VertIntMap vert_ind_map;
    int count = 0;
    for (Vertex_handle vd : vertices(halftooth))
    {
        Point proj = geodesic_aabb_tree.closest_point(vd->point());
        Face_location proj_loc = geodesic_tree.locate<AABBTraits>(proj, geodesic_aabb_tree);
        double dist = geodesic_tree.shortest_distance_to_source_points(proj_loc.first, proj_loc.second).first;
        if (dist < 0.)
            vert_geodesic_map.insert({ vd, -width - 1. });
        else
            vert_geodesic_map.insert({ vd, -dist });
        vert_ind_map.insert({ vd, count });
        count++;
    }

    // Extract surface
    PointList poly_points;
    extract_inside_vertices_from_function(halftooth, poly_points);

    HEdgeIntMap he_point_map;
    extract_isovertices_from_function(halftooth, vert_geodesic_map, -width, poly_points, he_point_map);

    PointList bounds;
    for (Halfedge_iterator hedge = halftooth.halfedges_begin(); hedge != halftooth.halfedges_end(); ++hedge)
    {
        if (hedge->is_border_edge() && he_point_map.find(hedge) != he_point_map.end())
        {
            Point p = poly_points[he_point_map[hedge]];
            if (CGAL::squared_distance(p, cut_plane.projection(p)) > 0.5)
                bounds.push_back(p);
        }
    }

    FacetList candidate_facets;
    for (Facet_handle fd : faces(halftooth))
        candidate_facets.push_back(fd);

    PolyList  poly_faces;
    extract_positive_faces_from_function(vert_geodesic_map, -width, vert_ind_map, he_point_map, candidate_facets, poly_faces);
    CGAL::Polygon_mesh_processing::remove_isolated_points_in_polygon_soup(poly_points, poly_faces);
    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, halfwire);

    if (bounds.size() < 2)
    {
        std::cout << "Wrong boundary points!" << std::endl;
        return false;
    }

    std::sort(bounds.begin(), bounds.end(), [&](Point A, Point B) -> bool { return A.z() > B.z(); });
    upper_bound = bounds[0];
    lower_bound = bounds[bounds.size()-1];
    return true;
}

void Teeth::find_vnormal_from_adjacent_mesh(Polyhedron& halfwire, VertexList& vert_handles, VectorList& vert_normals, double threshold)
{
    for (Vertex_handle vd : vertices(halfwire))
    {
        for (int i = 0; i < vert_handles.size(); i++)
        {
            if (CGAL::squared_distance(vd->point(), vert_handles[i]->point()) < threshold)
                vert_normals[i] = CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, halfwire);
        }
    }
}

//----------------------------- surface expansion pad Function -----------------------------//
void Teeth::compute_minimum_surface_pad(int index, double surface_pad)
{
    if (m_tooth[index]->get_init_flag())
    {
        if (index == 0)
        {
            Polyhedron& wire_outlines = m_wires[0]->get_wire_right_shape();
            Polyhedron& wire_shape = m_wires[0]->get_wire_extrude_shape();
            m_tooth[index]->compute_minimum_surface_pad(surface_pad);
        }
        else if (index == 5)
        {
            Polyhedron& wire_outlines = m_wires[4]->get_wire_left_shape();
            Polyhedron& wire_shape = m_wires[4]->get_wire_extrude_shape();
            m_tooth[index]->compute_minimum_surface_pad(surface_pad);
        }
        else
        {
            Polyhedron& left_wire_outlines = m_wires[index - 1]->get_wire_left_shape();
            Polyhedron& right_wire_outlines = m_wires[index]->get_wire_right_shape();
            Polyhedron& left_wire_shape = m_wires[index - 1]->get_wire_extrude_shape();
            Polyhedron& right_wire_shape = m_wires[index]->get_wire_extrude_shape();
            m_tooth[index]->compute_minimum_surface_pad(surface_pad);
        }

    }
}
//----------------------------- Inside / Convexity Help Function -----------------------------//

void Teeth::compute_smoothed_field(int smooth_range, VertDoubleMap& init_map, VertDoubleMap& value_map)
{
    // initialize value map
    value_map.clear();
    for (auto& elem : init_map)
        value_map.insert({ elem.first, elem.second });

    // Compute smoothed intersection map
    for (int i = 0; i < smooth_range; i++)
    {
        VertDoubleMap new_value_map;

        for (auto& elem : value_map)
        {
            Vertex_handle vd = elem.first;
            double count = 1., value = elem.second;
            Halfedge_vertex_circulator he = vd->vertex_begin();
            do {
                count = count + 1.;
                value = value + value_map[he->opposite()->vertex()];
            } while (++he != vd->vertex_begin());
            double new_value = value / count;
            new_value_map.insert({ vd, new_value });
        }

        value_map.swap(new_value_map);
        new_value_map.clear();
    }
}

double Teeth::locate_and_evaluate_function(Point& query, Facet_handle face, VertDoubleMap& vert_func_map)
{
    Halfedge_facet_circulator he = face->facet_begin();
    Vertex_handle va = he->vertex();
    Vertex_handle vb = he->next()->vertex();
    Vertex_handle vc = he->next()->next()->vertex();

    double a, b, c;
    barycentric_coordinates(query, va->point(), vb->point(), vc->point(), a, b, c);
    double value = a * vert_func_map[va] + b * vert_func_map[vb] + c * vert_func_map[vc];

    return value;
}

void Teeth::barycentric_coordinates(Point& p, Point& pa, Point& pb, Point& pc, double& a, double& b, double& c)
{
    double area_abc = std::sqrt(CGAL::squared_area(pa, pb, pc));
    a = std::sqrt(CGAL::squared_area(p, pb, pc)) / area_abc;
    b = std::sqrt(CGAL::squared_area(p, pa, pc)) / area_abc;
    c = std::sqrt(CGAL::squared_area(p, pa, pb)) / area_abc;
}

} // namespace OrthoLab
