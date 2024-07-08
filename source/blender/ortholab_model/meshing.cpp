#include "meshing.h"

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>
#include <CGAL/Polygon_mesh_processing/intersection.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/alpha_wrap_3.h>
#include "types.h"
#include "base_helpers.h"


namespace OrthoLab {

void extract_inside_vertices_from_function(Polyhedron& poly, PointList& poly_points)
{
    for(Vertex_handle vd: vertices(poly))
        poly_points.push_back(vd->point());
}

bool find_level_set_point(Point& point_a, double value_a, Point point_b, double value_b, double value, Point& p)
{
    double f1 = value_a - value;
    double f2 = value_b - value;

    // v1 and v2 have opposite function values while at least one is not 0.
    if((f1 * f2 <= 0) && !((std::abs(f1) < 1e-8) && (std::abs(f2) < 1e-8)))
    {
        if(std::abs(f1) < 1e-8) 
            p = point_a;
        else if(std::abs(f2) < 1e-8) 
            p = point_b;
        else if(f1 > 0 && f2 < 0) {
            double ratio = (0. - f1) / (f2 - f1);
            p = point_a + ratio * (point_b - point_a);
        }
        else {
            double ratio = (0. - f2) / (f1 - f2);
            p = point_b + ratio * (point_a - point_b);
        }
        return true;
    }
    return false;
}

int extract_isovertices_from_function(  Polyhedron& poly, 
                                        VertDoubleMap& vert_value_map, 
                                        double isovalue, 
                                        PointList& poly_points, 
                                        HEdgeIntMap& he_point_map)
{
    int index = (int)poly_points.size();

    for(Halfedge_iterator hedge = poly.halfedges_begin(); hedge != poly.halfedges_end(); ++hedge) 
    {
        Vertex_handle vert_a = hedge->vertex();
        Vertex_handle vert_b = hedge->opposite()->vertex();
        double value_a = vert_value_map[vert_a];
        double value_b = vert_value_map[vert_b];

        if(value_a > isovalue && value_b < isovalue)
        {
            Point point;
            bool flag = find_level_set_point(   vert_a->point(), 
                                                value_a, 
                                                vert_b->point(), 
                                                value_b, 
                                                isovalue, 
                                                point);
            if(flag)
            {
                poly_points.push_back(point);
                he_point_map.insert({hedge, index});
                index++;
            }
        }
    }

    return index;
}
// process and extract positive faces from a given implicit function based on a specific value threshold.
// The function operates on 3D mesh data, distinguishing mesh facets (faces) based on their vertices' 
// values relative to the given threshold (value). This extraction process is typically used in 
// applications such as isosurface extraction, mesh segmentation, or similar areas. 
// Extract a selected positive part from an implicit function
void extract_positive_faces_from_function(VertDoubleMap &vert_value_map,// Maps vertices to their corresponding values in the implicit function.
                                          double value,                 // The threshold value used to categorize vertices.
                                          VertIntMap &vert_ind_map,     // Maps vertices to their indices.
                                          HEdgeIntMap &he_point_map,    // Maps half-edges to point indices, presumably used for intermediate vertices created during the mesh processing. 
                                          FacetList &candidate_facets,  // A list of facets (faces) considered for processing.
                                          PolyList &poly_faces)         // The output list of polygons (faces) that represent the positive parts of the mesh
{
// Processing Logic: The function iterates through each candidate facet, 
// categorizing its vertices based on whether their associated value is above or below the threshold. 
// It distinguishes between four cases based on the count of positive and negative vertices:

// All negative (---): Skips the facet.
// All positive (+++): Adds the facet directly to the output.
// One positive, two negatives (+--): Processes to create a new face that captures the positive portion.
// Two positives, one negative (++-): Splits into new faces to represent the positive parts.
    for(Facet_handle fd: candidate_facets)
    {
        // Count positive and negative number
        HEdgeList hedges;
        int pos_count = 0, neg_count = 0;
        Halfedge_facet_circulator he = fd->facet_begin();
        do {
            hedges.push_back(he);
            if(vert_value_map[he->vertex()] >= value)
                pos_count++;
            else
                neg_count++;
        } while (++he != fd->facet_begin());

        // case ---
        if(neg_count == 3)
            continue;

        // case +++
        if(pos_count == 3)   
        {
            IntList face;
            for(int i = 0; i < 3; i++)
                face.push_back(vert_ind_map[hedges[i]->vertex()]);
            poly_faces.push_back(face);
        }
        
        // case +--
        if(pos_count == 1 && neg_count == 2)   
        {
            IntList face;
            for(int i = 0; i < 3; i++)
            {
                Halfedge_handle hedge = hedges[i];
                double va = vert_value_map[hedge->vertex()];
                double vb = vert_value_map[hedge->opposite()->vertex()];
                if(va >= value && vb < value)
                    face.push_back(vert_ind_map[hedge->vertex()]);
                else if(va < value && vb >= value)
                    face.push_back(he_point_map[hedge->opposite()]);
                else
                    face.push_back(he_point_map[hedge->next()]);
            }
            poly_faces.push_back(face);
        }
        
        // case ++-
        if(pos_count == 2 && neg_count == 1)   
        {
            IntList face_1(3, 0), face_2(3, 0);
            for(int i = 0; i < 3; i++)
            {
                Halfedge_handle hedge = hedges[i];
                double va = vert_value_map[hedge->vertex()];
                double vb = vert_value_map[hedge->opposite()->vertex()];
                if(va >= value && vb >= value)
                {
                    int index = vert_ind_map[hedge->vertex()];
                    face_1[2] = index;
                    face_2[2] = index;
                }
                else if(va >= value && vb < value)
                {
                    int index = he_point_map[hedge];
                    face_1[1] = index;
                    face_2[0] = index;
                    face_2[1] = vert_ind_map[hedge->vertex()];
                }
                else
                    face_1[0] = he_point_map[hedge->opposite()];
            }
            poly_faces.push_back(face_1);
            poly_faces.push_back(face_2);
        }
    }
}

void copy_polyehdral_surface(Polyhedron& poly_in, Polyhedron& poly_out)
{
    poly_out.clear();
    PointList poly_points;
    PolyList  poly_faces;
    VertIntMap vert_index_map;
    int count = 0;

    for(Vertex_handle vd: vertices(poly_in))
    {
        Point p = vd->point();
        poly_points.push_back(p);
        vert_index_map.insert({vd, count});
        count++;
    }

    for(Facet_iterator face = poly_in.facets_begin(); face != poly_in.facets_end(); ++face) 
    {
        Halfedge_facet_circulator he = face->facet_begin();
        IntList face_inds;
        do {
            Vertex_handle vd = he->vertex();
            face_inds.push_back(vert_index_map[vd]);
        } while (++he != face->facet_begin());
        poly_faces.push_back(face_inds);
    }

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, poly_out);
}

void merge_polyehdral_surface(Polyhedron& poly_in_1, Polyhedron& poly_in_2, Polyhedron& poly_out)
{
    PointList poly_points;
    PolyList  poly_faces;
    VertIntMap vert_index_map;
    int count = 0;

    for(Vertex_handle vd: vertices(poly_in_1))
    {
        Point p = vd->point();
        poly_points.push_back(p);
        vert_index_map.insert({vd, count});
        count++;
    }

    for(Vertex_handle vd: vertices(poly_in_2))
    {
        Point p = vd->point();
        poly_points.push_back(p);
        vert_index_map.insert({vd, count});
        count++;
    }

    for(Facet_iterator face = poly_in_1.facets_begin(); face != poly_in_1.facets_end(); ++face) 
    {
        Halfedge_facet_circulator he = face->facet_begin();
        IntList face_inds;
        do {
            Vertex_handle vd = he->vertex();
            face_inds.push_back(vert_index_map[vd]);
        } while (++he != face->facet_begin());
        poly_faces.push_back(face_inds);
    }

    for(Facet_iterator face = poly_in_2.facets_begin(); face != poly_in_2.facets_end(); ++face) 
    {
        Halfedge_facet_circulator he = face->facet_begin();
        IntList face_inds;
        do {
            Vertex_handle vd = he->vertex();
            face_inds.push_back(vert_index_map[vd]);
        } while (++he != face->facet_begin());
        poly_faces.push_back(face_inds);
    }

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(poly_points, poly_faces, poly_out);
}

bool solve_selfintersection(Polyhedron &poly_in, Polyhedron &poly_out)
{
    if (CGAL::Polygon_mesh_processing::does_self_intersect(poly_in))
    {
        std::cout << "      Find self-intersection!" << std::endl;
        alpha_wrap_mesh(poly_in, poly_out);
        std::cout << "      Alpha wrapper done!" << std::endl;
        return true;
    }
    copy_polyehdral_surface(poly_in, poly_out);
    return false;
}

// bool solve_is_closed(Polyhedron &poly_in, Polyhedron &poly_out)
// {
//     if (CGAL::Polygon_mesh_processing::is_closed(poly_in))
//     {
//         std::cout << "      Find self-intersection!" << std::endl;
//         alpha_wrap_mesh(poly_in, poly_out);
//         std::cout << "      Alpha wrapper done!" << std::endl;
//         return true;
//     }
//     return false;
// }


// manipulate a 3D mesh by expanding it uniformly in all directions based on the normals of its surfaces. 
// This can be useful in various applications, including graphics, 3D printing, and computational geometry,
// to adjust the size of an object or create a slightly enlarged version of it. 
// void expand_mesh(Polyhedron &poly_in, double distance)
// {
//     // Step 1: Compute the normals for each face and vertex displacements
//     // Compute Normals and Vertex Displacements:

// // Iterates over all faces (Facet_iterator) of the input mesh poly_in.
// // Computes the normal of each face using the CGAL::Polygon_mesh_processing::compute_face_normal function.
// // Normalizes this normal vector.
// // Updates the displacement for each vertex of the face by adding the normalized face normal to 
// // a map (vertexDisplacements) that tracks the total displacement for each vertex.
//     std::map<Vertex_handle, Vector> vertexDisplacements;
//     for (Facet_iterator f = poly_in.facets_begin(); f != poly_in.facets_end(); ++f)
//     {
//         Vector normal = CGAL::Polygon_mesh_processing::compute_face_normal(f, poly_in);
//         normal = normal / std::sqrt(normal * normal); // Normalize
//         // Iterate through the vertices of the face and update the displacement
//         Polyhedron::Halfedge_around_facet_circulator h = f->facet_begin();
//         do
//         {
//             vertexDisplacements[h->vertex()] += normal;
//         } while (++h != f->facet_begin());
//     }

//     // Normalize the displacement vectors and displace vertices
//     for (auto &vd : vertexDisplacements)
//     {
//         Vector displacement = vd.second / std::sqrt(vd.second.squared_length());
//         vd.first->point() = vd.first->point() + displacement * distance; // Displace by x mm
//     }

//     //alpha wrapping
//     //alpha_wrap_mesh(poly_in, poly_in);
// }

void compute_enlarge_msh(Polyhedron &poly_in, Polyhedron &poly_out, double dist)
{
    copy_polyehdral_surface(poly_in, poly_out);

    VertVecMap vert_normals;

    for(Vertex_handle vd: vertices(poly_out))
    {
        Vector normal = CGAL::Polygon_mesh_processing::compute_vertex_normal(vd, poly_out);
        vert_normals.insert({vd, normal});
    }
    
    for(Vertex_handle vd: vertices(poly_out))
    {
        Point p = vd->point() + vert_normals[vd] * dist;
        vd->point() = p;
    }
}


void alpha_wrap_mesh(Polyhedron &poly_in, Polyhedron &poly_out)
{
    Bbox bbox = CGAL::Polygon_mesh_processing::bbox(poly_in);
    double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                   CGAL::square(bbox.ymax() - bbox.ymin()) +
                                   CGAL::square(bbox.zmax() - bbox.zmin()));

    double alpha = diag_length / 200;
    double offset = diag_length / 2000;

    CGAL::alpha_wrap_3(poly_in, alpha, offset, poly_out);
}
double _compute_length_for_hole(Halfedge_handle he, Polyhedron& mesh)
{
    double length = 0.;

    for (Halfedge_handle hc : CGAL::halfedges_around_face(he, mesh))
        length += std::sqrt(CGAL::squared_distance(hc->vertex()->point(), hc->opposite()->vertex()->point()));

    return length;
}
void hole_filling(Polyhedron &m_tooth_poly, const bool shell_solid)
{
 // hole filling
    HEdgeList border_cycles;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m_tooth_poly, std::back_inserter(border_cycles));

    if (border_cycles.size() > (int)shell_solid) // contain holes
    {
        DoubleList border_lengths;
        for (int i = 0; i < border_cycles.size(); i++)
        {
            double length = _compute_length_for_hole(border_cycles[i], m_tooth_poly);
            border_lengths.push_back(length);
        }

        unsigned int nb_holes = 0;
        std::vector<int> indices(border_lengths.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](int A, int B) -> bool { return border_lengths[A] > border_lengths[B]; });

        for (int i = (int)shell_solid ; i < indices.size(); i++)
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
}

void copy_faces_to_triangles(Polyhedron& tooth_poly, TriangleList& triangles)
{
    for (Facet_iterator face = tooth_poly.facets_begin(); face != tooth_poly.facets_end(); ++face)
    {
        Halfedge_facet_circulator he = face->facet_begin();
        PointList triangle;
        do {
            triangle.push_back(he->vertex()->point());
        } while (++he != face->facet_begin());

        if (triangle.size() == 3)
        {
            Triangle tri(triangle[0], triangle[1], triangle[2]);
            triangles.push_back(tri);
        }
    }
}

Point_2 to_2d(Plane& plane, Point& center, Vector& base_1, Vector& base_2, Point& point_3d)
{
    Vector offset = plane.projection(point_3d) - center;
    double x_2d = offset * base_1;
    double y_2d = offset * base_2;
    return Point_2(x_2d, y_2d);
}

Point to_3d(Point& center, Vector& base_1, Vector& base_2, Point_2& point_2d)
{
    Point p_3d = center + point_2d.x() * base_1 + point_2d.y() * base_2;
    return p_3d;
}

Vector compute_average_normal(const Polyhedron& polyhedron) 
{
    Vector sum(0, 0, 0); // Initialize sum of normals
    for (Polyhedron::Facet_const_iterator fit = polyhedron.facets_begin(); fit != polyhedron.facets_end(); ++fit) {
        Polyhedron::Halfedge_around_facet_const_circulator hfc = fit->facet_begin();
        // Assuming triangular facets for simplicity
        Point a = hfc->vertex()->point();
        Point b = (++hfc)->vertex()->point();
        Point c = (++hfc)->vertex()->point();
        
        // Compute the normal of the triangle
        Vector normal = CGAL::cross_product(b - a, c - a);
        normal = normal / std::sqrt(normal * normal); // Normalize
        sum = sum + normal; // Add to the sum
    }
    
    if (sum.squared_length() > 0) {
        sum = sum / std::sqrt(sum * sum); // Normalize the sum
    }
    return sum;
}

// for a shell like polyhedron mesh, 
// 1.find its linear least squares fitting plane
// 2.find the plane normal direction which is in the same direction of the shell normal
// 3.find the center point of the shell projection on the fitting plan
bool shell_extract_projection_situation(Polyhedron& poly_in, Plane& plane,
 Vector& normal, Point& plane_projected_center)
{
    if (poly_in.empty())
    {
        std::cout << "    Error: The input polyhedron is empty!" << std::endl;
        return false;
    }

    TriangleList candidate_faces;
    copy_faces_to_triangles(poly_in, candidate_faces);
     
    CGAL::linear_least_squares_fitting_3(candidate_faces.begin(), 
            candidate_faces.end(), plane, CGAL::Dimension_tag<2>());

    normal = plane.orthogonal_vector();
    normal = normal / std::sqrt(normal.squared_length());

    // decide the direction of the normal
    Vector avg = compute_average_normal(poly_in);
    if (avg * normal < 0)
    {
        normal = -normal;
    }

    // find the tooth center and the projection to the fitting plan
    Point tooth_center = CGAL::centroid(candidate_faces.begin(), candidate_faces.end());
    plane_projected_center = plane.projection(tooth_center);

    return true;
}


void construct_cone(Polyhedron &cone, Point base, Vector base_direction, Vector base_normal, 
                    Vector cone_direction, double radius, double height, int num_segments)
{ 
    std::vector<Point> base_vertices;

    // Calculate base vertices
    Vector base1 = base_direction;
    Vector base2 = CGAL::cross_product(base_normal, base1);
    base1 = CGAL::cross_product(base2, cone_direction);
    base1 = base1 / std::sqrt(base1.squared_length()); // Normalize
    base2 = base2 / std::sqrt(base2.squared_length()); // Normalize
    for (int i = 0; i < num_segments; ++i)
    {
        double angle = 2 * M_PI * i / num_segments;
        Vector offset = radius * (std::cos(angle) * base1 + std::sin(angle) * base2);
        base_vertices.push_back(base + offset);
    }

    // apex point of the cone
    cone_direction = cone_direction / std::sqrt(cone_direction.squared_length()); // Normalize
    Point tip = base + cone_direction * height;
    // Begin constructing the cone
    Polyhedron::Halfedge_handle h = cone.make_triangle(base_vertices[0], base_vertices[1], tip);

    // Create the base and side faces
    for (int i = 0; i < num_segments; ++i) {
        // Add a triangle for the base
        cone.make_triangle(base_vertices[(i + 1) % num_segments],
                           base_vertices[i],
                           base);

        // Add a triangle for the side
        cone.make_triangle(base_vertices[i],
                           base_vertices[(i + 1) % num_segments],
                           tip);
    }
}

void compute_plane_mesh(double scale, Plane& plane, Point& center, Point& right, Polyhedron& plane_mesh)
{
    Point plane_center = plane.projection(center);
    Point plane_right = plane.projection(right);
    Vector horizontal = plane_right - plane_center;
    horizontal = horizontal / std::sqrt(horizontal.squared_length());
    Vector normal = plane.orthogonal_vector();
    Vector vertical = CGAL::cross_product(horizontal, normal);
    vertical = vertical / std::sqrt(vertical.squared_length());

    PointList corners;
    Point left_down = plane_center - scale * horizontal - scale * vertical;
    Point left_up = plane_center - scale * horizontal + scale * vertical;
    Point right_up = plane_center + scale * horizontal + scale * vertical;
    Point right_down = plane_center + scale * horizontal - scale * vertical;
    corners.push_back(left_down);
    corners.push_back(left_up);
    corners.push_back(right_up);
    corners.push_back(right_down);

    PolyList faces;
    faces.push_back(IntList({ 0, 1, 3 }));
    faces.push_back(IntList({ 1, 2, 3 }));

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(corners, faces, plane_mesh);
}

EMatrix3d matrix_from_two_vectors(const Vector& v1, const Vector& v2)
{
    EVector3d eigen_v1(v1.x(), v1.y(), v1.z());
    EVector3d eigen_v2(v2.x(), v2.y(), v2.z());
    EMatrix3d rotation_matrix =   
        EQuaterniond().setFromTwoVectors(eigen_v1,eigen_v2).toRotationMatrix();
    return rotation_matrix;
}
Vector matrix_transform(const EMatrix3d matrix, const Vector& v)
{
    EVector3d eigen_v(v.x(), v.y(), v.z());
    EVector3d eigen_v_transformed = matrix * eigen_v; 
    return Vector(eigen_v_transformed.x(), eigen_v_transformed.y(), eigen_v_transformed.z());
}
Transformation compute_alignment_transform(const Plane &first_plane,
                                           const Segment &first_segment, 
                                           const Plane &second_plane,
                                           const Segment &second_segment)
{
    auto first_p = first_plane.projection(first_segment.source());
    auto second_p = second_plane.projection(second_segment.source());

    // std::cout << "compute_alignment_transform" << std::endl;
    // Vector v1 = first_p - first_segment.source();
    // // std::cout << "" << v1.squared_length() << std::endl;
    // Vector v2 = second_p - second_segment.source();
    // std::cout << "" << v2.squared_length() << std::endl;

    // Step 1: Compute Normal Vectors
    Vector normal1 = first_plane.orthogonal_vector();
    Vector normal2 = second_plane.orthogonal_vector();
 
    // Compute rotation matrix using Eigen
    EMatrix3d rotation_matrix1 = matrix_from_two_vectors(normal2, normal1);

    //test
    Vector n2_t = matrix_transform(rotation_matrix1, normal2);
    // std::cout << "test normal 1:" << normal1 << std::endl;
    // std::cout << "test normal 2:" << normal2 << std::endl;
    // std::cout << "test transformed normal 2:" << n2_t << std::endl;
    
    
    // step2 rotate the second segment 
    Vector seg2_dir = second_segment.target() - second_segment.source();
    seg2_dir = seg2_dir / std::sqrt(seg2_dir.squared_length());
    Vector seg1_dir = first_segment.target() - first_segment.source();
    seg1_dir = seg1_dir / std::sqrt(seg1_dir.squared_length());
    Vector seg2_dir_rotated = matrix_transform(rotation_matrix1,seg2_dir); 
    EMatrix3d rotation_matrix2 = matrix_from_two_vectors(seg2_dir_rotated, seg1_dir);
    
    EMatrix3d fm = rotation_matrix2*rotation_matrix1;
    // Vector seg2_t = matrix_transform(final_rotation_matrix, seg2_dir);
    // std::cout << "test seg 1 dir:" << seg1_dir << std::endl;
    // std::cout << "test seg2_dir:" << seg2_dir << std::endl;
    // std::cout << "test transformed seg 2:" << seg2_t << std::endl;
    // seg2_t = matrix_transform(rotation_matrix2, matrix_transform(rotation_matrix1, seg2_dir));
    // std::cout << "test transformed seg 2:" << seg2_t << std::endl;
    // translate the segment 
    auto seg2_s = second_segment.source(); 
    auto seg1_t =  first_segment.target();
    Vector tvec = Vector(seg1_t.x(), seg1_t.y(), seg1_t.z())- 
                  matrix_transform(fm, Vector(seg2_s.x(), seg2_s.y(), seg2_s.z()));

    Transformation translation_transform(CGAL::TRANSLATION, tvec);
    // Create a Transformation object
    Transformation rotate_transform(fm(0,0), fm(0,1), fm(0,2), 0,
                                    fm(1,0), fm(1,1), fm(1,2), 0,
                                    fm(2,0), fm(2,1), fm(2,2), 0, 1);

    // Vector seg2_t =  rotate_transform (seg2_dir);
    // std::cout << "test seg 1 dir:" << seg1_dir << std::endl;
    // std::cout << "test seg2_dir:" << seg2_dir << std::endl;
    // std::cout << "test transformed seg 2:" << seg2_t << std::endl; 
    // translate the segment 


    Transformation final_transform =  translation_transform*rotate_transform;

    Point p1 = first_segment.target();
    Point p2 = second_segment.source();
    Point p2_t = final_transform.transform(p2); 

    // std::cout << "test seg 1 target:" << p1 << std::endl;
    // std::cout << "test seg 2 source:" << p2 << std::endl;
    // std::cout << "test seg 2 source t:" << p2_t << std::endl;
 
    Point p0 = first_segment.source();
    Point p3 = second_segment.target();
    Point p3_t = final_transform.transform(p3); 

    // std::cout << "test seg 1 dir:" << seg1_dir << std::endl;
    // std::cout << "test seg2_dir:" << seg2_dir << std::endl;
    Vector seg2_dir_t = p3_t - p2_t;
    seg2_dir_t = seg2_dir_t / std::sqrt(seg2_dir_t.squared_length());
    // std::cout << "test seg2_dir_t:" << seg2_dir_t << std::endl;

    return final_transform;
}

void transform_mesh(Polyhedron &poly_in_out, Transformation &trans)
{
    for (Polyhedron::Vertex_iterator it = poly_in_out.vertices_begin();
         it != poly_in_out.vertices_end(); ++it)
    {
        it->point() = trans.transform(it->point());
    }
}

} // namespace OrthoLab
