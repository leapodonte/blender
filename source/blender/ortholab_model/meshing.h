#ifndef MY_MESHING_H
#define MY_MESHING_H

#include "types.h"


namespace OrthoLab {

void extract_inside_vertices_from_function(Polyhedron& poly, PointList& poly_points);

bool find_level_set_point(Point& point_a, double value_a, Point point_b, 
     double value_b, double value, Point& p);

int extract_isovertices_from_function(Polyhedron& poly,
    VertDoubleMap& vert_value_map,
    double isovalue,
    PointList& poly_points,
    HEdgeIntMap& he_point_map);

// Extract a selected positive part from an implicit function
void extract_positive_faces_from_function(VertDoubleMap &vert_value_map,
                                          double value,
                                          VertIntMap &vert_ind_map,
                                          HEdgeIntMap &he_point_map,
                                          FacetList &candidate_facets,
                                          PolyList &poly_faces);

void copy_polyehdral_surface(Polyhedron& poly_in, Polyhedron& poly_out);

void merge_polyehdral_surface(Polyhedron& poly_in_1, Polyhedron& poly_in_2, Polyhedron& poly_out);

bool solve_selfintersection(Polyhedron& poly_in, Polyhedron& poly_out);

void alpha_wrap_mesh(Polyhedron &poly_in, Polyhedron &poly_out);

void expand_mesh(Polyhedron &poly_in, double distance);
// fill the holes of all polyhedron
// shell = 1 hole ==> 1
// solid = 0 hole ==> 0
void compute_enlarge_msh(Polyhedron &poly_in, Polyhedron &poly_out, double ratio);

void hole_filling(Polyhedron &m_tooth_poly, const bool shell_solid);

void copy_faces_to_triangles(Polyhedron& tooth_poly, TriangleList& triangles);

// --------------------- projection 2D/3D

Point_2 to_2d(Plane& plane, Point& center, Vector& base_1, Vector& base_2, Point& point_3d);

Point to_3d(Point& center, Vector& base_1, Vector& base_2, Point_2& point_2d);

// ---------------------shell form related services-------------------------
Vector compute_average_normal(const Polyhedron& polyhedron);
// for a shell like polyhedron mesh, 
// 1.find its linear least squares fitting plane
// 2.find the plane normal direction which is in the same direction of the shell normal
// 3.find the center point of the shell projection on the fitting plan
bool shell_extract_projection_situation(Polyhedron& poly_in, Plane& plane, Vector& normal, Point& plane_projected_center);


// --------------------- cylinder
void compute_cylinder_surface(PointList& polylines, double radius);

void construct_cone(Polyhedron &cone, Point base, Vector base_direction, Vector base_normal, 
                    Vector cone_direction, double radius, double height, int num_segments);

// --------------------- plane
void compute_plane_mesh(double scale, Plane& plane, Point& center, Point& right, Polyhedron& plane_mesh);


Transformation compute_alignment_transform(const Plane &first_plane,
                                           const Segment &first_segment, 
                                           const Plane &second_plane,
                                           const Segment &second_segment);

void transform_mesh(Polyhedron& poly_in_out, Transformation& trans);                  
} // namespace OrthoLab

#endif
