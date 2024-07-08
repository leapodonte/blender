#include "extrusion.h"

#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>


#include "types.h"
#include "meshing.h"

namespace OrthoLab {

Bottom::Bottom() {}
void Bottom::operator() (Vertex_handle vin, Vertex_handle vout) const {}

Top::Top(VertVec_property_map nmap, double vlen)
    : nmap(nmap), vlen(vlen)
{}
void Top::operator()(Vertex_handle vin, Vertex_handle vout) const
{
    Vector normal = nmap[vin];
    vout->point() = vin->point() + vlen * normal;
}

void compute_extrusion(Polyhedron& poly_in, Polyhedron& poly_out, double height)
{
    std::cout << "      Begin!" << std::endl;

    VertVecMap vnormals;
    CGAL::Polygon_mesh_processing::compute_vertex_normals(poly_in, boost::make_assoc_property_map(vnormals));
    Bottom bottom;
    Top top(vnormals, height);
    CGAL::Polygon_mesh_processing::extrude_mesh(poly_in, poly_out, bottom, top);
    CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(poly_out);

    std::cout << "      Extrusion done!" << std::endl;

    Polyhedron poly_out_wrap;
    if (solve_selfintersection(poly_out, poly_out_wrap))
    {
        copy_polyehdral_surface(poly_out_wrap, poly_out);
    }
}

void compute_extrusion(Polyhedron& poly_in, Polyhedron& poly_out, Vector& normal, double height)
{
    std::cout << "      Begin!" << std::endl;

    VertVecMap vnormals;
    for(Vertex_handle vd: vertices(poly_in))
        vnormals.insert({vd, normal});
    
    Bottom bottom;
    Top top(vnormals, height);
    CGAL::Polygon_mesh_processing::extrude_mesh(poly_in, poly_out, bottom, top);
    CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(poly_out);

    std::cout << "      Extrusion done!" << std::endl;

    Polyhedron poly_out_wrap;
    if (solve_selfintersection(poly_out, poly_out_wrap))
    {
        copy_polyehdral_surface(poly_out_wrap, poly_out);
    }
}

void compute_extrusion(Polyhedron& poly_in, Polyhedron& poly_out, Vector& normal, VertDoubleMap& heights)
{
    std::cout << "      Begin!" << std::endl;

    VertVecMap vnormals;
    for(Vertex_handle vd: vertices(poly_in))
    {
        double height = heights[vd];
        Vector vec = normal * height;
        vnormals.insert({vd, vec});
    }
    
    Bottom bottom;
    Top top(vnormals, 1.);
    CGAL::Polygon_mesh_processing::extrude_mesh(poly_in, poly_out, bottom, top);
    CGAL::Polygon_mesh_processing::orient_to_bound_a_volume(poly_out);

    std::cout << "      Extrusion done!" << std::endl;

    Polyhedron poly_out_wrap;
    if (solve_selfintersection(poly_out, poly_out_wrap))
    {
        copy_polyehdral_surface(poly_out_wrap, poly_out);
    }
}

void compute_swept_extrusion(Polyhedron& poly_in, Polyhedron& poly_out, Vector& normal, double height)
{
    std::cout << "      Begin!" << std::endl;
    Vector direction = normal * height;
    CGAL::Polygon_mesh_processing::extrude_mesh(poly_in, poly_out, direction);
    std::cout << "      Extrusion done!" << std::endl;

    Polyhedron poly_out_wrap;
    if (solve_selfintersection(poly_out, poly_out_wrap))
    {
        copy_polyehdral_surface(poly_out_wrap, poly_out);
    }
}

} // namespace OrthoLab
