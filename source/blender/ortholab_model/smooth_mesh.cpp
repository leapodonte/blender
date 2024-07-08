#include "smooth_mesh.h"

#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>

#include "types.h"
#include "meshing.h"
#include "mesh_io.h"
#include "cmesh.h"
#include "vcg/complex/algorithms/smooth.h"

namespace OrthoLab {

void smooth_mesh(Polyhedron& poly_in_out, int iteration, Polyhedron& fix_boundary, const std::string& algo)
{
    std::cout << "smooth_mesh PAD_SMOOTH_ITERATIONS: " << iteration << std::endl;
    
    // remesh first
    Polyhedron poly_remeshed;
    alpha_wrap_mesh(poly_in_out, poly_remeshed);
    // transfert to vcg mesh
    CMeshO temp_mesh;
    cgal2vcg(poly_remeshed, temp_mesh);
    // apply the fix boundary
    // create AABB tree for the boundary
    AABBTree boundary_aabb_tree;
    boundary_aabb_tree.insert(faces(fix_boundary).first, faces(fix_boundary).second, fix_boundary);
    boundary_aabb_tree.accelerate_distance_queries();
    auto vi = poly_remeshed.vertices_begin();
    for (auto v = temp_mesh.vert.begin(); v != temp_mesh.vert.end(); ++v)
    {
        double distance = boundary_aabb_tree.squared_distance(vi->point());
        if (distance > 0.0001) // smooth only points outside of boundary
        {   
            v->SetS();
        }
        vi++;
    }

    // smooth
    if(algo == "LaplacianHC")
        vcg::tri::Smooth<CMeshO>::VertexCoordLaplacianHC(temp_mesh, iteration, true);
    else if(algo == "Taubin")
        vcg::tri::Smooth<CMeshO>::VertexCoordTaubin(temp_mesh, iteration,0.9, -0.53, true);
    // transfert back
    copy_polyehdral_surface(poly_remeshed, poly_in_out);

    update_polymesh_from_cmesh(poly_in_out, temp_mesh);
    // mesh_export(temp_mesh, "smoothed_pad.stl");
}

} // namespace OrthoLab
