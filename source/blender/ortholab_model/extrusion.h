#pragma once

#include "types.h"

namespace OrthoLab {

struct Bottom
{
    Bottom();
    void operator() (Vertex_handle vin, Vertex_handle vout) const;
};

struct Top
{
    Top(VertVec_property_map nmap, double vlen);
    void operator()(Vertex_handle vin, Vertex_handle vout) const;
    VertVec_property_map nmap;
    double vlen;
};

void compute_extrusion(Polyhedron& poly_in, Polyhedron& poly_out, double height);

void compute_extrusion(Polyhedron& poly_in, Polyhedron& poly_out, Vector& normal, double height);

void compute_extrusion(Polyhedron& poly_in, Polyhedron& poly_out, Vector& normal, VertDoubleMap& heights);

void compute_swept_extrusion(Polyhedron& poly_in, Polyhedron& poly_out, Vector& normal, double height);

} // namespace OrthoLab