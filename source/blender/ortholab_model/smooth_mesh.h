#pragma once

#include "types.h"
namespace OrthoLab {

void smooth_mesh(Polyhedron& poly_in_out, int iteration, Polyhedron& fix_boundary, const std::string& algo="LaplacianHC");

}// namespace OrthoLab