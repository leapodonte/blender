#pragma once

#include "types.h"

#ifdef max
#undef max
#endif

#ifdef min
#undef min
#endif

#include "cmesh.h"


namespace OrthoLab
{

    bool mesh_import(const std::string &stl_filename, Polyhedron &poly_in_out);

    void mesh_export(const Polyhedron &poly_in, const std::string &stl_filename);

    void mesh_export(const CMeshO &poly_in, const std::string &stl_filename);

    // void exportPolyMesh(Polyhedron&  polyhedron, QJsonObject &meshObject);

    // create CMesh from Poly mesh
    void cgal2vcg(const Polyhedron &poly_in, CMeshO &poly_out);

    // Synch Poly Mesh with CMesh
    void update_polymesh_from_cmesh(Polyhedron &poly_in_out,
                                    const CMeshO &poly_in);

    void savePointsToXYZ(const PointList &points, const std::string &filename);

    void readPointsFromXYZ(const std::string &filename, PointList &points);

} // namespace OrthoLab