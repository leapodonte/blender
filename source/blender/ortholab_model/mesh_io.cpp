#include "mesh_io.h"
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/bbox.h>
#include <CGAL/IO/STL.h>

#define MESHLAB_SCALAR double
#include <vcg/complex/complex.h>
#include "io_mask.h"
#include "export_stl.h"
#include "types.h"
#include "meshing.h"

#include "spdlog/spdlog.h"
//#include "../include/odt.h"

unsigned long long hashFile(const std::string& input_path)
{
    unsigned long long result = 0;
    std::ifstream infile;
    infile.open(input_path);

    std::string line;
    for (; getline(infile, line, '\n'); )
    {
        size_t len = line.length();
        const char* it = line.c_str();
        const char* last = it + len;
        unsigned long long val = 0;
        for (; it != last; ++it)
            val = (val << 8 | val >> 24) + *it;
        result = val;
    }
    return result;
}

namespace OrthoLab
{

    bool mesh_import(const std::string &stl_filename,Polyhedron &poly_in_out)
    {
        if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(stl_filename,
                                                                  poly_in_out))
        {
            //spdlog::error("mesh_import Invalid mesh:{}", stl_filename);
            std::cout << "  Invalid mesh: " << stl_filename
                      << "." << std::endl;
                      return false;
        }
        return true;
    }



    // read mesh CGAL::Polygon_mesh_processing::IO::read_polygon_mesh
    void mesh_export(const Polyhedron &polyhedron, const std::string &stl_filename)
    {
        std::cout << "mesh_export file: " << stl_filename << std::endl;
        
        std::ofstream out(stl_filename);

        VertIntMap vert_index_map;
        PointList pointArray;

        for (Vertex_handle vd : vertices(polyhedron))
            pointArray.push_back(vd->point());

        std::vector<std::vector<int>> facets;

        int ivertices = 0;
        for (Vertex_handle v : vertices(polyhedron))
        {
            vert_index_map.insert({v, ivertices++});
        }

        for (Facet_handle face : faces(polyhedron))
        {
            std::vector<int> face_index;
            Halfedge_facet_circulator he = face->facet_begin();
            do
            {
                face_index.push_back(vert_index_map[he->vertex()]);
            } while (++he != face->facet_begin());
            facets.push_back(face_index);
        }

        std::cout << "mesh_export points: " << pointArray.size() << std::endl;
        std::cout << "mesh_export facets: " << facets.size() << std::endl;
        
        CGAL::IO::write_STL(out, pointArray, facets);
        out.close();

        spdlog::trace("mesh_export file hash: {0}", hashFile(stl_filename));
    }

    void mesh_export(const CMeshO &poly_in, const std::string &stl_filename)
    {
         std::cout << "mesh_export: " << stl_filename;
        int result = vcg::tri::io::ExporterSTL<CMeshO>::Save(poly_in, 
        stl_filename.c_str(), false, 0, NULL);

        spdlog::trace("mesh_export file hash: {0}", hashFile(stl_filename));
    }

    // void exportPolyMesh(Polyhedron&  polyhedron, QJsonObject &meshObject)
    // {
    //     VertIntMap vert_index_map;
    //     QJsonArray verticesArray;
    //     int ivertices = 0;
    //     for(Vertex_handle v: vertices(polyhedron))
    //     {
    //         //
    // vertexArray;
    //         verticesArray.append(QJsonValue(v->point().x()));
    //         verticesArray.append(QJsonValue(v->point().y()));
    //         verticesArray.append(QJsonValue(v->point().z()));
    //         //verticesArray.append(vertexArray);
    //         vert_index_map.insert({v, ivertices++});
    //     }
    //     std::cout << "exportPolyMesh: " << polyhedron.size_of_vertices()<< " vertext exported "<< std::endl;
        
    //     QJsonArray facesArray;
    //     for(Facet_handle face: faces(polyhedron)){
    //         //QJsonArray faceArray;
    //         Halfedge_facet_circulator he = face->facet_begin();
    //         do { 
    //             facesArray.append(QJsonValue(vert_index_map[he->vertex()]));
    //         } while (++he != face->facet_begin()); 
    //         //facesArray.append(faceArray);
    //     }
    //     std::cout << "exportPolyMesh: " << polyhedron.size_of_facets()<< " faces exported "<< std::endl;

    //     // Add the vertices and faces arrays to the meshObject
    //     meshObject["vertices"] = verticesArray;
    //     meshObject["faces"] = facesArray;
    // }

    void cgal2vcg(const Polyhedron& cgalPolyhedron, CMeshO& vcgMesh)
    {
        // Clear the VCG mesh to start fresh.
        vcgMesh.Clear();

        // Reserve space for vertices and faces to improve performance.
        vcgMesh.vert.reserve(cgalPolyhedron.size_of_vertices());
        vcgMesh.face.reserve(cgalPolyhedron.size_of_facets());
        CMeshO::VertexIterator vi=vcg::tri::Allocator<CMeshO>::AddVertices(vcgMesh,cgalPolyhedron.size_of_vertices());
        // Map to keep track of the corresponding vertices between CGAL and VCG.
        std::map<Polyhedron::Vertex_const_handle, CVertexO *> vMap;

        // Copy vertices.
        for (auto v = cgalPolyhedron.vertices_begin(); v != cgalPolyhedron.vertices_end(); ++v)
        {
            (*vi).P().Import(Point3m(v->point().x(), v->point().y(), v->point().z()));
            vMap[v] = &(*vi);
            ++vi;
        }

        // Copy faces.
        for (auto f = cgalPolyhedron.facets_begin(); f != cgalPolyhedron.facets_end(); ++f)
        {
            std::vector<CVertexO *> faceVertices;
            faceVertices.reserve(3); // Assuming triangular faces.

            // Iterate around the face to get the vertices.
            auto h = f->facet_begin();
            do
            {
                faceVertices.push_back(vMap[h->vertex()]);
            } while (++h != f->facet_begin());

            // Add the face to the VCG mesh.
            vcg::tri::Allocator<CMeshO>::AddFace(vcgMesh, faceVertices[0], faceVertices[1], faceVertices[2]);
        }

        // Update the normals and bounding box of the VCG mesh.
        // vcg::tri::UpdateNormals<vcgMesh>::PerVertexNormalizedPerFace(vcgMesh);
        // vcg::tri::UpdateNormals<vcgMesh>::PerFaceNormalized(vcgMesh);
        // vcg::tri::UpdateBounding<vcgMesh>::Box(vcgMesh);
    }

    // sync vectect position from CMesh to Poly mesh (after smoothing)
    void update_polymesh_from_cmesh(Polyhedron &poly_in_out,
                                    const CMeshO &poly_in)
    {

        if (poly_in_out.size_of_vertices() == poly_in.vert.size())
        {
            auto vit = poly_in.vert.begin();
            auto cit = poly_in_out.vertices_begin();

            for (; vit != poly_in.vert.end() && cit != poly_in_out.vertices_end(); ++vit, ++cit)
            {
                // Copy the position from VCG vertex to CGAL vertex
                cit->point() = Point(vit->P()[0], vit->P()[1], vit->P()[2]);
            }
        }
        else
        {
            std::cerr << "Meshes do not have the same topology." << std::endl;
        }
        return;
    }


void savePointsToXYZ(const PointList &points, const std::string &filename)
{
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << " for writing." << std::endl;
        return;
    }

    for (const Point& point : points) {
        file << point.x() << " " << point.y() << " " << point.z() << "\n";
    }

    file.close();

    spdlog::trace("savePointsToXYZ file hash: {0}", hashFile(filename));
}

void readPointsFromXYZ(const std::string &filename, PointList &points)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << " for reading." << std::endl;
        return;
    }

    points.clear();

    std::string line;
    double x, y, z;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        if (!(iss >> x >> y >> z)) {
            // Handle error or ignore malformed line
            continue;
        }
        points.push_back(Point(x, y, z));
    }

    file.close();
    spdlog::trace("readPointsFromXYZ file hash: {0}", hashFile(filename));
}

} // namespace OrthoLab
