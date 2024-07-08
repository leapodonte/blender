#include "pad.h"
#include "meshing.h"
#include "tooth.h"
#include "mesh_io.h"
#include "file_services.h"

namespace OrthoLab
{
    Pad::Pad(const std::string &id, const Tooth &tooth) : 
        m_tooth(tooth),
        SceneObject(id, tooth )
    {

    }

    Pad::~Pad()
    {
        reset();
    }

    void Pad::reset()
    {
        SceneObject::reset();
        m_polyline.clear();
        m_left_shape.clear();
        m_right_shape.clear(); 
        m_pad_extrude_shape.clear();
        m_pad_min_shape.clear(); 
    }

    void Pad::update()
    {
        SceneObject::update();

        
    }

    void Pad::save_to_db()
    {
    }

    void Pad::load_from_db()
    {
        std::string pad_path = m_tooth.get_base_path();
        pad_path = ensure_directory_exists(pad_path);
    }

    void Pad::set_cut_plane(Plane &cut_plane, Point &center, Point &right)
    {
        m_cut_plane.m_plane = cut_plane;
        m_cut_plane.m_center = center;
        m_cut_plane.m_right = right;
    }

    void Pad::set_polylines(PointList &intersection)
    {
        m_polyline.clear();
        for (int i = 0; i < intersection.size() - 1; i++)
            m_polyline.push_back(Segment(intersection[i], intersection[i + 1]));
    }

    void Pad::set_left_shape(Polyhedron &left_shape)
    {
        copy_polyehdral_surface(left_shape, m_left_shape);
    }
 
    void Pad::set_right_shape(Polyhedron &right_shape)
    {
        copy_polyehdral_surface(right_shape, m_right_shape);
    }

    // void Pad::set_extrude_shape(Polyhedron &extrude_shape)
    //{
    //     copy_polyehdral_surface(extrude_shape, m_pad_extrude_shape);
    // }
} // namespace OrthoLab
