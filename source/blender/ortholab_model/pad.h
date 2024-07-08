#ifndef MY_PAD_H
#define MY_PAD_H

#include "types.h"
#include "scene_object.h"

namespace OrthoLab
{
    class Tooth;
    /**
     * use to represent a generated pad**/
    class Pad : public SceneObject
    {
    private:
        bool                    m_pad_split_flag;
        MyPlane                 m_cut_plane;
        SegmentList             m_polyline;
        Polyhedron              m_left_shape;
        Polyhedron              m_right_shape;
        Polyhedron              m_pad_extrude_shape;
        Polyhedron              m_pad_min_shape;
        double                  m_minimun_surface;
        const Tooth&            m_tooth;
    public:
        Pad(const std::string& id, const Tooth& tooth);

        virtual ~Pad();

        void reset();

        void update();

        void save_to_db();
        void load_from_db();

        void set_cut_plane(Plane &cut_plane, Point &center, Point &right);

        void set_polylines(PointList &intersection);

        void set_left_shape(Polyhedron &left_shape);

        Polyhedron& get_left_shape() { return m_left_shape; };

        void set_right_shape(Polyhedron &right_shape);

        Polyhedron& get_right_shape() { return m_right_shape; };

        void set_pad_min_shape(Polyhedron &extrude_shape);

        Polyhedron& get_pad_min_shape() { return m_pad_min_shape; };

        // void set_extrude_shape(Polyhedron &extrude_shape);

        // Polyhedron& get_extrude_shape() { return m_pad_extrude_shape; };
    }; // end of class Pad
} // namespace OrthoLab

#endif
