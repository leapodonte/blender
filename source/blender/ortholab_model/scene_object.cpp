#include "scene_object.h"
#include "types.h"
// #include "meshing.h"
#include "visurep_base.h"
#include "project.h"


namespace OrthoLab
{
    SceneObject::SceneObject(const std::string& id, const Project& prj):
    m_id(id),
    m_prj(prj),
    m_parent(*this)
    {
    }
    SceneObject::SceneObject(const std::string& id, const SceneObject &parent) : 
    m_id(id),
    m_prj(parent.get_project()),
    m_parent(parent)
    {
    }

    SceneObject::~SceneObject()
    {
        reset();
    }

    void SceneObject::reset()
    {
        Object::reset();
        m_shape.clear();
        m_shape_data.clear();
        m_is_updated = false;
    }

    void SceneObject::update()
    {
        set_updated(true);
    };

    void SceneObject::set_shape(Polyhedron& shape)
    {
        copy_polyehdral_surface(shape, m_shape);
    }

    VisuRepBase *SceneObject::get_visu_rep()
    {
        return m_visu_rep;
    }

    void SceneObject::set_visu_rep(VisuRepBase* visu_rep)
    {
        m_visu_rep = visu_rep;
    }
} // namespace OrthoLab
