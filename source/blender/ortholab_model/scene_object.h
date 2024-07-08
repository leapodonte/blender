#ifndef MY_SCENE_OBJECT_H
#define MY_SCENE_OBJECT_H

#include "types.h"
#include "object.h"

namespace OrthoLab
{
    class VisuRepBase;
    class Project;
    /**
     * class to represent a visible object in the scene**/
    class SceneObject : public Object
    {
    private: 
        bool                    m_is_updated = false;
        bool                    m_is_visible = false;
        Polyhedron              m_shape;
        std::vector<DataList>   m_shape_data; 
        VisuRepBase*            m_visu_rep = nullptr;
        const Project&          m_prj;
        const SceneObject&      m_parent;
        const std::string       m_id;
    public:
        SceneObject(const std::string& id, const Project& prj);
        SceneObject(const std::string& id, const SceneObject &parent);

        virtual ~SceneObject();

        virtual void reset();

        virtual void update();
        
        virtual void save_to_db(){};

        virtual void load_from_db(){};

        const std::string& get_id() const {return m_id;}

        bool is_updated() { return m_is_updated; }

        void set_updated(const bool is_updated) { m_is_updated = is_updated; }

        bool is_visible() { return m_is_visible; }

        void set_visible(const bool is_visible) { m_is_visible = is_visible; }

        Polyhedron& get_shape() { return m_shape; };

        void set_shape(Polyhedron& shape);

        std::vector<DataList>& get_shape_data() { return m_shape_data; };

        DataList& get_shape_data(const int index) { return m_shape_data[index]; };

        void set_shape_data(int index, const DataList &data) { m_shape_data[index] = data; };

        void add_shape_data(const DataList& data) { m_shape_data.push_back(data); };

        VisuRepBase* get_visu_rep();

        void set_visu_rep(VisuRepBase* visu_rep);

        const Project &get_project() const { return m_prj; };

    }; // end of class SceneObject
} // namespace OrthoLab

#endif
