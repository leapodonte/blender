#ifndef VISUREPBASE_H_
#define VISUREPBASE_H_

#include <memory>

namespace OrthoLab { 
	class SceneObject;
	class VisuRepBase
	{
	private:
		bool m_visible;
		bool are_buffers_initialized;
		SceneObject* m_object;

	public:
		VisuRepBase(SceneObject* object);
		virtual ~VisuRepBase();

		virtual void initialize_buffers();
		bool buffers_initialized() { return are_buffers_initialized; }
		void invalidate_buffers(){are_buffers_initialized = false;};

		SceneObject* getSceneObject()
		{
			return m_object;
		}
		void setVisible(bool visible)
		{
			m_visible = visible;
		}
		
		void toggle_view();

		bool isVisible()
		{
			return m_visible;
		}
	};
} // end namespace
#endif // VISUMANAGER_H_
