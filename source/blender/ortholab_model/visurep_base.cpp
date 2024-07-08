#include "visurep_base.h"

namespace OrthoLab
{

	VisuRepBase::VisuRepBase(SceneObject* object) : m_object(object)
	{
	}
	VisuRepBase::~VisuRepBase(){

	};

	void VisuRepBase::initialize_buffers()
	{
		are_buffers_initialized = true;
	}
} // end namespace
