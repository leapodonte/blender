#ifndef MY_MANU_ASM_GUIDE_H
#define MY_MANU_ASM_GUIDE_H

#include "types.h"
#include "scene_object.h"

namespace OrthoLab
{
    class Teeth;
    /**
     * class to represent a assembly guide for manufacturing**/
    class ManuAsmGuide : public SceneObject
    {
    private:
        Teeth& m_teeth;

    public:
        ManuAsmGuide(Teeth& teeth);

        virtual ~ManuAsmGuide();

        void reset();
        void update2();
        void update(bool lower_teeth, const int start_index=0, const int end_index=6, const bool flattened=true);
    }; // end of class ManuAsmGuide
} // namespace OrthoLab

#endif


