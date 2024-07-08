#ifndef MY_OBJECT_H
#define MY_OBJECT_H

namespace OrthoLab
{
    /**
     * class to represent a object in ortholab name space**/
    class Object
    {
    private:
        bool m_init = false;

    public:
        Object(){}

        virtual ~Object(){reset();};

        virtual void reset() { m_init = false; }

        bool is_inited() const { return m_init; }

        void set_inited(bool init) { m_init = init; }
    }; // end of class Object
} // namespace OrthoLab

#endif
