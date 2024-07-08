#ifndef RAMP_H
#define RAMP_H

namespace OrthoLab {

class Ramp
{
private:
    unsigned char   m_colors[4][256];
    int             m_nodes[256];
    double          min_value, max_value;
    bool            cyclic;
    bool            saturation;

public:
    // life cycle
    Ramp();

    ~Ramp();

    enum ColorMap
    {
        RAINBOW,
        THERMAL,
        FIRE,
        INFRARED,
        LIGHT_RAINBOW,
        MULTICOLOR,
        RANDOM_PASTEL,
        MENTHOL
    };

    unsigned char* color(unsigned int index);
    unsigned char r(unsigned int index) const;
    unsigned char g(unsigned int index) const;
    unsigned char b(unsigned int index) const;

    unsigned int range_index(double value) const;

    void get_color(const double value, double& mr, double& mg, double& mb);

    void set_cyclic(bool b);
    void set_saturation(bool b);
    void set_range(double minv, double maxv);
    void get_range(double& minv, double& maxv) const;

private:

    void rebuild();

    void reset();

public:

    void build_default();

    void add_node(unsigned int index,
        unsigned char r,
        unsigned char g,
        unsigned char b);

    void build_colormap(ColorMap colormap);

    void build_menthol();

    void build_fire();

    void build_rainbow();

    void build_infrared();

    void build_light_rainbow();

    void build_thermal();

    void build_multicolor();

    void build_random_pastel();
};

} // namespace OrthoLab


#endif
