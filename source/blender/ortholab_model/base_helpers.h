#pragma once

#include <iostream>
#include <cstdlib>
#include <string>

namespace OrthoLab
{
#define BONDING_ZONE_DEPTH "BONDING_ZONE_DEPTH" // 0.15
#define BONDING_ZONE_MARGINE "BONDING_ZONE_MARGINE" // 0.25
#define BONDING_GRID_SIZE "BONDING_GRID_SIZE" // 0.5
#define BONDING_CONE_BASE_RADIUS "BONDING_CONE_BASE_RADIUS" // 0.15
#define BONDING_CONE_HEIGHT "BONDING_CONE_HEIGHT" // 1.5
#define BONDING_CONE_SEGMENTS "BONDING_CONE_SEGMENTS" // 20
#define BONDING_CONE_SURFACE_TELERENCE "BONDING_CONE_SURFACE_TELERENCE" // 0.05
#define BONDING_DISCRETIZE_STEP "BONDING_DISCRETIZE_STEP" // 10
#define WIRE_EXPANSION_PAD_RATIO "WIRE_EXPANSION_PAD_RATIO" // 2.0
#define PAD_SMOOTH_ITERATIONS "PAD_SMOOTH_ITERATIONS" // 10
#define EXTRUSION_ALPHA_WRAP_ALPHA_PARAM "EXTRUSION_ALPHA_WRAP_ALPHA_PARAM" // 200  
#define EXTRUSION_ALPHA_WRAP_OFFSET_PARAM "EXTRUSION_ALPHA_WRAP_OFFSET_PARAM"  // 2000

    float get_variable_float(const char *v_name, float default_value);

    bool set_env_variable(const std::string &name, const std::string &v);

    static int pad_geodesic_smooth_range = 10;
} // namespace OrthoLab
