#include "manu_asm_guide.h"

#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_mesh_processing/repair_polygon_soup.h>
#include <CGAL/convex_hull_3.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/alpha_wrap_3.h>
#include <CGAL/Polygon_mesh_processing/clip.h>

#include "teeth.h"
#include "tooth.h"
#include "extrusion.h"
#include "meshing.h"
#include "mesh_io.h"
#include "spdlog/spdlog.h"


namespace OrthoLab
{

    ManuAsmGuide::ManuAsmGuide(Teeth&  teeth) : 
        SceneObject(teeth.get_id(), teeth), 
        m_teeth(teeth)
    {
    }

    ManuAsmGuide::~ManuAsmGuide()
    {
        reset();
    }

    void ManuAsmGuide::reset()
    {
        SceneObject::reset(); 
    }


    void ManuAsmGuide::update2()
    {
        double guide_height = 20.0;
        double guide_cover_thickness = 1.0;
        double guide_bounding_box_max = 100.0;
        double guide_cover_overlap_height = 5.0;
        double alpha_wramp_alpha = 0.5;
        double alpha_wramp_offset = 0.05; 
        double guide_cover_hole_diameter = 4.0;
        
        Polyhedron* flatted_pads[10];
        Polyhedron* flatted_outlines[10];
        Vector pad_extrusion_direction[10];
        Polyhedron extruded_tooth[10];
		int total_pad = 0;
		int i = 0; int j = 0;
		//for (i = 0; i < 10;)
		//{
        // tooth 0
        flatted_outlines[i] = &m_teeth.get_tooth(j).get_pad_outline_poly(true);
        pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction();
		flatted_pads[i++] = &m_teeth.get_tooth(j++).get_pad_shape_poly(true);
        
        // tooth 1
        flatted_outlines[i] = &m_teeth.get_tooth(j).get_left_pad_outline_poly(true);
        pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction(1);
		flatted_pads[i++] = &m_teeth.get_tooth(j).get_left_pad_shape_poly(true);
        flatted_outlines[i] = &m_teeth.get_tooth(j).get_right_pad_outline_poly(true);
        pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction(2);
		flatted_pads[i++] = &m_teeth.get_tooth(j++).get_right_pad_shape_poly(true);

        // tooth 2
        flatted_outlines[i] = &m_teeth.get_tooth(j).get_left_pad_outline_poly(true);
        pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction(1);
		flatted_pads[i++] = &m_teeth.get_tooth(j).get_left_pad_shape_poly(true);
        flatted_outlines[i] = &m_teeth.get_tooth(j).get_right_pad_outline_poly(true);
        pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction(2);
		flatted_pads[i++] = &m_teeth.get_tooth(j++).get_right_pad_shape_poly(true);

        // tooth 3
        flatted_outlines[i] = &m_teeth.get_tooth(j).get_left_pad_outline_poly(true);
        pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction(1);
		flatted_pads[i++] = &m_teeth.get_tooth(j).get_left_pad_shape_poly(true);
        flatted_outlines[i] = &m_teeth.get_tooth(j).get_right_pad_outline_poly(true);
        pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction(2);
		flatted_pads[i++] = &m_teeth.get_tooth(j++).get_right_pad_shape_poly(true);

        // tooth 4
        flatted_outlines[i] = &m_teeth.get_tooth(j).get_left_pad_outline_poly(true);
        pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction(1);
		flatted_pads[i++] = &m_teeth.get_tooth(j).get_left_pad_shape_poly(true);
		if (flatted_pads[i - 1]->is_empty()) {
			// case tooth 4 is the final one
			flatted_outlines[i - 1] = &m_teeth.get_tooth(j).get_pad_outline_poly(true);
			pad_extrusion_direction[i - 1] = m_teeth.get_tooth(j).get_pad_extrusion_direction(1);
			flatted_pads[i - 1] = &m_teeth.get_tooth(j).get_pad_shape_poly(true);
            spdlog::trace("flatted_pads num{}", i);
		}
		else
		{
            flatted_outlines[i] = &m_teeth.get_tooth(j).get_right_pad_outline_poly(true);
            pad_extrusion_direction[i] = m_teeth.get_tooth(j).get_pad_extrusion_direction(2);
			flatted_pads[i++] = &m_teeth.get_tooth(j++).get_right_pad_shape_poly(true);

            // tooth 5
            flatted_outlines[i] = &m_teeth.get_tooth(j).get_pad_outline_poly(true);
			flatted_pads[i++] = &m_teeth.get_tooth(j++).get_pad_shape_poly(true);
            spdlog::trace("flatted_pads num{}", i);
		}
		//}
        total_pad = i;
        int tooth_count = total_pad;
        Point projection_center;
        Vector Projection_normal;
        Plane projection_plane;
        std::vector<Point> point_cloud, convex_hull_points;
        Vector teeth_direction[10];
        
        Vector teeth_centers[10];
        Point wire_extreem_point[6];
        double teeth_width[10];
        Vector global_direction(0, 0, 0);
        Vector teeth_center(0,0,0); 
        j = 0;
        Polyhedron tooth_outline_extrusion_collection;
        Polyhedron tooth_pad_collection;
        for (i = 0; i < tooth_count; i++)
        { 
            shell_extract_projection_situation(*flatted_pads[i], projection_plane, Projection_normal,
            projection_center);

            teeth_direction[i] = Projection_normal;

            point_cloud.push_back(projection_center);
            point_cloud.push_back(projection_center + guide_bounding_box_max*Projection_normal);
            Vector v(projection_center.x(), projection_center.y(), projection_center.z());
            teeth_center += v;
            teeth_centers[i] = v;
            global_direction += Projection_normal;
            for (auto v_it = flatted_pads[i]->vertices_begin(); v_it != flatted_pads[i]->vertices_end(); ++v_it)
            {
                convex_hull_points.push_back(v_it->point());
            }

            //if (i == 0 || i == tooth_count - 1)
            {
                // extrude tooths
                //Vector vec = m_teeth.get_tooth(i).get_pad_extrusion_direction();
                shell_extract_projection_situation(*flatted_outlines[i], projection_plane, Projection_normal,
                    projection_center);

                std::cout << "get_pad_extrusion_direction" << Projection_normal << std::endl;
                
                compute_extrusion(*flatted_outlines[i], extruded_tooth[i], -Projection_normal, 10);

				/*char buffer[10];
				_itoa(i, buffer, 10);
                std::string  filename = "extruded_tooth_pad_outline";
				std::string filename_i = filename + buffer + ".stl";
                mesh_export(extruded_tooth[i], filename_i);*/
                CGAL::Polygon_mesh_processing::corefine_and_compute_union(extruded_tooth[i], tooth_outline_extrusion_collection, tooth_outline_extrusion_collection);
                CGAL::Polygon_mesh_processing::corefine_and_compute_union(*flatted_pads[i], tooth_pad_collection, tooth_pad_collection);

			}
            //else if (i == tooth_count - 1)
            {

            }
        }          

        // teeth center
        Point teeth_center_point(teeth_center.x() / tooth_count, teeth_center.y() / tooth_count, teeth_center.z() / tooth_count);
        // global direction is the guide direction
        global_direction = global_direction / std::sqrt(global_direction.squared_length());
        // find the guide plan
        Plane guide_plane;
        CGAL::linear_least_squares_fitting_3(point_cloud.begin(), point_cloud.end(), guide_plane, CGAL::Dimension_tag<0>());
        // project the center point on the guide plan
        teeth_center_point = guide_plane.projection(teeth_center_point);
        // guide target point
        Point guide_target_point = teeth_center_point + global_direction;
        guide_target_point = guide_plane.projection(guide_target_point);

        Vector guide_direction = guide_target_point - teeth_center_point;
        guide_direction = guide_direction / std::sqrt(guide_direction.squared_length());

        // compute the guide top limit plan (withdraw 10mm to avoid collision with the teeth)
        Vector guid_plane_normal = guide_plane.orthogonal_vector();
        guid_plane_normal = guid_plane_normal / std::sqrt(guid_plane_normal.squared_length());

        Vector guid_left_direction = CGAL::cross_product(guid_plane_normal, guide_direction);
        guid_left_direction = guid_left_direction / std::sqrt(guid_left_direction.squared_length());

        Plane guide_top_limit_plane(teeth_center_point-guide_direction*guide_bounding_box_max, guide_direction);

        // compute the guide left limit plan
        Plane guide_left_limit_plane(teeth_center_point-guid_left_direction*guide_bounding_box_max, guid_left_direction);

        // compute the rear limit plan
        Plane guide_rear_limit_plane(teeth_center_point-guid_plane_normal*guide_bounding_box_max, guid_plane_normal);
        // Iterate over the vertices

        double max_distance_top = 0; double min_distance_top = 100000;
        double max_distance_left = 0; double min_distance_left = 100000;
        double max_distance_rear = 0; double min_distance_rear = 100000;
        for (int i = 0; i < tooth_count; i++)
        {
            double wire_max_distance_left = 0; double wire_min_distance_left = 100000;
            Point left_most_point; Point right_most_point;
            for (auto v = flatted_pads[i]->vertices_begin(); v != flatted_pads[i]->vertices_end(); ++v)
            {
                double distance_top = (CGAL::squared_distance(v->point(), guide_top_limit_plane));
                double distance_left = (CGAL::squared_distance(v->point(), guide_left_limit_plane));
                double distance_rear = (CGAL::squared_distance(v->point(), guide_rear_limit_plane));
                if (distance_top > max_distance_top)
                {
                    max_distance_top = distance_top;
                }
                if (distance_top < min_distance_top)
                {
                    min_distance_top = distance_top;
                }
                if (distance_left > max_distance_left)
                {
                    max_distance_left = distance_left;
                    right_most_point = v->point();
                }
                if(distance_left < min_distance_left) 
                {
                    min_distance_left = distance_left;
                    left_most_point = v->point();
                }
                if(distance_rear > max_distance_rear)
                {
                    max_distance_rear = distance_rear;
                }
                if(distance_rear < min_distance_rear)
                {
                    min_distance_rear = distance_rear;
                }
            }
            teeth_width[i]=CGAL::squared_distance(right_most_point, left_most_point);
            teeth_width[i]=sqrt(teeth_width[i]);
        }
        max_distance_top = sqrt(max_distance_top);
        min_distance_top = sqrt(min_distance_top);
        max_distance_left = sqrt(max_distance_left);
        min_distance_left = sqrt(min_distance_left);
        max_distance_rear = sqrt(max_distance_rear);
        min_distance_rear = sqrt(min_distance_rear);
        

        std::cout << "max_distance_top: " << max_distance_top << std::endl;
        std::cout << "min_distance_top: " << min_distance_top << std::endl;
        std::cout << "size_top: " << max_distance_top - min_distance_top << std::endl;
        std::cout << "max_distance_left: " << max_distance_left << std::endl;
        std::cout << "min_distance_left: " << min_distance_left << std::endl;
        std::cout << "size_left: " << max_distance_left - min_distance_left << std::endl;
        std::cout << "max_distance_rear: " << max_distance_rear << std::endl;
        std::cout << "min_distance_rear: " << min_distance_rear << std::endl;
        std::cout << "size_rear: " << max_distance_rear - min_distance_rear << std::endl;

        // find the majar direction of the teeth ==> left / right
        double bounding_box_left = std::max(abs(guide_bounding_box_max - max_distance_left), abs(guide_bounding_box_max - min_distance_left));
        double bounding_box_rear = std::max(abs(guide_bounding_box_max - max_distance_rear), abs(guide_bounding_box_max - min_distance_rear));
        
        Point guide_bottom_center = teeth_center_point + guide_direction * guide_height;
        Point guide_bottom_p1 = guide_bottom_center + guid_left_direction * bounding_box_left + guid_plane_normal * bounding_box_rear;
        Point guide_bottom_p2 = guide_bottom_center - guid_left_direction * bounding_box_left + guid_plane_normal * bounding_box_rear;
        Point guide_bottom_p3 = guide_bottom_center - guid_left_direction * bounding_box_left - guid_plane_normal * bounding_box_rear;
        Point guide_bottom_p4 = guide_bottom_center + guid_left_direction * bounding_box_left - guid_plane_normal * bounding_box_rear;
        convex_hull_points.push_back(guide_bottom_p1);
        convex_hull_points.push_back(guide_bottom_p2);
        convex_hull_points.push_back(guide_bottom_p3);
        convex_hull_points.push_back(guide_bottom_p4);

        Polyhedron final_convex_hull;
        CGAL::convex_hull_3(convex_hull_points.begin(), convex_hull_points.end(), final_convex_hull);    
        bool corefine_result;
        corefine_result = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(final_convex_hull, tooth_outline_extrusion_collection, final_convex_hull);
        if (corefine_result == false)
            std::cerr << "corefine_and_compute_difference error tooth_outline_extrusion_collection: " << std::endl;
        Polyhedron tooth_pad_collection_wrap;
        solve_selfintersection(tooth_pad_collection, tooth_pad_collection_wrap);
        corefine_result = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(final_convex_hull, tooth_pad_collection_wrap, final_convex_hull);
        if (corefine_result == false)
            std::cerr << "corefine_and_compute_difference error tooth_pad_collection: " << std::endl;

        mesh_export(tooth_pad_collection, "tooth_pad_collection.stl");
        mesh_export(final_convex_hull, "final_convex_hull_pads.stl");
        // Polyhedron wire_pad_union;
        // Polyhedron lowerpad[6];
        // filename = "upperpad";
        // if (lower_teeth) filename = "lowerpad";

        // std::vector<Point> teeth_wire_extreem_points;
        

        // // load the pads 
        // for (int i = start_index; i < end_index; i++)
        // {
        //     char buffer[10];
        //     _itoa(i, buffer, 10);
        //     std::string filename_i = filename + buffer + ".stl";
        //     if (!mesh_import(filename_i, lowerpad[i]))
        //         continue;

        //     if(lowerpad[i].empty()) continue;
        //     std::cout << "pad loaded: " << filename_i << std::endl;
        //     CGAL::Polygon_mesh_processing::corefine_and_compute_union(lowerpad[i], wire_pad_union, wire_pad_union);
        //     std::cout << "pad union: " << filename_i << std::endl;
        // }
        
        // // load the wires
        // Polyhedron lowerwire[5];
        // filename = "upperwire";
        // if (lower_teeth) filename = "lowerwire";

        // for (int i = start_index; i < end_index - 1; i++)
        // {
        //     double wire_max_distance_left = 0; double wire_min_distance_left = 100000;
        //     Point left_most_point; Point right_most_point;
        //     char buffer[10];
        //     _itoa(i, buffer, 10);
        //     std::string filename_i = filename + buffer + ".stl";
        //     if (!mesh_import(filename_i, lowerwire[i]))
        //         continue;
        //     if(lowerwire[i].empty()) continue;
        //     std::cout << "wire loaded: " << filename_i << std::endl;
        //     CGAL::alpha_wrap_3(lowerwire[i], alpha_wramp_alpha, alpha_wramp_offset, lowerwire[i]);
        //     std::cout << "wire wrapped: " << filename_i << std::endl;
        //     for (auto v_it = lowerwire[i].vertices_begin(); v_it != lowerwire[i].vertices_end(); ++v_it)
        //     {
        //         convex_hull_points.push_back(v_it->point());
        //         // find extreem points of wire
        //         double distance_left = (CGAL::squared_distance(v_it->point(), guide_left_limit_plane));
        //         if (distance_left > wire_max_distance_left)
        //         {
        //             wire_max_distance_left = distance_left;
        //             right_most_point = v_it->point();
        //         }
        //         if(distance_left < wire_min_distance_left) 
        //         {
        //             wire_min_distance_left = distance_left;
        //             left_most_point = v_it->point();
        //         }
        //     }
        //     wire_extreem_point[i] = left_most_point;
        //     wire_extreem_point[i+1] = right_most_point;
        //     max_distance_left = sqrt(max_distance_left);
        //     min_distance_left = sqrt(min_distance_left);

        //     CGAL::Polygon_mesh_processing::corefine_and_compute_union(lowerwire[i], wire_pad_union, wire_pad_union);
        //     std::cout << "wire union: " << filename_i << std::endl;

        //     //            
               
        // }

        

        // Polyhedron final_convex_hull;
        // CGAL::convex_hull_3(convex_hull_points.begin(), convex_hull_points.end(), final_convex_hull);
        // mesh_export(final_convex_hull, "final_convex_hull.stl");
        // // Polyhedron& closed_teeth = m_teeth.get_closed_teeth_poly(); 
        // if (lower_teeth)
        //     filename = "lowerteeth.stl";
        // else
        //     filename = "upperteeth.stl";
        // m_teeth.load_teeth_mesh(filename);
        // std::cout << "teeth model loaded: " << filename << std::endl;
        // Polyhedron& closed_teeth = m_teeth.get_closed_teeth_poly();
        // std::cout << "teeth model closed: " << filename << std::endl;
        // mesh_export(closed_teeth, "closed_teeth.stl");
        // Polyhedron final_convex_hull_sub_teeth;
        // bool result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(final_convex_hull, closed_teeth, final_convex_hull_sub_teeth);
        // solve_selfintersection(final_convex_hull_sub_teeth,final_convex_hull_sub_teeth);
        // mesh_export(final_convex_hull_sub_teeth, "final_convex_hull_sub_teeth.stl");
        // std::cout << "exported mesh: " << "final_convex_hull_sub_teeth.stl" << std::endl;
        // // 3. substract the teeth from this convex hull
        // Polyhedron final_convex_hull_sub_teeth_sub_pads;
        // solve_selfintersection(wire_pad_union,wire_pad_union);
        // result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(final_convex_hull_sub_teeth, wire_pad_union, final_convex_hull_sub_teeth_sub_pads);
        // mesh_export(final_convex_hull_sub_teeth_sub_pads, "final_convex_hull_sub_teeth_sub_pads.stl");
        // std::cout << "exported mesh: " << "final_convex_hull_sub_teeth_sub_pads.stl" << std::endl;
        // // 4 construct the guide cover
        // // move final_convex_hull with to the upper direction with 
        // // for (Polyhedron::Vertex_iterator it = final_convex_hull.vertices_begin(); it != final_convex_hull.vertices_end(); ++it)
        // // {
        // //     it->point() = it->point() - guide_direction * guide_cover_thickness;
        // // }
        // // expand the convex hull with the guide cover thickness  
        // alpha_wrap_mesh(final_convex_hull,final_convex_hull);
        // Polyhedron expanded_hull;
        // compute_enlarge_msh(final_convex_hull, expanded_hull, guide_cover_thickness);
        // mesh_export(expanded_hull, "final_convex_hull_expanded.stl");
        // solve_selfintersection(expanded_hull,expanded_hull);
        
        // Polyhedron final_convex_cover;
        // result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(expanded_hull, final_convex_hull_sub_teeth, final_convex_cover);
        // mesh_export(final_convex_cover, "final_convex_cover.stl");
        // std::cout << "exported mesh: " << "final_convex_cover.stl" << std::endl;
        // // clip the guide cover with plan
        // Plane guide_cover_clip_plane(teeth_center_point + 
        // guide_direction * (guide_bounding_box_max - min_distance_top + guide_cover_overlap_height), guide_direction);
        
        // solve_selfintersection(final_convex_cover,final_convex_cover);
        // // CGAL::alpha_wrap_3(final_convex_cover, alpha_wramp_alpha, alpha_wramp_offset, final_convex_cover);

        // CGAL::Polygon_mesh_processing::clip(final_convex_cover, guide_cover_clip_plane , CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
        // hole_filling(final_convex_cover, false);
        // mesh_export(final_convex_cover, "final_convex_cover_clip.stl");
        // std::cout << "exported mesh: " << "final_convex_cover_clip.stl" << std::endl;

        // // compute cone of intervention to be extracted from the cover
        // //Vector teeth_direction[6];
        // //Point wire_extreem_point[6];
        // Polyhedron  cones_collection;
        // for (int i = start_index; i < end_index; i++)
        // {
        //     Vector tooth_dir ((teeth_centers[i]).x()-wire_extreem_point[i].x(),teeth_centers[i].y()-wire_extreem_point[i].y(),
        //                         teeth_centers[i].z()-wire_extreem_point[i].z());
        //     Vector dir2 = CGAL::cross_product(tooth_dir, teeth_direction[i]);
        //     Polyhedron cone;
        //     construct_cone(cone, wire_extreem_point[i], dir2, -teeth_direction[i], -teeth_direction[i], guide_cover_hole_diameter*.5, 100, 20);
        //     CGAL::Polygon_mesh_processing::corefine_and_compute_union(cone, cones_collection, cones_collection);

        //     char buffer[10];
        //     _itoa(i, buffer, 10);
        //     std::string filename("final_convex_cover_holes_cone");
        //     std::string filename_i = filename + buffer + ".stl";
        //     mesh_export(cone, filename_i);
        // // 
        // }
        // alpha_wrap_mesh(cones_collection,cones_collection);
        // result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(final_convex_cover, cones_collection, final_convex_cover);
        // mesh_export(final_convex_cover, "final_convex_cover_with_holes.stl");
        // 
        // 1. generate a convex hull that contains all the teeth pads
        // Your polyhedra
        // Polyhedron poly1, poly2; // Add more as needed
        // fill_polyhedron(poly1);
        // fill_polyhedron(poly2);

        // // Extract vertices from each polyhedron and add them to the points vector
        // for (Polyhedron::Vertex_iterator v = poly1.vertices_begin(); v != poly1.vertices_end(); ++v)
        // {
        //     points.push_back(v->point());
        // }
        // for (Polyhedron::Vertex_iterator v = poly2.vertices_begin(); v != poly2.vertices_end(); ++v)
        // {
        //     points.push_back(v->point());
        // }

        // // Prepare an object to hold the convex hull
        // Polyhedron convex_hull;

        // // Compute the convex hull of the points
        // CGAL::convex_hull_3(points.begin(), points.end(), convex_hull);

        // 2. add some marges to this convex hull

        // 3. substract the teeth from this convex hull

        // 4. substract the pads from this convex hull

        // 5. substract the wire "tunnel" from this convex hull
    }

    void ManuAsmGuide::update(bool lower_teeth, const int start_index, const int end_index, 
                              const bool flattened)
    {
        if(flattened) return update2();

        double guide_height = 20.0;
        double guide_cover_thickness = 1.0;
        double guide_bounding_box_max = 100.0;
        double guide_cover_overlap_height = 5.0;
        double alpha_wramp_alpha = 0.5;
        double alpha_wramp_offset = 0.05; 
        double guide_cover_hole_diameter = 4.0;
        SceneObject::update();
        //0. prepare data
        Polyhedron lowertooth[6];
        Vector teeth_direction[6];
        Vector teeth_centers[6];
        Point wire_extreem_point[6];
        double teeth_width[6];
        std::string filename = "uppertooth";
        if (lower_teeth) filename = "lowertooth";
        std::string file_extension = ".stl";
        if (end_index) file_extension = "_t.stl";
        Point projection_center; Vector Projection_normal; Plane projection_plane;
        std::vector<Point> point_cloud, convex_hull_points;

        Vector global_direction(0,0,0);         
        Vector teeth_center(0,0,0);
        const int tooth_count = end_index - start_index;
        for (int i = start_index; i < end_index; i++)
        { 
            char buffer[10];
            _itoa(i, buffer, 10);
            std::string filename_i = filename + buffer + ".stl";
            if (!mesh_import(filename_i, lowertooth[i]))
                continue;
            if(lowertooth[i].empty()) continue;

            shell_extract_projection_situation(lowertooth[i], projection_plane, Projection_normal,
            projection_center);

            teeth_direction[i] = Projection_normal;

            point_cloud.push_back(projection_center);
            point_cloud.push_back(projection_center + guide_bounding_box_max*Projection_normal);
            Vector v(projection_center.x(), projection_center.y(), projection_center.z());
            teeth_center += v;
            teeth_centers[i] = v;
            global_direction += Projection_normal;
            for (auto v_it = lowertooth[i].vertices_begin(); v_it != lowertooth[i].vertices_end(); ++v_it)
            {
                convex_hull_points.push_back(v_it->point());
            }
        }
        // // temp dump
        // Polyhedron point_cloud_convex_hull;
        // CGAL::convex_hull_3(point_cloud.begin(), point_cloud.end(), point_cloud_convex_hull);
        // mesh_export(point_cloud_convex_hull, "point_cloud_convex_hull.stl");
        // // temp dump

        // teeth center
        Point teeth_center_point(teeth_center.x() / tooth_count, teeth_center.y() / tooth_count, teeth_center.z() / tooth_count);
        // global direction is the guide direction
        global_direction = global_direction / std::sqrt(global_direction.squared_length());
        // find the guide plan
        Plane guide_plane;
        CGAL::linear_least_squares_fitting_3(point_cloud.begin(), point_cloud.end(), guide_plane, CGAL::Dimension_tag<0>());
        // project the center point on the guide plan
        teeth_center_point = guide_plane.projection(teeth_center_point);
        // guide target point
        Point guide_target_point = teeth_center_point + global_direction;
        guide_target_point = guide_plane.projection(guide_target_point);

        Vector guide_direction = guide_target_point - teeth_center_point;
        guide_direction = guide_direction / std::sqrt(guide_direction.squared_length());

        // compute the guide top limit plan (withdraw 10mm to avoid collision with the teeth)
        Vector guid_plane_normal = guide_plane.orthogonal_vector();
        guid_plane_normal = guid_plane_normal / std::sqrt(guid_plane_normal.squared_length());

        Vector guid_left_direction = CGAL::cross_product(guid_plane_normal, guide_direction);
        guid_left_direction = guid_left_direction / std::sqrt(guid_left_direction.squared_length());

        Plane guide_top_limit_plane(teeth_center_point-guide_direction*guide_bounding_box_max, guide_direction);

        // compute the guide left limit plan
        Plane guide_left_limit_plane(teeth_center_point-guid_left_direction*guide_bounding_box_max, guid_left_direction);

        // compute the rear limit plan
        Plane guide_rear_limit_plane(teeth_center_point-guid_plane_normal*guide_bounding_box_max, guid_plane_normal);
        // Iterate over the vertices

        double max_distance_top = 0; double min_distance_top = 100000;
        double max_distance_left = 0; double min_distance_left = 100000;
        double max_distance_rear = 0; double min_distance_rear = 100000;
        for (int i = start_index; i < end_index; i++)
        {
            double wire_max_distance_left = 0; double wire_min_distance_left = 100000;
            Point left_most_point; Point right_most_point;
            for (auto v = lowertooth[i].vertices_begin(); v != lowertooth[i].vertices_end(); ++v)
            {
                double distance_top = (CGAL::squared_distance(v->point(), guide_top_limit_plane));
                double distance_left = (CGAL::squared_distance(v->point(), guide_left_limit_plane));
                double distance_rear = (CGAL::squared_distance(v->point(), guide_rear_limit_plane));
                if (distance_top > max_distance_top)
                {
                    max_distance_top = distance_top;
                }
                if (distance_top < min_distance_top)
                {
                    min_distance_top = distance_top;
                }
                if (distance_left > max_distance_left)
                {
                    max_distance_left = distance_left;
                    right_most_point = v->point();
                }
                if(distance_left < min_distance_left) 
                {
                    min_distance_left = distance_left;
                    left_most_point = v->point();
                }
                if(distance_rear > max_distance_rear)
                {
                    max_distance_rear = distance_rear;
                }
                if(distance_rear < min_distance_rear)
                {
                    min_distance_rear = distance_rear;
                }
            }
            teeth_width[i]=CGAL::squared_distance(right_most_point, left_most_point);
            teeth_width[i]=sqrt(teeth_width[i]);
        }
        max_distance_top = sqrt(max_distance_top);
        min_distance_top = sqrt(min_distance_top);
        max_distance_left = sqrt(max_distance_left);
        min_distance_left = sqrt(min_distance_left);
        max_distance_rear = sqrt(max_distance_rear);
        min_distance_rear = sqrt(min_distance_rear);
        

        std::cout << "max_distance_top: " << max_distance_top << std::endl;
        std::cout << "min_distance_top: " << min_distance_top << std::endl;
        std::cout << "size_top: " << max_distance_top - min_distance_top << std::endl;
        std::cout << "max_distance_left: " << max_distance_left << std::endl;
        std::cout << "min_distance_left: " << min_distance_left << std::endl;
        std::cout << "size_left: " << max_distance_left - min_distance_left << std::endl;
        std::cout << "max_distance_rear: " << max_distance_rear << std::endl;
        std::cout << "min_distance_rear: " << min_distance_rear << std::endl;
        std::cout << "size_rear: " << max_distance_rear - min_distance_rear << std::endl;

        // find the majar direction of the teeth ==> left / right
        double bounding_box_left = std::max(abs(guide_bounding_box_max - max_distance_left), abs(guide_bounding_box_max - min_distance_left));
        double bounding_box_rear = std::max(abs(guide_bounding_box_max - max_distance_rear), abs(guide_bounding_box_max - min_distance_rear));
        
        Point guide_bottom_center = teeth_center_point + guide_direction * guide_height;
        Point guide_bottom_p1 = guide_bottom_center + guid_left_direction * bounding_box_left + guid_plane_normal * bounding_box_rear;
        Point guide_bottom_p2 = guide_bottom_center - guid_left_direction * bounding_box_left + guid_plane_normal * bounding_box_rear;
        Point guide_bottom_p3 = guide_bottom_center - guid_left_direction * bounding_box_left - guid_plane_normal * bounding_box_rear;
        Point guide_bottom_p4 = guide_bottom_center + guid_left_direction * bounding_box_left - guid_plane_normal * bounding_box_rear;
        convex_hull_points.push_back(guide_bottom_p1);
        convex_hull_points.push_back(guide_bottom_p2);
        convex_hull_points.push_back(guide_bottom_p3);
        convex_hull_points.push_back(guide_bottom_p4);

        Polyhedron wire_pad_union;
        Polyhedron lowerpad[6];
        filename = "upperpad";
        if (lower_teeth) filename = "lowerpad";

        std::vector<Point> teeth_wire_extreem_points;
        

        // load the pads 
        for (int i = start_index; i < end_index; i++)
        {
            char buffer[10];
            _itoa(i, buffer, 10);
            std::string filename_i = filename + buffer + ".stl";
            if (!mesh_import(filename_i, lowerpad[i]))
                continue;

            if(lowerpad[i].empty()) continue;
            std::cout << "pad loaded: " << filename_i << std::endl;
            CGAL::Polygon_mesh_processing::corefine_and_compute_union(lowerpad[i], wire_pad_union, wire_pad_union);
            std::cout << "pad union: " << filename_i << std::endl;
        }
        
        // load the wires
        Polyhedron lowerwire[5];
        filename = "upperwire";
        if (lower_teeth) filename = "lowerwire";

        for (int i = start_index; i < end_index - 1; i++)
        {
            double wire_max_distance_left = 0; double wire_min_distance_left = 100000;
            Point left_most_point; Point right_most_point;
            char buffer[10];
            _itoa(i, buffer, 10);
            std::string filename_i = filename + buffer + ".stl";
            if (!mesh_import(filename_i, lowerwire[i]))
                continue;
            if(lowerwire[i].empty()) continue;
            std::cout << "wire loaded: " << filename_i << std::endl;
            CGAL::alpha_wrap_3(lowerwire[i], alpha_wramp_alpha, alpha_wramp_offset, lowerwire[i]);
            std::cout << "wire wrapped: " << filename_i << std::endl;
            for (auto v_it = lowerwire[i].vertices_begin(); v_it != lowerwire[i].vertices_end(); ++v_it)
            {
                convex_hull_points.push_back(v_it->point());
                // find extreem points of wire
                double distance_left = (CGAL::squared_distance(v_it->point(), guide_left_limit_plane));
                if (distance_left > wire_max_distance_left)
                {
                    wire_max_distance_left = distance_left;
                    right_most_point = v_it->point();
                }
                if(distance_left < wire_min_distance_left) 
                {
                    wire_min_distance_left = distance_left;
                    left_most_point = v_it->point();
                }
            }
            wire_extreem_point[i] = left_most_point;
            wire_extreem_point[i+1] = right_most_point;
            max_distance_left = sqrt(max_distance_left);
            min_distance_left = sqrt(min_distance_left);

            CGAL::Polygon_mesh_processing::corefine_and_compute_union(lowerwire[i], wire_pad_union, wire_pad_union);
            std::cout << "wire union: " << filename_i << std::endl;

            //            
               
        }

        

        Polyhedron final_convex_hull;
        CGAL::convex_hull_3(convex_hull_points.begin(), convex_hull_points.end(), final_convex_hull);
        mesh_export(final_convex_hull, "final_convex_hull.stl");
        // Polyhedron& closed_teeth = m_teeth.get_closed_teeth_poly(); 
        if (lower_teeth)
            filename = "lowerteeth.stl";
        else
            filename = "upperteeth.stl";
        m_teeth.load_teeth_mesh(filename);
        std::cout << "teeth model loaded: " << filename << std::endl;
        Polyhedron& closed_teeth = m_teeth.get_closed_teeth_poly();
        std::cout << "teeth model closed: " << filename << std::endl;
        mesh_export(closed_teeth, "closed_teeth.stl");
        Polyhedron final_convex_hull_sub_teeth;
        bool result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(final_convex_hull, closed_teeth, final_convex_hull_sub_teeth);
        solve_selfintersection(final_convex_hull_sub_teeth,final_convex_hull_sub_teeth);
        mesh_export(final_convex_hull_sub_teeth, "final_convex_hull_sub_teeth.stl");
        std::cout << "exported mesh: " << "final_convex_hull_sub_teeth.stl" << std::endl;
        // 3. substract the teeth from this convex hull
        Polyhedron final_convex_hull_sub_teeth_sub_pads;
        solve_selfintersection(wire_pad_union,wire_pad_union);
        result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(final_convex_hull_sub_teeth, wire_pad_union, final_convex_hull_sub_teeth_sub_pads);
        mesh_export(final_convex_hull_sub_teeth_sub_pads, "final_convex_hull_sub_teeth_sub_pads.stl");
        std::cout << "exported mesh: " << "final_convex_hull_sub_teeth_sub_pads.stl" << std::endl;
        // 4 construct the guide cover
        // move final_convex_hull with to the upper direction with 
        // for (Polyhedron::Vertex_iterator it = final_convex_hull.vertices_begin(); it != final_convex_hull.vertices_end(); ++it)
        // {
        //     it->point() = it->point() - guide_direction * guide_cover_thickness;
        // }
        // expand the convex hull with the guide cover thickness  
        alpha_wrap_mesh(final_convex_hull,final_convex_hull);
        Polyhedron expanded_hull;
        compute_enlarge_msh(final_convex_hull, expanded_hull, guide_cover_thickness);
        mesh_export(expanded_hull, "final_convex_hull_expanded.stl");
        solve_selfintersection(expanded_hull,expanded_hull);
        
        Polyhedron final_convex_cover;
        result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(expanded_hull, final_convex_hull_sub_teeth, final_convex_cover);
        mesh_export(final_convex_cover, "final_convex_cover.stl");
        std::cout << "exported mesh: " << "final_convex_cover.stl" << std::endl;
        // clip the guide cover with plan
        Plane guide_cover_clip_plane(teeth_center_point + 
        guide_direction * (guide_bounding_box_max - min_distance_top + guide_cover_overlap_height), guide_direction);
        
        solve_selfintersection(final_convex_cover,final_convex_cover);
        // CGAL::alpha_wrap_3(final_convex_cover, alpha_wramp_alpha, alpha_wramp_offset, final_convex_cover);

        CGAL::Polygon_mesh_processing::clip(final_convex_cover, guide_cover_clip_plane , CGAL::Polygon_mesh_processing::parameters::clip_volume(true));
        hole_filling(final_convex_cover, false);
        mesh_export(final_convex_cover, "final_convex_cover_clip.stl");
        std::cout << "exported mesh: " << "final_convex_cover_clip.stl" << std::endl;

        // compute cone of intervention to be extracted from the cover
        //Vector teeth_direction[6];
        //Point wire_extreem_point[6];
        Polyhedron  cones_collection;
        for (int i = start_index; i < end_index; i++)
        {
            Vector tooth_dir ((teeth_centers[i]).x()-wire_extreem_point[i].x(),teeth_centers[i].y()-wire_extreem_point[i].y(),
                                teeth_centers[i].z()-wire_extreem_point[i].z());
            Vector dir2 = CGAL::cross_product(tooth_dir, teeth_direction[i]);
            Polyhedron cone;
            construct_cone(cone, wire_extreem_point[i], dir2, -teeth_direction[i], -teeth_direction[i], guide_cover_hole_diameter*.5, 100, 20);
            CGAL::Polygon_mesh_processing::corefine_and_compute_union(cone, cones_collection, cones_collection);

            char buffer[10];
            _itoa(i, buffer, 10);
            std::string filename("final_convex_cover_holes_cone");
            std::string filename_i = filename + buffer + ".stl";
            mesh_export(cone, filename_i);
        // 
        }
        alpha_wrap_mesh(cones_collection,cones_collection);
        result_code = CGAL::Polygon_mesh_processing::corefine_and_compute_difference(final_convex_cover, cones_collection, final_convex_cover);
        mesh_export(final_convex_cover, "final_convex_cover_with_holes.stl");
        // 
        // 1. generate a convex hull that contains all the teeth pads
        // Your polyhedra
        // Polyhedron poly1, poly2; // Add more as needed
        // fill_polyhedron(poly1);
        // fill_polyhedron(poly2);

        // // Extract vertices from each polyhedron and add them to the points vector
        // for (Polyhedron::Vertex_iterator v = poly1.vertices_begin(); v != poly1.vertices_end(); ++v)
        // {
        //     points.push_back(v->point());
        // }
        // for (Polyhedron::Vertex_iterator v = poly2.vertices_begin(); v != poly2.vertices_end(); ++v)
        // {
        //     points.push_back(v->point());
        // }

        // // Prepare an object to hold the convex hull
        // Polyhedron convex_hull;

        // // Compute the convex hull of the points
        // CGAL::convex_hull_3(points.begin(), points.end(), convex_hull);

        // 2. add some marges to this convex hull

        // 3. substract the teeth from this convex hull

        // 4. substract the pads from this convex hull

        // 5. substract the wire "tunnel" from this convex hull
    }
} // namespace OrthoLab
