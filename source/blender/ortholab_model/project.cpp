#include "project.h"

#include "wire.h"
#include "save.h"
#include "timer.h" 
#include "file_services.h"
#include "kvstore.h"
#include "spdlog/spdlog.h"

namespace OrthoLab {
	inline bool ends_with(std::string const& value, std::string const& ending)
	{
		if (ending.size() > value.size()) return false;
		return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
	}

	MeshOptions Project::m_mesh_options;
	Project::Project(const std::string& case_name) :
		m_case_name(case_name),
		m_func(*this)
	{
		set_case_name(case_name);
		// std::cout << "Project address:" << this << std::endl;
		for (int i = 0; i < 2 * teeth_count; i++)
		{
			upper_seg_attribus.push_back(DataList());
			lower_seg_attribus.push_back(DataList());
		}

		for (int i = 0; i < teeth_count; i++)
		{
			upper_geodesic_attribus.push_back(DataList());
			lower_geodesic_attribus.push_back(DataList());
		}

		for (int i = 0; i < 2 * teeth_count; i++)
		{
			upper_pad_outlines_attribus.push_back(DataList());
			lower_pad_outlines_attribus.push_back(DataList());
		}

		for (int i = 0; i < teeth_count; i++)
		{
			upper_tooth_bounds_attribus.push_back(DataList());
			lower_tooth_bounds_attribus.push_back(DataList());
			upper_tooth_cut_points_attribus.push_back(DataList());
			lower_tooth_cut_points_attribus.push_back(DataList());
		}

		for (int i = 0; i < 2 * teeth_count; i++)
		{
			upper_tooth_cut_planes_attribus.push_back(DataList());
			lower_tooth_cut_planes_attribus.push_back(DataList());
		}

		for (int i = 0; i < 2 * pair_count; i++)
		{
			upper_pair_fitting_planes_attribus.push_back(DataList());
			lower_pair_fitting_planes_attribus.push_back(DataList());

			upper_pair_cut_planes_attribus.push_back(DataList());
			lower_pair_cut_planes_attribus.push_back(DataList());
		}

		for (int i = 0; i < pair_count; i++)
		{
			upper_pair_fitting_lines_attribus.push_back(DataList());
			lower_pair_fitting_lines_attribus.push_back(DataList());

			upper_pair_fitting_points_attribus.push_back(DataList());
			lower_pair_fitting_points_attribus.push_back(DataList());

			upper_pair_wires_attribus.push_back(DataList());
			lower_pair_wires_attribus.push_back(DataList());
		}

		for (int i = 0; i < 2 * pair_count; i++)
		{
			upper_wire_shapes_attribus.push_back(DataList());
			lower_wire_shapes_attribus.push_back(DataList());
			upper_wire_extrude_shapes_attribus.push_back(DataList());
			lower_wire_extrude_shapes_attribus.push_back(DataList());
		}

		for (int i = 0; i < 2 * teeth_count; i++)
		{
			upper_pad_shapes_attribus.push_back(DataList());
			lower_pad_shapes_attribus.push_back(DataList());
		}
	}

	Project::~Project()
	{
		if (m_kvstore) delete m_kvstore;
	}

	/******************* I/O **********************/


	void Project::save_to_db()
	{
#ifdef USE_OPENMP
#pragma omp parallel
		{
#pragma omp sections
			{
#pragma omp section
				{
#endif
					m_func.m_upper_teeth.save_to_db();
#ifdef USE_OPENMP
				}
#pragma omp section
				{
#endif
					m_func.m_lower_teeth.save_to_db();
#ifdef USE_OPENMP
				}
			}
		}
#endif
	}

	void Project::load_from_db()
	{
#ifdef USE_OPENMP
#pragma omp parallel
		{
#pragma omp sections
			{
#pragma omp section
				{
#endif
					m_func.m_upper_teeth.load_from_db();
#ifdef USE_OPENMP
				}
#pragma omp section
				{
#endif
					m_func.m_lower_teeth.load_from_db();
#ifdef USE_OPENMP
				}
			}
		}
#endif
		
	}


	bool Project::load_upper_teeth(const std::string& filename)
	{
		bool flag = false;

		if (ends_with(filename,".stl"))
			flag = m_func.load_upper_teeth_mesh(filename);
		else
			std::cout << "Format not supported!" << std::endl;

		if (flag)
			compute_elements(FlagMesh | FlagUpper);
		return flag;
	}

	bool Project::load_lower_teeth(const  std::string& filename)
	{
		bool flag = false;

		if (ends_with(filename, ".stl"))
			flag = m_func.load_lower_teeth_mesh(filename);
		else
			std::cout << "Format not supported!" << std::endl;

		if (flag)
			compute_elements(FlagMesh | FlagLower);
		return flag;
	}

	bool Project::load_upper_lower_teeth(const std::string &filename1,
										 const std::string &filename2)
	{
		bool flag = false;

		if (ends_with(filename1, ".stl") && ends_with(filename2, ".stl"))
		{
			flag = m_func.load_upper_lower_teeth(filename1, filename2);
			if (flag)
				compute_elements(FlagMesh | FlagLower | FlagUpper);
		}
		return flag;
	}
	void Project::import_segment(const std::string &filename, bool is_upper)
	{
		m_func.import_segment(filename, is_upper);
	}
	void Project::get_bounding_shape(double& x, double& y, double& z, double& r)
	{
		Point center = m_func.get_center();
		r = m_func.get_radius();
		x = center.x();
		y = center.y();
		z = center.z();
	}
	void Project::closeTeethMesh()
	{
		m_func.m_upper_teeth.close_teeth_mesh();
		m_func.m_lower_teeth.close_teeth_mesh();
	}

	void Project::expandTeethMesh(double dist)
	{
		m_func.m_upper_teeth.expandTeethMesh(dist);
		m_func.m_lower_teeth.expandTeethMesh(dist);
	}

	void Project::flattenWires(double dist)
	{
		m_func.m_upper_teeth.flattenWires(dist);
		m_func.m_lower_teeth.flattenWires(dist);
	}

	int Project::set_case_name(const std::string &case_name)
	{
		m_case_name = case_name;
		std::string db_path = ensure_directory_exists(m_case_name);

		if (db_path.length() == 0)
		{
			std::cerr << "Cannot open datastore path" << db_path << std::endl;
			return -1;
		}
		
		// Open a new database
		std::string db_file_name = db_path + PATH_SEPERATION + case_name + DB_EXT;
		KVStore* tmp_m_kvstore = new KVStore(db_file_name);

		if (nullptr == tmp_m_kvstore)
		{
			std::cerr << "Cannot open datastore" << db_file_name << std::endl;
			return -1;
		}
		else
		{
			if (m_kvstore) delete m_kvstore;
			m_kvstore = tmp_m_kvstore;
		}
		return 0;
	}
	/******************* Update **********************/
	void Project::set_mouse_pos(double x, double y, double z, double vx, double vy, double vz)
	{
		m_mouse_pos = OrthoLab::Point(x, y, z);
		m_mouse_vec = OrthoLab::Vector(vx, vy, vz);
	}

	bool Project::update_mouse_selection()
	{
		bool flag = m_func.update_selection(m_mouse_pos, m_mouse_vec, m_select_pos, selected_points);
		return flag;
		//if(flag)
			 //are_buffers_initialized = false;
	}

	void Project::add_selected_point(const double x, const double y, const double z)
	{
		bool flag = m_func.add_selected_point(x, y, z, m_select_pos, selected_points);
	}

	void Project::compute_elements(unsigned int flags)
	{
		clear_elements(flags);

		if ((flags & FlagMesh) && (flags & FlagUpper))
		{
			m_func.update_teeth_mesh(upper_teeth_facets, upper_teeth_normals, false);
		}

		if ((flags & FlagMesh) && (flags & FlagLower))
		{
			m_func.update_teeth_mesh(lower_teeth_facets, lower_teeth_normals, true);
		}

		if ((flags & FlagInside))
		{
			m_func.update_inside_function(upper_inside_ray_facets,
				upper_inside_ray_normals,
				upper_inside_ray_colors,
				upper_inside_func_colors,
				upper_inside_facets,
				upper_inside_normals,
				false);
			m_func.update_inside_function(lower_inside_ray_facets,
				lower_inside_ray_normals,
				lower_inside_ray_colors,
				lower_inside_func_colors,
				lower_inside_facets,
				lower_inside_normals,
				true);

		}

		if ((flags & FlagConvex))
		{
			m_func.update_convex_function(upper_local_convex_facets,
				upper_local_convex_normals,
				upper_local_convex_colors,
				upper_convex_func_colors,
				false);
			m_func.update_convex_function(lower_local_convex_facets,
				lower_local_convex_normals,
				lower_local_convex_colors,
				lower_convex_func_colors,
				true);
		}

		if ((flags & FlagConvexIso))
		{
			m_func.update_convexity_isovalue(upper_convex_iso_edges,
				m_mesh_options.convex_isovalue,
				false);
			m_func.update_convexity_isovalue(lower_convex_iso_edges,
				m_mesh_options.convex_isovalue,
				true);
		}

		if ((flags & FlagSegment))
		{
			m_func.update_segmentation(selected_seg_facets, selected_seg_normals);
		}

		if ((flags & FlagValid) && (flags & FlagUpper))
		{
			for (int i = 0; i < 6; i++)
				m_func.update_validated_segmentation(i, false, upper_seg_attribus[2 * i], upper_seg_attribus[2 * i + 1]);
		}

		if ((flags & FlagValid) && (flags & FlagLower))
		{
			for (int i = 0; i < 6; i++)
				m_func.update_validated_segmentation(i, true, lower_seg_attribus[2 * i], lower_seg_attribus[2 * i + 1]);
		}

		if ((flags & FlagPadGeodesic))
		{
			for (int i = 0; i < 6; i++)
			{
				if (m_mesh_options.upper_flags[i])
					m_func.update_tooth_geodesic_function(i, false, upper_geodesic_attribus[i]);

				if (m_mesh_options.lower_flags[i])
					m_func.update_tooth_geodesic_function(i, true, lower_geodesic_attribus[i]);
			}
		}

		if ((flags & FlagPadOutline))
		{
			for (int i = 0; i < 6; i++)
			{
				if (m_mesh_options.upper_flags[i])
					m_func.update_pad_outline(i, false, upper_pad_outlines_attribus[2 * i], upper_pad_outlines_attribus[2 * i + 1]);

				if (m_mesh_options.lower_flags[i])
					m_func.update_pad_outline(i, true, lower_pad_outlines_attribus[2 * i], lower_pad_outlines_attribus[2 * i + 1]);
			}
		}

		if ((flags & FlagToothCut))
		{
			for (int i = 0; i < teeth_count; i++)
			{
				if (m_mesh_options.upper_flags[i])
				{
					m_func.update_tooth_cut_plane(i, false,
						upper_tooth_bounds_attribus[i],
						upper_tooth_cut_points_attribus[i],
						upper_tooth_cut_planes_attribus[2 * i],
						upper_tooth_cut_planes_attribus[2 * i + 1]);
				}


				if (m_mesh_options.lower_flags[i])
				{
					m_func.update_tooth_cut_plane(i, true,
						lower_tooth_bounds_attribus[i],
						lower_tooth_cut_points_attribus[i],
						lower_tooth_cut_planes_attribus[2 * i],
						lower_tooth_cut_planes_attribus[2 * i + 1]);
				}
			}
		}

		if ((flags & FlagPairCut))
		{
			for (int i = 0; i < pair_count; i++)
			{
				if (m_mesh_options.upper_flags[i] && m_mesh_options.upper_flags[i + 1])
				{
					m_func.m_upper_teeth.get_wire(i).update_pair_cut_plane(
						upper_pair_fitting_planes_attribus[2 * i],
						upper_pair_fitting_planes_attribus[2 * i + 1],
						upper_pair_fitting_lines_attribus[i],
						upper_pair_fitting_points_attribus[i],
						upper_pair_cut_planes_attribus[2 * i],
						upper_pair_cut_planes_attribus[2 * i + 1],
						upper_pair_wires_attribus[i]);
				}

				if (m_mesh_options.lower_flags[i] && m_mesh_options.lower_flags[i + 1])
				{
					m_func.m_lower_teeth.get_wire(i).update_pair_cut_plane(
						lower_pair_fitting_planes_attribus[2 * i],
						lower_pair_fitting_planes_attribus[2 * i + 1],
						lower_pair_fitting_lines_attribus[i],
						lower_pair_fitting_points_attribus[i],
						lower_pair_cut_planes_attribus[2 * i],
						lower_pair_cut_planes_attribus[2 * i + 1],
						lower_pair_wires_attribus[i]);
				}
			}
		}

		if ((flags & FlagPairWire))
		{
			for (int i = 0; i < pair_count; i++)
			{
				if (m_mesh_options.upper_flags[i] && m_mesh_options.upper_flags[i + 1])
				{
					m_func.update_pair_wire_shape(i, false,
						upper_wire_shapes_attribus[2 * i],
						upper_wire_shapes_attribus[2 * i + 1],
						upper_wire_extrude_shapes_attribus[2 * i],
						upper_wire_extrude_shapes_attribus[2 * i + 1]);
				}

				if (m_mesh_options.lower_flags[i] && m_mesh_options.lower_flags[i + 1])
				{
					m_func.update_pair_wire_shape(i, true,
						lower_wire_shapes_attribus[2 * i],
						lower_wire_shapes_attribus[2 * i + 1],
						lower_wire_extrude_shapes_attribus[2 * i],
						lower_wire_extrude_shapes_attribus[2 * i + 1]);
				}
			}
		}

		if ((flags & FlagPadShape))
		{
			for (int i = 0; i < teeth_count; i++)
			{
				if (m_mesh_options.upper_flags[i])
				{
					m_func.update_pad_shape(i, false,
						upper_pad_shapes_attribus[2 * i],
						upper_pad_shapes_attribus[2 * i + 1]);
				}

				if (m_mesh_options.lower_flags[i])
				{
					m_func.update_pad_shape(i, true,
						lower_pad_shapes_attribus[2 * i],
						lower_pad_shapes_attribus[2 * i + 1]);
				}
			}
		}

		if ((flags & FlagTeethDist))
		{
			m_func.update_teeth_distance_colors(false, m_mesh_options.teeth_dist_thresh, upper_distance_colors);
			m_func.update_teeth_distance_colors(true, m_mesh_options.teeth_dist_thresh, lower_distance_colors);
		}

		//are_buffers_initialized = false;
	}

	void Project::clear_elements(unsigned int flags)
	{
		if ((flags & FlagMesh) && (flags & FlagUpper))
		{
			upper_teeth_facets.clear();
			upper_teeth_normals.clear();
		}

		if ((flags & FlagMesh) && (flags & FlagLower))
		{
			lower_teeth_facets.clear();
			lower_teeth_normals.clear();
		}

		if ((flags & FlagInside))
		{
			upper_inside_ray_facets.clear();
			upper_inside_ray_normals.clear();
			upper_inside_ray_colors.clear();
			upper_inside_func_colors.clear();
			upper_inside_facets.clear();
			upper_inside_normals.clear();
			lower_inside_ray_facets.clear();
			lower_inside_ray_normals.clear();
			lower_inside_ray_colors.clear();
			lower_inside_func_colors.clear();
			lower_inside_facets.clear();
			lower_inside_normals.clear();
		}

		if ((flags & FlagConvex))
		{
			upper_local_convex_facets.clear();
			upper_local_convex_normals.clear();
			upper_local_convex_colors.clear();
			upper_convex_func_colors.clear();
			lower_local_convex_facets.clear();
			lower_local_convex_normals.clear();
			lower_local_convex_colors.clear();
			lower_convex_func_colors.clear();

		}

		if ((flags & FlagConvexIso))
		{
			upper_convex_iso_edges.clear();
			lower_convex_iso_edges.clear();
		}

		if ((flags & FlagSegment))
		{
			selected_seg_facets.clear();
			selected_seg_normals.clear();
		}

		if ((flags & FlagValid) && (flags & FlagUpper))
		{
			for (int i = 0; i < upper_seg_attribus.size(); i++)
				upper_seg_attribus[i].clear();
		}

		if ((flags & FlagValid) && (flags & FlagLower))
		{
			for (int i = 0; i < lower_seg_attribus.size(); i++)
				lower_seg_attribus[i].clear();
		}

		if ((flags & FlagPadGeodesic))
		{
			for (int i = 0; i < 6; i++)
			{
				if (m_mesh_options.upper_flags[i])
				{
					upper_geodesic_attribus[i].clear();
				}

				if (m_mesh_options.lower_flags[i])
				{
					lower_geodesic_attribus[i].clear();
				}
			}
		}

		if ((flags & FlagPadOutline))
		{
			for (int i = 0; i < 6; i++)
			{
				if (m_mesh_options.upper_flags[i])
				{
					upper_pad_outlines_attribus[2 * i].clear();
					upper_pad_outlines_attribus[2 * i + 1].clear();
				}

				if (m_mesh_options.lower_flags[i])
				{
					lower_pad_outlines_attribus[2 * i].clear();
					lower_pad_outlines_attribus[2 * i + 1].clear();
				}
			}
		}

		if ((flags & FlagToothCut))
		{
			for (int i = 0; i < teeth_count; i++)
			{
				if (m_mesh_options.upper_flags[i])
				{
					upper_tooth_bounds_attribus[i].clear();
					upper_tooth_cut_points_attribus[i].clear();
					upper_tooth_cut_planes_attribus[2 * i + 0].clear();
					upper_tooth_cut_planes_attribus[2 * i + 1].clear();
				}

				if (m_mesh_options.lower_flags[i])
				{
					lower_tooth_bounds_attribus[i].clear();
					lower_tooth_cut_points_attribus[i].clear();
					lower_tooth_cut_planes_attribus[2 * i + 0].clear();
					lower_tooth_cut_planes_attribus[2 * i + 1].clear();
				}
			}
		}

		if ((flags & FlagPairCut))
		{
			for (int i = 0; i < pair_count; i++)
			{
				if (m_mesh_options.upper_flags[i] && m_mesh_options.upper_flags[i + 1])
				{
					upper_pair_fitting_planes_attribus[2 * i + 0].clear();
					upper_pair_fitting_planes_attribus[2 * i + 1].clear();
					upper_pair_fitting_lines_attribus[i].clear();
					upper_pair_fitting_points_attribus[i].clear();
					upper_pair_cut_planes_attribus[2 * i + 0].clear();
					upper_pair_cut_planes_attribus[2 * i + 1].clear();
					upper_pair_wires_attribus[i].clear();
				}

				if (m_mesh_options.lower_flags[i] && m_mesh_options.lower_flags[i + 1])
				{
					lower_pair_fitting_planes_attribus[2 * i + 0].clear();
					lower_pair_fitting_planes_attribus[2 * i + 1].clear();
					lower_pair_fitting_lines_attribus[i].clear();
					lower_pair_fitting_points_attribus[i].clear();
					lower_pair_cut_planes_attribus[2 * i + 0].clear();
					lower_pair_cut_planes_attribus[2 * i + 1].clear();
					lower_pair_wires_attribus[i].clear();
				}
			}
		}

		if ((flags & FlagPairWire))
		{
			for (int i = 0; i < pair_count; i++)
			{
				if (m_mesh_options.upper_flags[i] && m_mesh_options.upper_flags[i + 1])
				{
					upper_wire_shapes_attribus[2 * i + 0].clear();
					upper_wire_shapes_attribus[2 * i + 1].clear();
					upper_wire_extrude_shapes_attribus[2 * i + 0].clear();
					upper_wire_extrude_shapes_attribus[2 * i + 1].clear();
				}

				if (m_mesh_options.lower_flags[i] && m_mesh_options.lower_flags[i + 1])
				{
					lower_wire_shapes_attribus[2 * i + 0].clear();
					lower_wire_shapes_attribus[2 * i + 1].clear();
					lower_wire_extrude_shapes_attribus[2 * i + 0].clear();
					lower_wire_extrude_shapes_attribus[2 * i + 1].clear();
				}
			}
		}

		if ((flags & FlagPadShape))
		{
			for (int i = 0; i < teeth_count; i++)
			{
				if (m_mesh_options.upper_flags[i])
				{
					upper_pad_shapes_attribus[2 * i + 0].clear();
					upper_pad_shapes_attribus[2 * i + 1].clear();
				}

				if (m_mesh_options.lower_flags[i])
				{
					lower_pad_shapes_attribus[2 * i + 0].clear();
					lower_pad_shapes_attribus[2 * i + 1].clear();
				}
			}
		}

		//are_buffers_initialized = false;
	}

	void Project::clear_data()
	{
		m_func.reset();
		clear_elements(FlagAll);
	}

	/******************* Function **********************/

	void Project::compute_laplacian_based_inside()
	{
		spdlog::trace("compute_laplacian_based_inside");
		// m_func.compute_laplacian_based_inside(m_mesh_options.inside_smooth_range,
		// 	m_mesh_options.inside_laplacian_param,
		// 	m_mesh_options.inside_isovalue);
		// Function::compute_laplacian_based_inside(int smooth_range, double lambda, double value)
		int smooth_range = m_mesh_options.inside_smooth_range;
		double lambda = m_mesh_options.inside_laplacian_param;
		double value = m_mesh_options.inside_isovalue;
		double enlarge_dist = m_mesh_options.inside_enlarge_dist;

		OrthoLab::Timer timer;
		if (!m_func.m_upper_init || !m_func.m_lower_init)
		{
			std::cout << "Please load both upper teeth and lower teeth!";
			return;
		}

		// Find camera position for teeth
		PointList lower_verts;
		m_func.m_lower_teeth.read_polygon_vertices(lower_verts);
		Bbox bounding_box = CGAL::bbox_3(lower_verts.begin(), lower_verts.end());
		double x = 0.5 * (bounding_box.xmin() + bounding_box.xmax());
		double y = 0.5 * (bounding_box.ymin() + bounding_box.ymax());
		double z = 0.5 * (bounding_box.zmin() + bounding_box.zmax());

		Point upper_center(x, y, bounding_box.zmax());
		Point lower_center(x, y, z);

		bool flag_upper = false;
		bool flag_lower = false;

		// Compute enlarged teeth tree;
		AABBTree enlarge_tree;
		
		Polyhedron& upper_poly = m_func.m_upper_teeth.get_teeth_poly();
		enlarge_tree.insert(faces(upper_poly).first, faces(upper_poly).second, upper_poly);

		Polyhedron enlarge_lower_poly;
        Polyhedron& lower_poly = m_func.m_lower_teeth.get_teeth_poly();
        compute_enlarge_msh(lower_poly, enlarge_lower_poly, enlarge_dist);
		enlarge_tree.insert(faces(enlarge_lower_poly).first, faces(enlarge_lower_poly).second, enlarge_lower_poly);

		enlarge_tree.accelerate_distance_queries();

		// Solve constrained laplacian for inside segmentation
#ifdef USE_OPENMP
#pragma omp parallel
		{
#pragma omp sections
			{
#pragma omp section
				{
#endif
					std::cout << "Solving laplacian-based inside for upper teeth..." << std::endl;
					flag_upper = m_func.m_upper_teeth.compute_laplacian_based_inside(
						enlarge_tree, upper_center, smooth_range, lambda);
#ifdef USE_OPENMP
				}
#pragma omp section
				{
#endif
					std::cout << "Solving laplacian-based inside for lower teeth..." << std::endl;
					flag_lower = m_func.m_lower_teeth.compute_laplacian_based_inside(
						enlarge_tree, lower_center, smooth_range, lambda);
#ifdef USE_OPENMP
				}
			}
		}
#endif

		if (!flag_upper || !flag_lower)
		{
			std::cout << "Solver failed!" << std::endl;
			return;
		}
#ifdef USE_OPENMP
#pragma omp parallel
		{
#pragma omp sections
			{
#pragma omp section
				{
#endif
					// Extract isosurface for inside segmentation
					// std::cout << "Extracting inside surface for upper teeth..." << std::endl;
					m_func.m_upper_teeth.extract_inside_from_function(value);
#ifdef USE_OPENMP
				}
#pragma omp section
				{
#endif
					// std::cout << "Extracting inside surface for lower teeth..." << std::endl;
					m_func.m_lower_teeth.extract_inside_from_function(value);
#ifdef USE_OPENMP
				}
			}
		}
#endif
		// spdlog::trace("compute_laplacian_based_inside Elapsed time: {} seconds", timer.elapsed());
		std::cout << "compute_laplacian_based_inside Elapsed time: " << timer.elapsed() << " seconds" << std::endl;
		compute_elements(FlagInside);
	}

	void Project::compute_laplacian_based_convexity()
	{
		// m_func.compute_laplacian_based_convexity(m_mesh_options.convex_neighbor_radius,
		// 	m_mesh_options.convex_smooth_range,
		// 	m_mesh_options.convex_laplacian_param);
		//compute_laplacian_based_convexity(double radius_ratio, int smooth_range, double lambda)
		double radius_ratio = m_mesh_options.convex_neighbor_radius;
		int smooth_range = m_mesh_options.convex_smooth_range;
		double lambda = m_mesh_options.convex_laplacian_param;

		if (!m_func.m_upper_init || !m_func.m_lower_init)
		{
			std::cout << "Please load both upper teeth and lower teeth!";
			return;
		}
		bool flag_upper = false;
		bool flag_lower = false;
#ifdef USE_OPENMP
#pragma omp parallel
		{
#pragma omp sections
			{
#pragma omp section
				{
#endif
					// Solve constrained laplacian for convexity field
					std::cout << "Solving laplacian-based convexity field for upper teeth..." << std::endl;
					flag_upper = m_func.m_upper_teeth.compute_laplacian_based_convexity(
						radius_ratio, m_func.m_radius, smooth_range, lambda);
#ifdef USE_OPENMP
				}
#pragma omp section
				{
#endif
					std::cout << "Solving laplacian-based convexity field for lower teeth..." << std::endl;
					flag_lower = m_func.m_lower_teeth.compute_laplacian_based_convexity(
						radius_ratio, m_func.m_radius, smooth_range, lambda);
#ifdef USE_OPENMP
				}
			}
		}
#endif
		if (!flag_upper || !flag_lower)
		{
			std::cout << "Solver failed!" << std::endl;
			return;
		}

		compute_elements(FlagConvex | FlagConvexIso);
	}

	void Project::update_convexity_function_isovalue()
	{
		compute_elements(FlagConvexIso);
	}

	void Project::compute_selected_segmentation()
	{
		m_func.compute_selected_segmentation(m_select_pos,
			m_mesh_options.convex_isovalue);

		compute_elements(FlagSegment);
	}

	void Project::clear_selected_points()
	{
		m_select_pos.clear();
		selected_points.clear();
		//are_buffers_initialized = false;
	}

	void Project::validate_selected_segmentation()
	{
		int status = m_func.validate_selected_segmentation();

		if (status == 0)
			return;

		if (status == 1)
			compute_elements(FlagSegment | FlagValid | FlagUpper);
		else if (status == -1)
			compute_elements(FlagSegment | FlagValid | FlagLower);
	}

	/******************* test **********************/
	void Project::save_validated_segmentation()
	{
		//int status = m_func.validate_selected_segmentation();

		//if (status == 0)
		//	return;

		//if (status == 1)
		m_func.m_upper_teeth.save_validated_segmentation("upper");
		// else if (status == -1)
		m_func.m_lower_teeth.save_validated_segmentation("lower");
	}

	void Project::load_validated_segmentation()
	{
		m_func.m_upper_teeth.load_validated_segmentation("upper");
		m_func.m_lower_teeth.load_validated_segmentation("lower");
	}

	void Project::save_pad_statistics()
	{
		m_func.save_pad_statistics();
	}
	
	void Project::close_mesh(const std::string &pad, bool shell_solid)
	{
		m_func.close_mesh(pad, shell_solid);
	}

	void Project::expand_mesh(const std::string &pad,double thickness)
	{
		m_func.Expand_mesh(pad,thickness);
	}

	void Project::compute_geodesic_pad_outlines()
	{
		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.upper_flags[i])
			{
				double geodesic_isovalue = m_mesh_options.geodesic_isovalue;
				if(m_mesh_options.geodesic_isovalues[i]>0)geodesic_isovalue = m_mesh_options.geodesic_isovalues[i];
				m_func.compute_geodesic_pad_outlines(i, false, m_mesh_options.geodesic_step, geodesic_isovalue);
			}
		}

		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.lower_flags[i])
			{
				double geodesic_isovalue = m_mesh_options.geodesic_isovalue;
				if(m_mesh_options.geodesic_isovalues[i]>0)geodesic_isovalue = m_mesh_options.geodesic_isovalues[i];
				m_func.compute_geodesic_pad_outlines(i, true, m_mesh_options.geodesic_step, geodesic_isovalue);
			}
		}

		compute_elements(FlagPadGeodesic | FlagPadOutline);
	}

	void Project::compute_circular_pad_outlines()
	{
		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.upper_flags[i])
			{
				double geodesic_isovalue = m_mesh_options.geodesic_isovalue;
				if(m_mesh_options.geodesic_isovalues[i]>0)geodesic_isovalue = m_mesh_options.geodesic_isovalues[i];
				m_func.compute_circular_pad_outlines(i, false, m_mesh_options.geodesic_step, geodesic_isovalue, m_mesh_options.pad_minimim_area);
			}
		}

		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.lower_flags[i])
			{
				double geodesic_isovalue = m_mesh_options.geodesic_isovalue;
				if(m_mesh_options.geodesic_isovalues[i]>0)geodesic_isovalue = m_mesh_options.geodesic_isovalues[i];
				m_func.compute_circular_pad_outlines(i, true, m_mesh_options.geodesic_step, geodesic_isovalue, m_mesh_options.pad_minimim_area);
			}
		}

		compute_elements(FlagPadGeodesic | FlagPadOutline);
	}

	void Project::recompute_geodesic_pad_outlines()
	{
		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.upper_flags[i])
			{
				double geodesic_isovalue = m_mesh_options.geodesic_isovalue;
				if(m_mesh_options.geodesic_isovalues[i]>0)geodesic_isovalue = m_mesh_options.geodesic_isovalues[i];
				m_func.recompute_geodesic_pad_outlines(i, geodesic_isovalue, false);
			}
		}

		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.lower_flags[i])
			{
				double geodesic_isovalue = m_mesh_options.geodesic_isovalue;
				if(m_mesh_options.geodesic_isovalues[i]>0)geodesic_isovalue = m_mesh_options.geodesic_isovalues[i];
				m_func.recompute_geodesic_pad_outlines(i, geodesic_isovalue, true);
			}
		}

		compute_elements(FlagPadOutline);
	}

	void Project::compute_tooth_split_planes()
	{
		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.upper_flags[i])
			{
				m_func.compute_tooth_cut_plane(i, false);
			}
		}

		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.lower_flags[i])
			{
				m_func.compute_tooth_cut_plane(i, true);
			}
		}

		compute_elements(FlagToothCut);
	}

	void Project::compute_paired_wires()
	{
		for (int i = 0; i < pair_count; i++)
		{
			if (m_mesh_options.upper_flags[i] && m_mesh_options.upper_flags[i + 1])
				m_func.compute_paired_wires(i, false, m_mesh_options.pair_distances[i], m_mesh_options.pair_tols[i], m_mesh_options.pair_max_bends[i], m_mesh_options.pair_radius, m_mesh_options.pair_offsets[i]);

			if (m_mesh_options.lower_flags[i] && m_mesh_options.lower_flags[i + 1])
				m_func.compute_paired_wires(i, true, m_mesh_options.pair_distances[i], m_mesh_options.pair_tols[i], m_mesh_options.pair_max_bends[i], m_mesh_options.pair_radius, m_mesh_options.pair_offsets[i]);
		}

		compute_elements(FlagPairCut);
	}

	void Project::compute_coplanar_wires()
	{
		if(m_mesh_options.upper_flags[6])
			m_func.compute_coplanar_wires(false, m_mesh_options.pair_distances[0]);

		if(m_mesh_options.lower_flags[6])
			m_func.compute_coplanar_wires(true, m_mesh_options.pair_distances[0]);

		compute_elements(FlagPairCut);
	}

	void Project::compute_wire_shapes()
	{
		#pragma omp parallel for
		for (int i = 0; i < pair_count; i++)
		{
			if (m_mesh_options.upper_flags[i] && m_mesh_options.upper_flags[i + 1])
			{
				m_func.compute_wire_shapes(i, false, m_mesh_options.geodesic_step, m_mesh_options.pair_radius, m_mesh_options.pair_enlarge_ratio);
			}
				

			if (m_mesh_options.lower_flags[i] && m_mesh_options.lower_flags[i + 1])
			{
				m_func.compute_wire_shapes(i, true, m_mesh_options.geodesic_step, m_mesh_options.pair_radius, m_mesh_options.pair_enlarge_ratio);
			}
				
		}

		compute_elements(FlagPairWire);
	}

	void Project::compute_pad_shapes()
	{
#pragma omp parallel for
		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.upper_flags[i])
			{
				bool flag_left = (i > 0) ? m_mesh_options.upper_flags[i-1] : false;
				bool flag_right = (i < teeth_count-1) ? m_mesh_options.upper_flags[i+1] : false;
				m_func.m_upper_teeth.compute_pad_shapes(i, flag_left, flag_right, m_mesh_options.pad_height, m_mesh_options.pad_wire_height, m_mesh_options.pair_radius);
			}
			if (m_mesh_options.lower_flags[i])
			{
				bool flag_left = (i > 0) ? m_mesh_options.lower_flags[i-1] : false;
				bool flag_right = (i < teeth_count-1) ? m_mesh_options.lower_flags[i+1] : false;
				m_func.m_lower_teeth.compute_pad_shapes(i, flag_left, flag_right, m_mesh_options.pad_height, m_mesh_options.pad_wire_height, m_mesh_options.pair_radius);
			}
		}

		compute_elements(FlagPadShape);
	}

	void Project::compute_minimum_surface_pad()
	{
		#pragma omp parallel for
		for (int i = 0; i < teeth_count; i++)
		{
			if (m_mesh_options.upper_flags[i])
				m_func.compute_minimum_surface_pad(i, false, m_mesh_options.pad_surface);

			if (m_mesh_options.lower_flags[i])
				m_func.compute_minimum_surface_pad(i, true, m_mesh_options.pad_surface);
		}

		compute_elements(FlagPadShape);
	}
	void Project::compute_manu_guide(bool flag_lower, const int start_index, const int end_index, const bool flattened)
	{
		if(!flag_lower)
			m_func.m_upper_teeth.compute_manu_guide(flag_lower,start_index,end_index, flattened);
		else
			m_func.m_lower_teeth.compute_manu_guide(flag_lower,start_index,end_index, flattened);
	}
	void Project::compute_bonding_surface(const std::string& pad, const std::string& padoutline, const std::string& padout)
	{
		m_func.m_upper_teeth.compute_bonding_surface(pad, padoutline, padout);
		// m_func.m_lower_teeth.compute_bonding_surface();
	}

	void Project::compute_teeth_distances()
	{
		m_func.compute_teeth_distances();

		compute_elements(FlagTeethDist);
	}

	/******************* Options **********************/

	void Project::setInput_point_size(double d)
	{
		m_mesh_options.input_point_size = (float)d;
	}

	void Project::setColor_point_size(double d)
	{
		m_mesh_options.color_point_size = (float)d;
	}

	void Project::setMesh_edge_width(double d)
	{
		m_mesh_options.mesh_edge_width = (float)d;
	}

	void Project::setVector_ratio(double d)
	{
		m_mesh_options.vector_ratio = (float)d;
	}

	void Project::setInside_laplacian_parameter(double d)
	{
		m_mesh_options.inside_laplacian_param = d;
	}

	void Project::setInside_smooth_range(int i)
	{
		m_mesh_options.inside_smooth_range = i;
	}

	void Project::setInside_isovalue(double d)
	{
		m_mesh_options.inside_isovalue = d;
	}

	void Project::setInside_enlarge_dist(double d)
	{
		m_mesh_options.inside_enlarge_dist = d;
	}

	void Project::setConvex_neighbor_radius(double d)
	{
		m_mesh_options.convex_neighbor_radius = d;
	}

	void Project::setConvex_smooth_range(int i)
	{
		m_mesh_options.convex_smooth_range = i;
	}

	void Project::setConvex_laplacian_thresh(double d)
	{
		m_mesh_options.convex_laplacian_param = d;
	}

	void Project::setConvex_isovalue(double d)
	{
		m_mesh_options.convex_isovalue = d;
	}

	void Project::setGeodesic_step(int i)
	{
		m_mesh_options.geodesic_step = i;
	}

	void Project::setGeodesic_isovalue(double d)
	{
		m_mesh_options.geodesic_isovalue = d;
	}

	void Project::setGeodesic_isovalue_index(double d, int index)
	{
		m_mesh_options.geodesic_isovalues[index] = d;
	}

	void Project::setUpperFlag(int ind, int status)
	{
		m_mesh_options.upper_flags[ind] = (bool)status;
	}

	void Project::setLowerFlag(int ind, int status)
	{
		m_mesh_options.lower_flags[ind] = (bool)status;
	}

	void Project::setPair_distance(double d)
	{
		for(int i = 0; i < m_mesh_options.pair_distances.size(); i++)
			m_mesh_options.pair_distances[i] = d;
	}

	void Project::setPair_tol(double d)
	{
		for(int i = 0; i < m_mesh_options.pair_tols.size(); i++)
			m_mesh_options.pair_tols[i] = d;
	}

	void Project::setPair_max_bend(int max_bend)
	{
		for(int i = 0; i < m_mesh_options.pair_max_bends.size(); i++)
			m_mesh_options.pair_max_bends[i] = max_bend;
	}

	void Project::setPair_offset(double d)
	{
		for(int i = 0; i < m_mesh_options.pair_offsets.size(); i++)
			m_mesh_options.pair_offsets[i] = d;
	}

	void Project::setPair_distance_index(int index, double d)
	{
		if(index < 0 || index >= m_mesh_options.pair_distances.size())
		{
			std::cout << "Wrong index for pair distance!" << std::endl;
			return;
		}
		m_mesh_options.pair_distances[index] = d;
	}

	void Project::setPair_tol_index(int index, double d)
	{
		if(index < 0 || index >= m_mesh_options.pair_tols.size())
		{
			std::cout << "Wrong index for pair tolerance!" << std::endl;
			return;
		}
		m_mesh_options.pair_tols[index] = d;
	}

	void Project::setPair_max_bend_index(int index, int i)
	{
		if(index < 0 || index >= m_mesh_options.pair_max_bends.size())
		{
			std::cout << "Wrong index for pair max bend!" << std::endl;
			return;
		}
		m_mesh_options.pair_max_bends[index] = i;
	}

	void Project::setPair_offset_index(int index, double d)
	{
		if(index < 0 || index >= m_mesh_options.pair_offsets.size())
		{
			std::cout << "Wrong index for pair offset!" << std::endl;
			return;
		}
		m_mesh_options.pair_offsets[index] = d;
	}

	void Project::setPair_enlarge_ratio(double d)
	{
		m_mesh_options.pair_enlarge_ratio = d;
	}

	void Project::setPair_radius(double d)
	{
		m_mesh_options.pair_radius = d;
	}

	void Project::setPair_width(double d)
	{
		m_mesh_options.pair_width = d;
	}

	void Project::setPad_height(double d)
	{
		m_mesh_options.pad_height = d;
	}

	void Project::setPad_wire_height(double d)
	{
		m_mesh_options.pad_wire_height = d;
	}

	void Project::setPad_minimum_area(double d)
	{
		m_mesh_options.pad_minimim_area = d;
	}

	void Project::setTeeth_dist_thresh(double d)
	{
		m_mesh_options.teeth_dist_thresh = d;
	}

static const uint32_t SAVE_CHECKSUM = 0x56879525;
static const float CURRENT_VERSION = 1.0f;


void Project::save(FILE* output)
{
	FileWriter file(output);

	file.write_header(SAVE_CHECKSUM, CURRENT_VERSION);

	file.write_key_value_cb("m_func", [&](FileWriter& writer) { m_func.save(writer); });
	file.write_key_value("m_select_pos", m_select_pos);
	file.write_key_value_cb("m_mesh_options", [&](FileWriter& writer) { m_mesh_options.save(writer); });
	//file.write_key_value("m_view", m_view, SceneObjectTypeCount);

	// Need to write again to write the total data size in the header.
	file.write_header(SAVE_CHECKSUM, CURRENT_VERSION);
}


void Project::load(FILE* input)
{
	FileReader reader(input);
	if (!reader.check_header(SAVE_CHECKSUM, CURRENT_VERSION)) return;

	while (!reader.empty)
	{
		reader.next_key();
		if (reader.empty)
			break;

		FileReader value = reader.get_value_reader();

		reader.read_key_value_cb("m_func", [&](FileReader& reader) { m_func.load(value); });
		reader.read_key_value("m_select_pos", m_select_pos);
		reader.read_key_value_cb("m_mesh_options", [&](FileReader& reader) { m_mesh_options.load(value); });
		//reader.read_key_value("m_view", m_view, SceneObjectTypeCount);
	}
}


#define MESH_OPTIONS_LIST() \
	MESH_OPTION(float, input_point_size) \
	MESH_OPTION(float, color_point_size) \
	MESH_OPTION(float, mesh_edge_width) \
	MESH_OPTION(float, vector_ratio) \
	MESH_OPTION(double, inside_laplacian_param) \
	MESH_OPTION(int, inside_smooth_range) \
	MESH_OPTION(double, inside_isovalue) \
	MESH_OPTION(double, inside_enlarge_dist) \
	MESH_OPTION(double, convex_neighbor_radius) \
	MESH_OPTION(int, convex_smooth_range) \
	MESH_OPTION(double, convex_laplacian_param) \
	MESH_OPTION(double, convex_isovalue) \
	MESH_OPTION(int, geodesic_step) \
	MESH_OPTION(double, geodesic_isovalue) \
	MESH_OPTION(std::vector<double>, pair_distances) \
	MESH_OPTION(double, pair_width) \
	MESH_OPTION(double, pad_height) \
	MESH_OPTION(double, pad_wire_height) \
	MESH_OPTION(std::vector<uint8_t>, upper_flags) \
	MESH_OPTION(std::vector<uint8_t>, lower_flags) \



void MeshOptions::save(FileWriter& file)
{
#undef MESH_OPTION
#define MESH_OPTION(type, name) file.write_key_value(#name, name);
	MESH_OPTIONS_LIST();

}

void MeshOptions::load(FileReader& reader)
{
	// Reset to default.
	*this = {};
	
	while (!reader.empty)
	{
		reader.next_key();

#undef MESH_OPTION
#define MESH_OPTION(type, name) else if (reader.read_key_value(#name, name)) {}

		if (false) {}
		MESH_OPTIONS_LIST();
	}
}



}
