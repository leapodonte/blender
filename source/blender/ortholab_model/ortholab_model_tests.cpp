#include "ortholab_model.h"
#include "stdio.h"
#include "stdlib.h"

#include "cgal_types.h"
// #include "db/db_services.h"
#include "vedis.h"
#include "kvstore.h"
#include "spdlog/spdlog.h"
#include <omp.h>
/***********************test kvstore********************/
int ortholab_model_test_db()
{
	// check openmp
	#pragma omp parallel
	spdlog::trace("OpenMP thread {0}, nthreads {1}", omp_get_thread_num(), omp_get_num_threads());

	// KVStore
	for (int test_iteration_i = 0; test_iteration_i < 10; test_iteration_i++)
	{
		std::string db_name = "test_output/db.db";
		OrthoLab::KVStore store(db_name);
		for (int test_iteration_j = 0; test_iteration_j < 100; test_iteration_j++)
		{
			store.store_int("int_key", 10);
			assert(store.retrieve_int("int_key") == 10);

			store.store_double("double_key", 10.0);
			assert(store.retrieve_double("double_key") == 10.0);

			std::vector<double> doubles3{ 0, 1, 2, 3, 4, 1.12, 142.244554, 3654.11254, 785.35466 };
			int rt = store.store_vector("double_list_test", doubles3);
			assert(rt == VEDIS_OK);
			std::vector<double> doubles4;
			rt = store.retrieve_vector("double_list_test", doubles4);
			assert(rt == VEDIS_OK);
			assert(doubles4.size() == doubles3.size());
			for (int i = 0; i < doubles3.size(); ++i)
				assert(doubles4[i] == doubles3[i]);

			std::vector<int> ints4{ 0, 1, 2, 3, 4, 5, 701, 457, 214, 7840 };
			rt = store.store_vector<int>("int_list_test", ints4);
			assert(rt == VEDIS_OK);
			std::vector<int> ints3;
			rt = store.retrieve_vector<int>("int_list_test", ints3);
			assert(rt == VEDIS_OK);
			assert(ints4.size() == ints3.size());
			for (int i = 0; i < ints3.size(); ++i)
				assert(ints3[i] == ints4[i]);
		}
	}
	return 0;
}

#include "project.h"
#include "function.h"
/***********************test Project********************/
int ortholab_model_test_prj()
{
	OrthoLab::Project prj("test_output");
	// prj.load_from_db();
	prj.load_upper_lower_teeth("../test_input/kerrouche_Mohamed/156502001_shell_occlusion_u.stl",
		"../test_input/kerrouche_Mohamed/156502001_shell_occlusion_l.stl");
	prj.compute_laplacian_based_inside();
	prj.setConvex_isovalue(-0.04);
	prj.compute_laplacian_based_convexity();
	prj.add_selected_point(10.8195, -45.0984, -0.856306);
	prj.add_selected_point(6.97666, -48.1888, -0.232999);
	prj.add_selected_point(2.27919, -48.0177, -0.0227968);
	prj.add_selected_point(-2.63259, -49.7204, 0.367291);
	prj.add_selected_point(-6.899, -48.8317, 0.506314);
	prj.add_selected_point(-10.7966, -45.9149, -0.807502);

	prj.compute_selected_segmentation();
	prj.validate_selected_segmentation();

	prj.setGeodesic_isovalue(0.5);
	prj.compute_geodesic_pad_outlines();
	prj.compute_tooth_split_planes();

	prj.setPair_enlarge_ratio(2);
	prj.setPair_tol(0.995);
	prj.setPair_max_bend(8);
	prj.setUpperFlag(6, 0);
	prj.setLowerFlag(6, 1);
	prj.setPair_distance(2.0);
	prj.setPair_distance_index(2, 2.2);
	prj.setPair_distance_index(1, 2.4);
	prj.setPair_distance_index(3, 2.4);
	prj.compute_paired_wires();
	prj.compute_wire_shapes();
	// pad thickness
	prj.setPad_height(0.4);
	// wire
	prj.setPad_wire_height(0.9);
	prj.compute_pad_shapes();

	prj.save_pad_statistics();
	prj.flattenWires(1);

	prj.save_to_db();
	return 0;
}
