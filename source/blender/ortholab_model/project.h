#ifndef PROJECT_H_
#define PROJECT_H_

// std
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include <filesystem>
#include <time.h>
#include <chrono>
#include <cstdint>

#include "types.h"
#include "function.h"
#include "meshing.h"
#include "kvstore.h"

#undef min
#undef max

class FileReader; 
namespace OrthoLab {

struct MeshOptions
{
	// View
	float   input_point_size = 2.5f;
	float   color_point_size = 5.f;
	float   mesh_edge_width = 0.5f;
	float	vector_ratio = 0.1f;
	// Inside Function
	double  inside_laplacian_param = 0.01;
	int 	inside_smooth_range = 3;
	double  inside_isovalue = 0.5;
	double  inside_enlarge_dist = 0.5;
	// Convexity Function
	double 	convex_neighbor_radius = 0.02;
	int 	convex_smooth_range = 3;
	double  convex_laplacian_param = 0.5;
	double  convex_isovalue = -0.065;
    // Geodesic Distance
    int     geodesic_step = 10;
    double  geodesic_isovalue = 0.5;
	std::vector<double> geodesic_isovalues = {-1, -1,-1, -1, -1, -1};
    // Tooth flags
    std::vector<uint8_t> upper_flags = {true, true, true, true, true, true, true};
    std::vector<uint8_t> lower_flags = {true, true, true, true, true, true, true};
	// Pair distance
	std::vector<double> pair_distances = {1.5, 1.5, 1.5, 1.5, 1.5};
	std::vector<double> pair_tols = {0.98, 0.98, 0.98, 0.98, 0.98};
	std::vector<int> pair_max_bends = {5, 5, 5, 5, 5};
	std::vector<double> pair_offsets = {1.3, 1.3, 1.3, 1.3, 1.3};
	double pair_enlarge_ratio = 1.5;
	double  pair_radius = 0.3;
	double 	pair_width = 0.1; 
	double  pad_height = 0.2;
	double  pad_wire_height = 0.25;
	double  pad_minimim_area = 12.;
	// Teeth distance
	double teeth_dist_thresh = 0.5;

    // pad minimun surface
	double pad_surface = 5; // 5 mm²

	void save(FileWriter& output);
	void load(FileReader& reader);
};

class Project
{
public:	
	typedef OrthoLab::Point									Point;
	typedef OrthoLab::PointList							    PointList;
	typedef OrthoLab::Vector								Vector;	
	typedef OrthoLab::DataList								DataList;
	
private:
	
	Point       m_mouse_pos;
	Vector      m_mouse_vec;
	PointList   m_select_pos;

	Function	m_func;
	// vedis*		m_store = nullptr;
	KVStore*	m_kvstore = nullptr;
	std::string m_case_name;
public:

	Function& getFunction() { return m_func; };

	Teeth& get_upper_teeth(){ return m_func.m_upper_teeth;};
    Teeth& get_lower_teeth(){ return m_func.m_lower_teeth;};

	static MeshOptions m_mesh_options;

	void save_to_db();

	void load_from_db();

	const int teeth_count = 6;
	const int pair_count = 5;

	DataList					upper_teeth_facets; 			// vao[0] buffers[0] vao[170]
	DataList					upper_teeth_normals;  			// vao[0] buffers[1] vao[170]
	DataList					lower_teeth_facets; 			// vao[1] buffers[2] vao[171]
	DataList					lower_teeth_normals;  			// vao[1] buffers[3] vao[171]
	DataList					upper_distance_colors; 			// vao[170] buffers[271]
	DataList					lower_distance_colors; 			// vao[171] buffers[272]

	DataList					upper_inside_ray_facets; 		// vao[2] buffers[4]   
	DataList					upper_inside_ray_normals; 	    // vao[2] buffers[5]   
	DataList					upper_inside_ray_colors; 		// vao[2] buffers[6]
	DataList					upper_inside_func_colors; 	    // vao[3] buffers[7]
	DataList					upper_inside_facets; 	 		// vao[4] buffers[8]
	DataList					upper_inside_normals;  	 		// vao[4] buffers[9]
	DataList					lower_inside_ray_facets; 		// vao[5] buffers[10]  
	DataList					lower_inside_ray_normals; 	    // vao[5] buffers[11]  
	DataList					lower_inside_ray_colors; 		// vao[5] buffers[12]
	DataList					lower_inside_func_colors; 	    // vao[6] buffers[13]
	DataList					lower_inside_facets; 	 		// vao[7] buffers[14]
	DataList					lower_inside_normals;  	 		// vao[7] buffers[15]

	DataList					upper_local_convex_facets; 	    // vao[8] buffers[16]  
	DataList					upper_local_convex_normals;     // vao[8] buffers[17]  
	DataList					upper_local_convex_colors; 	    // vao[8] buffers[18]
	DataList					upper_convex_func_colors; 	    // vao[9] buffers[19]
	DataList					lower_local_convex_facets; 	    // vao[10] buffers[20]  
	DataList					lower_local_convex_normals;     // vao[10] buffers[21]  
	DataList					lower_local_convex_colors; 	    // vao[10] buffers[22]
	DataList					lower_convex_func_colors; 	    // vao[11] buffers[23]
	DataList					upper_convex_iso_edges;			// vao[12] buffers[24]
	DataList					lower_convex_iso_edges;			// vao[13] buffers[25]

	DataList					selected_points;				// vao[14] buffers[26]
	DataList					selected_seg_facets;			// vao[15] buffers[27]
	DataList					selected_seg_normals;			// vao[15] buffers[28]

    std::vector<DataList>       upper_seg_attribus;             // vao[16-21] buffers[29-40]
    std::vector<DataList>       lower_seg_attribus;             // vao[22-27] buffers[41-52]

	std::vector<DataList>       upper_geodesic_attribus;        // vao[28-33] buffers[53-58]
    std::vector<DataList>       lower_geodesic_attribus;        // vao[34-39] buffers[59-64]

	std::vector<DataList>       upper_pad_outlines_attribus;    // vao[40-45] buffers[65-76]
    std::vector<DataList>       lower_pad_outlines_attribus;    // vao[46-51] buffers[77-88]

	std::vector<DataList>		upper_tooth_bounds_attribus;				// vao[52-57] buffers[89-94]
	std::vector<DataList>		lower_tooth_bounds_attribus;				// vao[58-63] buffers[95-100]

	std::vector<DataList>		upper_tooth_cut_points_attribus;			// vao[64-69] buffers[101-106]
	std::vector<DataList>		lower_tooth_cut_points_attribus;			// vao[70-75] buffers[107-112]

	std::vector<DataList>		upper_tooth_cut_planes_attribus;			// vao[76-81] buffers[113-124]
	std::vector<DataList>		lower_tooth_cut_planes_attribus;			// vao[82-87] buffers[125-136]

	std::vector<DataList>		upper_pair_fitting_planes_attribus;			// vao[88-92] buffers[137-146]
	std::vector<DataList>		lower_pair_fitting_planes_attribus;			// vao[93-97] buffers[147-156]

	std::vector<DataList>		upper_pair_fitting_lines_attribus;			// vao[98-102] buffers[157-161]
	std::vector<DataList>		lower_pair_fitting_lines_attribus;			// vao[103-107] buffers[162-166]

	std::vector<DataList>		upper_pair_fitting_points_attribus;			// vao[108-112] buffers[167-171]
	std::vector<DataList>		lower_pair_fitting_points_attribus;			// vao[113-117] buffers[172-176]

	std::vector<DataList>		upper_pair_cut_planes_attribus;				// vao[118-122] buffers[177-186]
	std::vector<DataList>		lower_pair_cut_planes_attribus;				// vao[123-127] buffers[187-196]

	std::vector<DataList>		upper_pair_wires_attribus;					// vao[128-132] buffers[197-201]
	std::vector<DataList>		lower_pair_wires_attribus;					// vao[133-137] buffers[202-206]

	std::vector<DataList>		upper_wire_shapes_attribus;					// vao[138-142] buffers[207-216]
	std::vector<DataList>		lower_wire_shapes_attribus;					// vao[143-147] buffers[217-226]

	std::vector<DataList>		upper_wire_extrude_shapes_attribus;			// vao[148-152] buffers[227-236]
	std::vector<DataList>		lower_wire_extrude_shapes_attribus;			// vao[153-157] buffers[237-246]

	std::vector<DataList>		upper_pad_shapes_attribus;					// vao[158-163] buffers[247-258]
	std::vector<DataList>		lower_pad_shapes_attribus;					// vao[164-169] buffers[259-270]

	Project(const std::string& case_name);
	~Project();

	/******************* I/O **********************/
	bool load_upper_teeth(const std::string& filename);
	bool load_lower_teeth(const std::string& filename);
	bool load_upper_lower_teeth(const std::string& filename1, 
	const std::string& filename2);
	void import_segment(const std::string &filename, bool is_upper);
	void get_bounding_shape(double& x, double& y, double& z, double& r);
	void closeTeethMesh();
	void expandTeethMesh(double dist);
	void flattenWires(double dist);
	int set_case_name(const std::string& case_name);
	const std::string& get_case_name() const {return m_case_name;}
	KVStore& get_store() const {return *m_kvstore;}
	/******************* Export mesh **********************/
	// void exportUpperTeethMesh(QJsonObject &meshObject)
	// {
	// 	m_func.m_upper_teeth.exportTeethMesh(meshObject);
	// }

	/******************* OpenGL **********************/
	//void init_gl();
	//void initialize_buffers();
	//void init_buffer(QOpenGLVertexArrayObject& vao, QOpenGLBuffer& buffer, QOpenGLShaderProgram& program, const DataList& data, const char* attribute_name);
	//void compile_shaders();
	//void attrib_buffers(CGAL::QGLViewer* viewer);
	void compute_elements(unsigned int flags);
	void clear_elements(unsigned int flags);

	/******************* Function **********************/
	void compute_laplacian_based_inside();
	void compute_laplacian_based_convexity();
	void compute_selected_segmentation();
	void update_convexity_function_isovalue();
    void validate_selected_segmentation();
	void clear_selected_points();
    void compute_geodesic_pad_outlines();
	void compute_circular_pad_outlines();
    void recompute_geodesic_pad_outlines();
	void compute_tooth_split_planes();
	void compute_paired_wires();
	void compute_coplanar_wires();
	void compute_wire_shapes();
	void compute_pad_shapes();
	void compute_minimum_surface_pad();
	void compute_manu_guide(bool flag_lower, const int start_index, const int end_index, const bool flattened);
	void compute_bonding_surface(const std::string& pad, const std::string& padoutline, const std::string& padout);
	void compute_teeth_distances();
	/******************* test **********************/
	void save_validated_segmentation();
	void load_validated_segmentation();
	void save_pad_statistics();
	void close_mesh(const std::string &pad, bool shell_solid);
	void expand_mesh(const std::string &pad,double thickness);
	/******************* Visualization **********************/
	// Teeth mesh
	void toggle_view(SceneObjectType type);

	//void render(CGAL::QGLViewer*);
	//void render_program_point(QOpenGLVertexArrayObject& vao, bool visible, const DataList data, QColor color);
 	//void render_program_line(QOpenGLVertexArrayObject& vao, bool visible, const DataList data, QColor color);
	//void render_program_mesh(QOpenGLVertexArrayObject& vao, bool visible, const DataList data, QColor color, QColor line_color);
	//void render_program_function(QOpenGLVertexArrayObject& vao, bool visible, const DataList data, QColor line_color);

	/******************* Update **********************/
	void set_mouse_pos(double x, double y, double z, double vx, double vy, double vz);
	bool update_mouse_selection();
	void add_selected_point(const double x, const double y, const double z);
	/******************* Clear **********************/
	void clear_data();

	/******************* Options **********************/
	void setInput_point_size(double d);
	void setColor_point_size(double d);

	void setMesh_edge_width(double d);
	void setVector_ratio(double d);

	void setInside_laplacian_parameter(double d);
	void setInside_smooth_range(int i);
	void setInside_isovalue(double d);
	void setInside_enlarge_dist(double d);

	void setConvex_neighbor_radius(double d);
	void setConvex_smooth_range(int i);
	void setConvex_laplacian_thresh(double d);
	void setConvex_isovalue(double d);

    void setGeodesic_step(int i);
    void setGeodesic_isovalue(double d);
	void setGeodesic_isovalue_index(double d, int index);

    void setUpperFlag(int ind, int status);
    void setLowerFlag(int ind, int status);

	void setPair_distance(double d);
	void setPair_tol(double d);
	void setPair_max_bend(int i);
	void setPair_radius(double d);
	void setPair_width(double d);
	void setPair_offset(double d);
	void setPair_distance_index(int index, double d);
	void setPair_tol_index(int index, double d);
	void setPair_max_bend_index(int index, int i);
	void setPair_offset_index(int index, double d);

	void setPair_enlarge_ratio(double d);
	void setPad_height(double d);
	void setPad_wire_height(double d);
	void setPad_minimum_area(double d);

	void setTeeth_dist_thresh(double d);

	void save(FILE* output);
	void load(FILE* input);
};
}
#endif // _PROJECT_H_
