#ifndef MY_SOLVER_H
#define MY_SOLVER_H

#include "types.h"

namespace OrthoLab {

inline void assemble_laplacian_matrix(Polyhedron& polymesh, ESTripleList& LTriplets, 
EVector& BVector, VertIntMap& vert_ind_map)
{
    for(Vertex_handle vd: vertices(polymesh))
    {
        int index_i = vert_ind_map[vd];
        double diag = 0.;

        Halfedge_vertex_circulator he_begin = vd->vertex_begin();
        do {
            Vector u = he_begin->next()->vertex()->point() - vd->point();
            Vector v = he_begin->next()->vertex()->point() - he_begin->prev()->vertex()->point();
            double cot_uv = 0.5 * (u * v) / ortho_max(1e-5, std::sqrt(CGAL::cross_product(u, v).squared_length()));
            diag += cot_uv;

            int index_j = vert_ind_map[he_begin->prev()->vertex()];
            LTriplets.emplace_back(index_i, index_j, cot_uv);
        } while(++he_begin != vd->vertex_begin());

        LTriplets.emplace_back(index_i, index_i, -diag);
        BVector[index_i] = 0.;
    }
}

inline void assemble_constraint_matrix(double lambda, ESTripleList& LTriplets, 
EVector& BVector, IntDoubleMap& vert_constraint_map, int total_vert_num)
{
    for(int index = 0; index < total_vert_num; index++)
    {
        if(vert_constraint_map.find(index) != vert_constraint_map.end())
        {
            LTriplets.emplace_back(total_vert_num + index, index, lambda);
            BVector[total_vert_num + index] = lambda * vert_constraint_map[index];
        }
        else
            BVector[total_vert_num + index] = 0.;
    }
}

inline bool solve_laplacian(ESMatrix& LMatrix, EVector& BVector, EVector& X)
{
    ESolver solver;
    ESMatrix left_hand = LMatrix.transpose() * LMatrix;
    EVector right_hand = LMatrix.transpose() * BVector;
    X = solver.compute(left_hand).solve(right_hand);

    if (solver.info() != Eigen::Success)
    {
        std::cout << "Solver failed." << std::endl;
        return false;
    }
    else {
        //std::cout << "xmin:" << X.minCoeff() << ", xmax: " << X.maxCoeff() << std::endl;
        return true;
    }
}

} // namespace OrthoLab

#endif
