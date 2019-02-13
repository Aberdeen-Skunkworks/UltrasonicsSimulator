#include "optimisation_data.hpp"



LaplacianOptimisationData::LaplacianOptimisationData(const std::vector<std::array<double, 3> > points, const double width, const double m, const double d):
    optimisation_points(points), laplacian_width(width), particle_mass(m), particle_diameter(d)
{}
