#include "optimisation_data.hpp"



LaplacianOptimisationData::LaplacianOptimisationData(const std::vector<std::array<double, 3> > points, const double width, const double m, const double d, const bool dx, const bool dy, const bool dz):
    optimisation_points(points),
    laplacian_width(width),
    particle_mass(m),
    particle_diameter(d),
    dx(dx),
    dy(dy),
    dz(dz)
{}
