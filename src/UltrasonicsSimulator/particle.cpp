#include "particle.hpp"



Particle::Particle(const Eigen::Vector3d pos, const double m, const double d):
    pos(pos), diameter(d), mass(m)
{}



Particle::Particle(const std::array<double, 3> pos, const double m, const double d):
    pos({pos[0], pos[1], pos[2]}), diameter(d), mass(m)
{}



const double Particle::volume() const {

    return M_PI * std::pow(diameter, 3) / 6;
    
}



const double Particle::density() const {

    return mass / volume();
    
}
