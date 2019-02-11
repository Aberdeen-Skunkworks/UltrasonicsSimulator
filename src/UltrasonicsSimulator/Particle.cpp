#include "Particle.hpp"



Particle::Particle(const Eigen::Vector3d pos, const double m, const double d):
    pos(pos), diameter(d), mass(m)
{}



double Particle::volume() const {

    return M_PI * std::pow(diameter, 3) / 6;
    
}



double Particle::density() const {

    return mass / volume();
    
}
