#pragma once
#include <array>
#include <Eigen/Dense>



class Particle {
    
 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    
    Particle(const Eigen::Vector3d, const double, const double);

    
    Particle(const std::array<double, 3>, const double, const double);


    const Eigen::Vector3d pos;

    
    const double diameter;

    
    const double mass;

    
    const double volume() const;

    
    const double density() const;
    
};
