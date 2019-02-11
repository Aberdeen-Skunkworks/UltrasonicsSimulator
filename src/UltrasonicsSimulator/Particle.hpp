#pragma once
#include <array>
#include <Eigen/Dense>



class Particle {

 public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    Particle(const Eigen::Vector3d, const double, const double);

    Eigen::Vector3d pos;

    double diameter;

    double mass;
    
    double volume() const;

    double density() const;
    
};
