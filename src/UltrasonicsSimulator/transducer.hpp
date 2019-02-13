#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <complex>

#include <stator/symbolic/symbolic.hpp>



class Transducer {

    double PktoPkA(const double) const;
public:

    Eigen::Vector3d pos;
    
    Eigen::Vector3d director;

    double k;

    double phi;

    const double wavelength = 0.00865;

    const double piston_radius = 0.0045;

    const double p0 = 0.364;        
    
    Transducer(const Eigen::Vector3d pos, const Eigen::Vector3d dir, const double phi=0);

    Transducer(const std::array<double, 3> pos, const std::array<double, 3> dir, const double phi=0);

    std::complex<double> pressure(const Eigen::Vector3d point, const double shift=0) const;
    
    Eigen::Vector3cd nablap(const Eigen::Vector3d point, const double shift=0) const;



    // // stator
    // auto pressure_equation() const;

    // std::complex<double> stator_pressure(const Eigen::Vector3d point, const double shift=0) const;
    
    // Eigen::Vector3cd stator_nablap(const Eigen::Vector3d point, const double shift=0) const;

};
