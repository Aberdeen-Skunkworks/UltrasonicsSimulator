#pragma once
#include <vector>
#include <complex>

#include "nlopt.hpp"
#include "LBFGS.h"

#include "transducer.hpp"
#include "particle.hpp"
#include "field.hpp"
#include "optimisation_data.hpp"



class Simulation {
    
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW


    const Eigen::Vector3d gravity;

    
    const double frequency = 40000;

    
    const double sound_speed_air = 346; // m/s

    
    const double sound_speed_particle = 2600; // m/s

    
    const double air_density = 1.2; // kg/m^3


    size_t optimisation_counter = 0;

    
    Simulation();

    
    std::vector<Transducer> transducers;    

    
    void add_transducer(const Transducer);

    
    void add_transducers(const std::vector<Transducer>);

    
    void clear();

    
    std::complex<double> pressure_sum(const Eigen::Vector3d, const double shift=0) const;

    
    Eigen::Vector3cd nablap_sum(const Eigen::Vector3d, const double shift=0) const;
    

    double Gorkov_potential(const Particle&) const;

    
    void focus(const Eigen::Vector3d);

    
    void focus(const std::array<double, 3>);

    
    double laplacian_sum(const std::vector<std::array<double, 3> >&, const double, const double, const double, const bool x, const bool y, const bool z) const;

    
    Field<double> Gorkov_potential_field(const std::array<size_t, 3>&, const std::array<double, 3>&, const std::array<double, 3>&, const double, const double) const;
    
    
    void optimise_Gorkov_laplacian(const std::vector<std::array<double, 3> >&, const double, const double, const double, const bool, const bool, const bool, const double, const size_t, const double, const std::string);

    
    double optimise_laplacian_function(const std::vector<double>& , std::vector<double>&, void*);
    
};



double optimise_laplacian_function_wrapper(const std::vector<double>& x, std::vector<double>& grad, void* wrapper_func_data);



class LBFGS_laplacian_optimiser {

    Simulation* sim;
    const std::vector<std::array<double, 3> > opt_points;
    const double width;
    const double mass;
    const double diameter;
    const bool dx;
    const bool dy;
    const bool dz;
    
public:

    LBFGS_laplacian_optimiser(Simulation*, const std::vector<std::array<double, 3> >, const double, const double, const double, const bool, const bool, const bool);
    
    double operator()(const Eigen::VectorXd&, Eigen::VectorXd&);
    
};
