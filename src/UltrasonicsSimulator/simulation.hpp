#pragma once
#include <vector>
#include <complex>

#include "nlopt.hpp"

#include "transducer.hpp"
#include "particle.hpp"
#include "field.hpp"
#include "optimisation_data.hpp"


typedef struct {
    Simulation* sim;
    LaplacianOptimisationData* data;
} Test;



class Simulation {

    // const double gravity = -9.81;
    Eigen::Vector3d gravity;
    
    const double frequency = 40000;

    const double sound_speed_air = 346; // m/s

    const double sound_speed_particle = 2600; // m/s

    const double air_density = 1.2; // kg/m^3

    size_t optimisation_counter=0;
    
public:

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
    
    void optimise_Gorkov_laplacian(const std::vector<std::array<double, 3> >&, const double, const double, const double, const bool x, const bool y, const bool z);

    double optimise_laplacian_function(const std::vector<double>& , std::vector<double>&, void*);
    
};

double owrapper(const std::vector<double>& x, std::vector<double>& grad, void* wrapper_func_data);
