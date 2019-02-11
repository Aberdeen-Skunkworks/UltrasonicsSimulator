#pragma once
#include <vector>
#include <complex>

#include "Transducer.hpp"
#include "Particle.hpp"



class ParticleSystem {

    const double gravity = -9.81;

    const double frequency = 40000;

    const double sound_speed_air = 346; // m/s

    const double sound_speed_particle = 2600; // m/s

    const double air_density = 1.2; // kg/m^3

public:

    ParticleSystem();

    std::vector<Transducer> transducers;
    
    void add_transducer(const Transducer);

    void add_transducers(const std::vector<Transducer>);

    void clear();

    std::complex<double> pressure_sum(const Eigen::Vector3d, const double shift=0) const;
    
    Eigen::Vector3cd nablap_sum(const Eigen::Vector3d, const double shift=0) const;

    double Gorkov_potential(const Particle&) const;
    
    void focus(const Eigen::Vector3d);
    
};
