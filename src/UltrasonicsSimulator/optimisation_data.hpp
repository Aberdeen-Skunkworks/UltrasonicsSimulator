#pragma once
#include <vector>
#include <array>



class Simulation;
    


class LaplacianOptimisationData {

public:
    
    std::vector<std::array<double, 3> > optimisation_points;

    
    double laplacian_width;

    
    double particle_mass;

    
    double particle_diameter;    

    
    bool dx;

    
    bool dy;

    
    bool dz;

    
    LaplacianOptimisationData(const std::vector<std::array<double, 3> >, const double, const double, const double, const bool, const bool, const bool);    
    
};



typedef struct {
    Simulation* sim;
    LaplacianOptimisationData* data;
} LaplacianSimCombo;
