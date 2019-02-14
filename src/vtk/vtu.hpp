#pragma once
#include <vector>
#include <sstream>
#include <fstream>

#include "../UltrasonicsSimulator/simulation.hpp"


namespace vtu {

    void dump(const Simulation&, const std::string);

    void dump(const std::vector<Transducer>&, const std::string);

    std::string header(const size_t);

    std::string bottom();

    std::string point_data(const std::vector<Transducer>&);
    
    std::string cell_data();

    std::string points(const std::vector<Transducer>&);
    
    std::string cells();

}
