#pragma once
#include <vector>
#include <sstream>
#include <fstream>

#include "../UltrasonicsSimulator/Transducer.hpp"


namespace vtu {

    void dump(const std::vector<Transducer>&, const std::string);

    std::string header(const std::array<size_t, 3>);

    std::string bottom();

    std::string point_data(const std::vector<Transducer>&);
    
    std::string cell_data();

    std::string points(const std::vector<Transducer>&);
    
    std::string cells();

}
