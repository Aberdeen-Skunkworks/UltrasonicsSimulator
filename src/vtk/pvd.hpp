#pragma once
#include <string>
#include <sstream>



namespace pvd {

    std::string header();

    
    std::string entry(const std::string filename, const size_t entry);

    
    std::string bottom();

}
