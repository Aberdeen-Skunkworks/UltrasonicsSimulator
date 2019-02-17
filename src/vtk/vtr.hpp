#pragma once
#include <string>
#include <sstream>
#include <fstream>
#include <sys/stat.h>
#include <Eigen/Dense>

#include "../UltrasonicsSimulator/field.hpp"
#include "pvd.hpp"



namespace vtr {

    template<class T>
    void dump(const Field<T>& f, const std::string);

    
    template<class T>
    void dump_series(const std::vector<Field<T> >&, const std::string);

    
    std::string header(const std::array<size_t, 3>);

    
    std::string bottom();

    
    std::string cell_data();

    
    std::string coordinates(const std::array<size_t, 3>, const std::array<double, 3>, const std::array<double, 3>);

    
    std::string point_data(const Field<double>&);

    
    std::string point_data(const Field<Eigen::Vector3d>&);    

    
    template<class T, typename std::enable_if<std::is_base_of<Eigen::EigenBase<T>, T>::value>::type>
    std::string point_data(const Field<T>&);
    
}



#include "vtr.ipp"
