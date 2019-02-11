#pragma once
#include <iostream>
#include <vector>
#include <array>


template<class T>
class Field {

    const std::array<size_t, 3> _N;

    const std::array<double, 3> _L;

    T& pos(std::array<size_t, 3>);
    const T& pos(std::array<size_t, 3>) const;
    
public:

    std::vector<T> values;

    Field(const std::array<size_t, 3>, const std::array<double, 3>, const std::array<double, 3>);

    T& operator() (std::array<size_t, 3>);
    const T& operator() (std::array<size_t, 3>) const;
    
    std::array<size_t, 3> N() const;
    
    size_t N(const size_t) const;

    std::array<double, 3> L() const;
    
    double L(const size_t) const;

    const std::array<double, 3> orig;
    
};



#include "Field.ipp"
