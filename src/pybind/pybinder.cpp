#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>

#include "../UltrasonicsSimulator/simulation.hpp"
#include "../UltrasonicsSimulator/transducer.hpp"
#include "../UltrasonicsSimulator/field.hpp"

namespace py = pybind11;


PYBIND11_MODULE(ultrasonics, m) {

#include "simulation_pybind.cpp"
#include "transducer_pybind.cpp"
#include "field_pybind.cpp"
    
}
