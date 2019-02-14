#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/operators.h>



#include "../UltrasonicsSimulator/simulation.hpp"
#include "../UltrasonicsSimulator/particle.hpp"
#include "../UltrasonicsSimulator/transducer.hpp"
#include "../UltrasonicsSimulator/field.hpp"
#include "../vtk/vtr.hpp"
#include "../vtk/vtu.hpp"

namespace py = pybind11;


PYBIND11_MODULE(ultrasonics, m) {

#include "simulation_pybind.cpp"
#include "transducer_pybind.cpp"
#include "field_pybind.cpp"
#include "vtk_pybind.cpp"
#include "particle_pybind.cpp"
    
}
