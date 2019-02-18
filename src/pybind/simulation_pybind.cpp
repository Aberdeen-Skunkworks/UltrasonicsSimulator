


{
    
    // Simulation class
    py::class_<Simulation>(m, "Simulation")
	.def(py::init<>())
	.def("add_transducer",
	     &Simulation::add_transducer,
	     "adds a transducer to the simulation",
	     py::arg("transducer"))
	.def("add_transducers",
	     &Simulation::add_transducers,
	     "adds a list of transducers to the simulation",
	     py::arg("transducers"))
	.def("clear",
	     &Simulation::clear,
	     "removes current transducers from the simulation")
	.def("Gorkov_potential_field",
	     &Simulation::Gorkov_potential_field,
	     "returns the Gorkov potential field",
	     py::arg("N"),
	     py::arg("L"),
	     py::arg("origin"),
	     py::arg("particle mass"),
	     py::arg("particle diameter"))
	.def("Gorkov_potential",
	     &Simulation::Gorkov_potential,
	     "returns the Gorkov potential given a particle at a position in space",
	     py::arg("p"))
	.def("optimise_Gorkov_laplacian",
	     &Simulation::optimise_Gorkov_laplacian,
	     "optimises the transducer phases to optimise the size of a laplacian at a point in space",
	     py::arg("optimisation points"),
	     py::arg("laplacian width"),
	     py::arg("particle mass"),
	     py::arg("particle diameter"),
	     py::arg("dx"),
	     py::arg("dy"),
	     py::arg("dz"),
	     py::arg("max_optimisation_time") = 1,
	     py::arg("xtol") = 1e-4)
	.def("focus",
	     py::overload_cast<const std::array<double, 3> >(&Simulation::focus),
	     "focusses the transducers on a point in space",
	     py::arg("point"))
	.def("transducer", [](Simulation& sim, const size_t index) {
		return sim.transducers[index];
	    })
	.def_readonly("gravity", &Simulation::gravity)
    	.def_readonly("frequency", &Simulation::frequency)
    	.def_readonly("air_density", &Simulation::air_density)
    	.def_readonly("c_air", &Simulation::sound_speed_air)
	.def_readonly("c_particle", &Simulation::sound_speed_particle);
    
}
