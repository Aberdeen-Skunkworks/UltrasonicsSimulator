
{
    
    // Simulation class
    py::class_<Simulation>(m, "Simulation")
	.def(py::init<>())
	.def("add_transducer", &Simulation::add_transducer, "adds a transducer to the simulation")
	.def("Gorkov_potential_field", &Simulation::Gorkov_potential_field, "returns the Gorkov potential field",
	     py::arg("N"), py::arg("L"), py::arg("origin"), py::arg("particle mass"), py::arg("particle diameter"))
	.def("Gorkov_potential", &Simulation::Gorkov_potential)
	.def("optimise_Gorkov_laplacian", &Simulation::optimise_Gorkov_laplacian)
	.def("focus", &Simulation::focus);
}
