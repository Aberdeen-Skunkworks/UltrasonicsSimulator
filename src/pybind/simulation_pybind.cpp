



    
    // Simulation class
    py::class_<Simulation>(m, "Simulation")
	.def(py::init<>())
	.def_readonly("gravity", &Simulation::gravity)
    	.def_readonly("frequency", &Simulation::frequency)
    	.def_readonly("air_density", &Simulation::air_density)
    	.def_readonly("c_air", &Simulation::sound_speed_air)
	.def_readonly("c_particle", &Simulation::sound_speed_particle)
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
	     py::arg("max_iterations") = 100,
	     py::arg("tol") = 1e-4,
	     py::arg("algorithm") = "GN_ESCH")
	.def("focus",
	     py::overload_cast<const std::array<double, 3> >(&Simulation::focus),
	     "focusses the transducers on a point in space",
	     py::arg("point"))
	.def("transducer", [](Simulation& sim, const size_t index) {
		if (sim.transducers.size() < index+1) {
		    std::cout << "Error: There aren't that many transducers!" << std::endl;
		    throw;
		}
		return sim.transducers[index];
	    })
	.def(py::pickle(
			[](const Simulation& sim) { // __getstate__
			    /* Return a tuple that fully encodes the state of the object */
			    std::vector<Eigen::Vector3d> posv;
			    std::vector<Eigen::Vector3d> dirv;
			    std::vector<double> phiv;
			    for (const Transducer t : sim.transducers) {
				posv.push_back(t.pos);
				dirv.push_back(t.director);
				phiv.push_back(t.phi);
			    }
			    return py::make_tuple(posv, dirv, phiv);
			},
			[](py::tuple t) { // __setstate__
			    if (t.size() != 3)
				throw std::runtime_error("Invalid state!");

			    /* Create a new C++ instance */
			    Simulation sim;

			    /* Assign additional transducer information */
			    std::vector<Eigen::Vector3d> posv = t[0].cast<std::vector<Eigen::Vector3d> >();
			    std::vector<Eigen::Vector3d> dirv = t[1].cast<std::vector<Eigen::Vector3d> >();
			    std::vector<double> phiv = t[2].cast<std::vector<double> >();
			    
			    for (size_t i = 0; i < posv.size(); i++) {
				sim.add_transducer(Transducer(posv[i], dirv[i], phiv[i]));
			    }

			    return sim;
			}
			));
    
