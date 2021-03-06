


{
    
    // Transducer class
    py::class_<Transducer>(m, "Transducer")
	.def(py::init<const std::array<double, 3>,
	     const std::array<double, 3>,
	     const double>(),
	     py::arg("position"),
	     py::arg("direction"),
	     py::arg("phi") = 0)
	.def_readonly("pos", &Transducer::pos)
    	.def_readonly("dir", &Transducer::director)
	.def_readonly("phi", &Transducer::phi);
    
}
