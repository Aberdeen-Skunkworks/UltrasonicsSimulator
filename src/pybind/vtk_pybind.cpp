



{
    
    // vtk functions
    
    m.def("dump",
	  py::overload_cast<const Field<double>&, const std::string>(&vtr::dump<double>));
    m.def("dump",
	  py::overload_cast<const Field<Eigen::Vector3d>&, const std::string>(&vtr::dump<Eigen::Vector3d>),
	  "creates vtr file containing the field data with specified name");

    m.def("dump",
	  py::overload_cast<const Simulation&, const std::string>(&vtu::dump));
    m.def("dump",
	  py::overload_cast<const std::vector<Transducer>&, const std::string>(&vtu::dump));

}
