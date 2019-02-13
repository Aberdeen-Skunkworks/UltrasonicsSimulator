
{

    // vtk functions
    
    m.def("dump", py::overload_cast<const Field<double>&, const std::string>(&vtr::dump<double>));
    m.def("dump", py::overload_cast<const Field<Eigen::Vector3d>&, const std::string>(&vtr::dump<Eigen::Vector3d>), "creates vtr file containing the field data with specified name");
    

}