
{
    
    // Particle class
    py::class_<Particle>(m, "Particle")
	.def(py::init<const Eigen::Vector3d, const double, const double>())
    	.def(py::init<const std::array<double, 3>, const double, const double>());
    
}
