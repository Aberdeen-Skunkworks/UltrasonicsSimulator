


{

    // Field class
    py::class_<Field<double>>(m, "Field")
	.def(py::init<const std::array<size_t, 3>,
	     const std::array<double, 3>,
	     const std::array<double, 3>>(),
	     py::arg("N"),
	     py::arg("L"),
	     py::arg("origin"))
	.def("N",
	     py::overload_cast<size_t>(&Field<double>::N, py::const_),
	     "returns the number of points in the specified dimension",
	     py::arg("dim"))
	.def("L",
	     py::overload_cast<size_t>(&Field<double>::L, py::const_),
	     "returns the Field length in the specified dimension",
	     py::arg("dim"))
	.def("N",
	     py::overload_cast<>(&Field<double>::N, py::const_),
	     "returns the number of points in all 3 dimensions of the Field")
	.def("L",
	     py::overload_cast<>(&Field<double>::L, py::const_),
	     "returns the length in all 3 dimensions of the Field")
	.def("__getitem__", [](const Field<double>& f, const std::array<size_t, 3> c) {
		return f(c);
	    });
    
}
